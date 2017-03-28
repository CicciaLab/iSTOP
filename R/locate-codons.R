# ---- Locate codons ----
#' Locate codons in transcripts
#'
#' Locate genome coordinates of a desired set of codons given a data.frame of
#' CDS coordinates (see [CDS]), and its corresponding
#'  [BSgenome][BSgenome::BSgenome] sequence database.
#'
#' @param cds A data.frame of CDS coordinates. See [CDS][CDS#Value] for details on required columns.
#' @param genome A [BSgenome][BSgenome::BSgenome] sequence database. CDS
#' coordinates in `cds` should correspond to this genome assembly.
#' @param codons A character vecor of codons to consider. Defaults to
#' `c('CAA', 'CAG', 'CGA', 'TGG', 'TGG')` which are all of the codons
#' with iSTOP targetable bases.
#' @param positions An integer vector of positions in `codons` to use for
#' returned coordinates. Should be the same length as `codons`. Defaults to
#' `c(1L, 1L, 1L, 2L, 3L)` which are the positions in `codons` that can be
#' targeted with iSTOP.
#' @param switch_strand A logical vector indicating whether the corresponding
#' `codon` is targeted on the opposite strand. Determines the resulting value of
#' `sg_strand`. Should be the same length as `codons`. Defaults to `c(F, F, F, T, T)`
#' @param cores Number of cores to use for parallel processing with [pblapply][pbapply::pblapply]
#'
#' @details Each transcript is processed independently based on the `tx_id`
#' column of the `cds` data.frame.
#'
#' **CDS validation** - Each CDS is validated by checking that the coordinates
#' result in a DNA sequence that begins with a start codon, ends with a stop codon,
#' has a length that is a multiple of 3 and has no internal in-frame stop codons.
#'
#' **No match?** - If a CDS passes validation, but does not have any of the
#' considered codons, the transcript will still be included in the resulting
#' data.frame, but coordinates will be missing (i.e. `NA`).
#'
#' @return A data.frame with the following columns where each row represents
#' a single targetable coordinate (codons with multiple targetable coordinates
#' will have a separate row for each):
#'
#'   - __COLUMN-NAME__ _`DATA-TYPE`_ DESCRIPTION
#'   - __tx__           _`chr`_ Transcript symbol
#'   - __gene__         _`chr`_ Gene symbol
#'   - __exon__         _`int`_ Exon rank in gene
#'   - __pep_length__   _`int`_ Lenth of peptide
#'   - __cds_length__   _`int`_ Length of CDS DNA
#'   - __chr__          _`chr`_ Chromosome
#'   - __strand__       _`chr`_ Strand (+/-)
#'   - __sg_strand__    _`chr`_ Guide strand (+/-)
#'   - __aa_target__    _`chr`_ One letter code of targeted amino acid
#'   - __codon__        _`chr`_ Three letter codon
#'   - __aa_coord__     _`int`_ Coordinate of targeted amino acid (start == 1)
#'   - __cds_coord__    _`int`_ Coordinate of targeted base in CDS (start == 1)
#'   - __genome_coord__ _`int`_ Coordinate of targeted base in genome
#'
#' @import dplyr
#' @import stringr
#' @import assertthat
#' @importFrom purrr map map2 pmap is_character is_integer is_logical
#' @importFrom tibble data_frame lst
#' @importFrom BSgenome getSeq
#' @importFrom Biostrings DNAStringSet translate
#' @importFrom pbapply pblapply
#' @export
#' @md

locate_codons <- function(cds,
                          genome,
                          codons = c('CAA', 'CAG', 'CGA', 'TGG', 'TGG'),
                          positions = c(1L, 1L, 1L, 2L, 3L),
                          switch_strand = c(F, F, F, T, T),
                          cores = getOption('mc.cores', 1L)) {

  message('Locating codon coordinates...')

  assert_that(
    cds %has_name% 'tx',  cds %has_name% 'gene',   cds %has_name% 'exon',
    cds %has_name% 'chr', cds %has_name% 'strand', cds %has_name% 'start',
    cds %has_name% 'end',
    length(codons) == length(positions),
    length(codons) == length(switch_strand),
    is_character(codons),
    is_integer(positions),
    is_logical(switch_strand),
    is_integer(cds$start),
    is_integer(cds$end),
    is.count(cores),
    all(cds$strand %in% c('+', '-')),
    class(genome) %in% c('BSgenome', 'DNAStringSet'),
    all(cds$chr %in% names(genome))
  )

  if (cores <= 1L) cores <- NULL

  # Progress bar nonsense
  pbo <- pbapply::pboptions(type = 'timer', char = '=')
  on.exit(pbapply::pboptions(pbo), add = TRUE)

  result <-
    cds %>%
    split(.$tx) %>%
    pbapply::pblapply(locate_codons_of_one_tx, genome, codons, positions, switch_strand, cl = cores) %>%
    bind_rows

  if (nrow(result) < 1) {
    warning('None of the searched CDS sequences met the criteria of:\n',
         '    1. Begin with ATG\n',
         '    2. End with TAA|TAG|TGA\n',
         '    3. Sequence length that is a multiple of 3\n',
         '    4. No internal in-frame TAA|TAG|TGA\n')
  }
  return(result)
}



locate_codons_of_one_tx <- function(cds, genome, codons, positions, switch_strand) {

  # Explicitely arrange exons by increasing exon number
  cds <- arrange(cds, exon)

  # Define constants
  opposite <- c('+' = '-', '-' = '+')

  # Construct CDS from exons
  cds_dna <-
    BSgenome::getSeq(genome, as(select(cds, chr, strand, start, end), 'GRanges')) %>%
    as.character %>%
    str_c(collapse = '')

  n_in_frame_stop <- length(which(!str_locate_all(cds_dna, 'TAA|TAG|TGA')[[1]][, 'end'] %% 3))

  # The CDS must be properly formed
  if (!str_detect(cds_dna, '^ATG') |           # Must begin with start codon
      !str_detect(cds_dna, 'TAA$|TAG$|TGA$') | # Must end with stop codon
      nchar(cds_dna) %% 3 |                    # Sequence must be multiple of 3
      n_in_frame_stop != 1) {                  # Must have only one in-frame stop codon
    return(NULL) # will NOT be counted as a transcript
  }

  # Indices of features at each position in CDS
  index_genome  <- with(cds, get_genomic_index(strand, start, end))
  index_chr     <- with(cds, rep(chr,     times = abs(end - start) + 1))
  index_strand  <- with(cds, rep(strand,  times = abs(end - start) + 1))
  index_exon    <- with(cds, rep(exon,    times = abs(end - start) + 1))


  # Find codon targets
  cds_coord <-
    tibble::lst(codons, positions, switch_strand) %>%
    pmap(function(codons, positions, switch_strand) {
      cds_position_of_first_base_in_codon <-
        str_locate_all(cds_dna, codons) %>%
        map(~.[.[, 'end'] %% 3 == 0, 'start']) %>% # extracts start position of in-frame codons
        unlist %>%
        unname

      cds_coord <- cds_position_of_first_base_in_codon + positions - 1L
      data_frame(cds_coord) %>% mutate(codon = codons, codon_position = positions, switch_strand = switch_strand)
    }) %>%
    bind_rows

  # Construct a basic results data frame
  result <- data_frame(
    tx           = unique(cds$tx),
    gene         = unique(cds$gene),
    cds_length   = nchar(cds_dna),
    pep_length   = as.integer(cds_length / 3)
  )

  # Exit if no codons were found
  if (nrow(cds_coord) < 1L) {
    return(result) # Will be counted as a valid transcript with no matches
  }

  # Using genome coordinates of desired base in codon matches to extract nearby sequence
  cds_coord %>%
    mutate(
      tx           = result$tx,
      gene         = result$gene,
      exon         = index_exon[cds_coord],
      pep_length   = result$pep_length,
      cds_length   = result$cds_length,
      chr          = index_chr[cds_coord],    # chr is exon dependent
      strand       = index_strand[cds_coord], # strand is exon dependent
      sg_strand    = ifelse(switch_strand, opposite[strand], strand),
      aa_target    = codon %>% Biostrings::DNAStringSet() %>% Biostrings::translate() %>% as.character,
      genome_coord = index_genome[cds_coord],
      aa_coord     = as.integer(ceiling(cds_coord / 3))
    ) %>%
    select(tx:aa_target, codon, aa_coord, cds_coord, genome_coord) %>%
    arrange(cds_coord)
}


# Given vectors of strand, start, and end, generate a vector of integers between each
# start and end point with respect to strand orientation
get_genomic_index <- function(strand, start, end) {
  tibble::lst(strand, start, end) %>%
    pmap(function(strand, start, end) switch(strand, '+' = start:end, '-' = end:start)) %>%
    unlist
}
