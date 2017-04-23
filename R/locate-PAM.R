# ---- Locate iSTOP ----
#' Locate PAM near a genomic coordinate
#'
#' Uses the results of [locate_codons] to determine whether there is an available
#' PAM with desired spacing from the genomic coordinate of the targeted base.
#'
#' @param codons A dataframe resulting from [locate_codons].
#' @param genome A [BSgenome][BSgenome::BSgenome] sequence database,
#' or a [Biostrings][readDNAStringSet]. Used to extract the genomic sequence context
#' for each genomic coordinate in `codons`.
#' @param spacing A vector of two integers. The range of PAM spacing. These
#' numbers correspond to the minimum and maximum allowable nucleotides between
#' the targeted base and the PAM. Defaults to `c(12, 16)`
#' @param PAM A named vector of PAM motifs to be considered. Names correspond to
#' resulting column names in the returned dataframe. Values correspond to regular
#' expressions used to identify the PAM. Default PAM patterns include `.GG`, `.GA`
#' `.GCG`, `.GAG`, `..G[AG][AG]T`, and `...[AG][AG]T`.
#' @param PAM_widths A vector with one integer for every `PAM`. These numbers correspond
#' to the width of the PAM.
#' @param flanking An integer specifying how much flanking genomic context to return
#' in the resulting dataframe. Defaults to 150.
#'
#' @export
#' @md

locate_PAM <- function(codons,
                       genome,
                       spacing = c(12, 16),
                       PAM = c(sgNGG    = '.GG',
                               sgNGA    = '.GA',
                               sgNGCG   = '.GCG',
                               sgNGAG   = '.GAG',
                               sgNNGRRT = '..G[AG][AG]T',
                               sgNNNRRT = '...[AG][AG]T'),
                       PAM_widths = c(3, 3, 4, 4, 6, 6),
                       flanking = 150) {

  assert_that(
    length(PAM) == length(PAM_widths),
    length(spacing) == 2,
    length(flanking) == 1
  )
  if (nrow(codons) < 1) return(invisible(codons))

  codons <-
    codons %>%
    group_by(gene) %>%
    mutate(n_tx_in_gene = length(unique(tx)))

  no_targetable_codons <- filter(codons,  is.na(genome_coord))
  targetable_codons    <- filter(codons, !is.na(genome_coord))

  sequences <-
    targetable_codons %>%
    # Determine number of targeted transcripts at each coordinate
    group_by(gene, chr, genome_coord) %>%
    mutate(
      n_tx = n(),
      percent_tx = (n() / n_tx_in_gene) * 100
    ) %>%
    ungroup %>%
    mutate(
      # Get flanking genomic sequence context
      searched = get_genomic_sequence(
        at = genome_coord,
        add_5prime = flanking,
        add_3prime = flanking,
        genome,
        chr,
        sg_strand
      ),
      # make targeted base lower case
      searched =  str_c(
        str_sub(searched, end   = flanking),                                          # LHS
        str_sub(searched, start = flanking + 1, end = flanking + 1) %>% str_to_lower, # C -> c
        str_sub(searched, start = flanking + 2)                                       # RHS
      )
    )

  # Add an sgSTOP column for each PAM
  for (i in 1:length(PAM)) {
    sequences[[names(PAM)[i]]] <- sgSTOP(sequences$searched, pattern = PAM[i], base_edit = 'c', spacing = spacing, width = 20 + PAM_widths[i])
  }

  # If not all of the PAM columns are NA then there is at least 1 match
  sequences$match_any <- !apply(apply(sequences[, names(PAM)], 2, is.na), 1, all)

  bind_rows(sequences, no_targetable_codons)
}

# Spacing is number of bases between edited base and PAM
sgSTOP <- function(sequence, pattern, base_edit, spacing, width) {
  regex  <- str_c(base_edit, '.{', spacing[1], ',', spacing[2], '}', pattern) # construct search pattern
  coords <- str_locate(sequence, regex)[,'end']
  return(str_sub(sequence, start = coords - (width - 1), end = coords))
}
