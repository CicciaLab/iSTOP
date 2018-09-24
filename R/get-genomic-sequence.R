# ---- Genomic Sequence ----
# Get Genomic Sequence
#
# From a [BSgenome::BSgenome] object access upstream and downstream sequence
# with respect to strand at a chromosomal coordinate. A '-' character will be
# added for each position outside of sequence boundary.
#
# @param at         <int:vector> A coordinate in a chromosome. First base is at 1.
# @param add_5prime <int:vector> Number of bases to add to 5' side of `at`
#                                with respect to `strand`.
# @param add_3prime <int:vector> Number of bases to add to 3' side of `at`
#                                with respect to `strand`.
# @param genome     <BSgenome>   A BSgenome object of sequence data.
# @param chr        <chr:vector> Chromosome names. Must match names in `BSgenome`.
# @param strand     <chr:vector> Strand orientation. All strand values must
#                                be '+' or '-'.
#' @importFrom purrr map2_dbl map_chr map

get_genomic_sequence <- function(at, add_5prime, add_3prime, genome, chr, strand) {

  assertthat::assert_that(
    length(add_5prime) == length(at) | length(add_5prime) == 1L,
    length(add_3prime) == length(at) | length(add_3prime) == 1L,
    length(chr)        == length(at) | length(chr)        == 1L,
    length(strand)     == length(at) | length(strand)     == 1L,
    class(genome) %in% c('BSgenome', 'DNAStringSet'),
    all(chr %in% names(genome)),
    all(strand %in% c('+', '-'))
  )

  at         <- as.integer(round(at))
  add_5prime <- as.integer(round(add_5prime))
  add_3prime <- as.integer(round(add_3prime))

  if (length(strand) != length(at)) strand <- rep(strand, length(at))

  # Determine the boundaries of this chromosome
  max_coord <- GenomeInfoDb::seqlengths(genome)[chr]

  # Add to 'at' in 5' and 3' direction which depends on the strand
  start <- ifelse(strand == '+', at - add_5prime, at - add_3prime)
  end   <- ifelse(strand == '+', at + add_3prime, at + add_5prime)

  # Start should always be smaller of the two
  st <- pmin(start, end)
  ed <- pmax(start, end)

  # Limit ranges to bounds of chromosome
  ranges <-
    data_frame(
      seqnames  = chr,
      strand = strand,
      start  = ifelse(st < 1, 1, st),
      end    = ifelse(ed > max_coord, max_coord, ed)
    ) %>%
    as('GRanges')

  seqs <- BSgenome::getSeq(genome, ranges) %>% as.character

  # For sequences that were out of bounds, pad with '-' such that 'at' will
  # always be in the same location of the string
  seqs %>%
    str_pad_if(st < 1,         with = '-', width = 1  - st,        side = ifelse(strand == '+', 'left',  'right')) %>%
    str_pad_if(ed > max_coord, with = '-', width = ed - max_coord, side = ifelse(strand == '+', 'right', 'left'))
}


str_pad_if <- function(string, test, with, width, side = c('left', 'right')) {
  # empty string by default
  left  <- test & side == 'left'
  right <- test & side == 'right'
  # pad with string at width

  pad_left  <- purrr::map(width[left],  ~rep(with, times = .)) %>% purrr::map_chr(stringr::str_c, collapse = '')
  pad_right <- purrr::map(width[right], ~rep(with, times = .)) %>% purrr::map_chr(stringr::str_c, collapse = '')
  string[left]  <- stringr::str_c(pad_left, string[left])
  string[right] <- stringr::str_c(string[right], pad_right)

  return(string)
}
