#' Search for candidate off target sequences.
#'
#' Search for candidate off target sequences allowing for a fixed-band to reduce computation time
#' and a set number of allowed mismatches.
#'
#' @param guides Character vector of sequences to search for fuzzy matches in the genome.
#' @param genome A BSgenome
#' @param chromosomes Character vector of chromosome names to consider when searching the genome.
#' @param fixed_start Beginning of fixed band (no tolerated mismatches). Larger band reduces search time.
#' @param fixed_end End of fixed band (no tolerated mismatches). Larger band reduces search time.
#' @param max_mismatch Number of tolerated mismatches outside of the fixed band.
#' @param secondary_filter A regular expression passed to [stringr::str_detect] used to filter for
#' off-targets that match this pattern (e.g. 'GG$' will require the matching sequence to end with 'GG').
#'
#' @export
#' @md

search_off_target <- function(guides, genome, chromosomes,
                              fixed_start = 9,
                              fixed_end = 20,
                              max_mismatch = 2,
                              secondary_filter = NULL,
                              cores = 1) {

  message('Counting possible off targets\n',
          '  Maximum mismatches: ', max_mismatch, '\n',
          '  Fixed band (no tolerated mismatch): ',
          fixed_start, '-', fixed_end)

  # Progress bar nonsense
  if (cores <= 1L) cores <- NULL
  pbo <- pbapply::pboptions(type = 'timer', char = '=')
  on.exit(pbapply::pboptions(pbo), add = TRUE)

  # Search each chromosome
  result <-
    chromosomes %>%
    purrr::set_names(chromosomes) %>%  # give the result list names
    pbapply::pblapply(function(chrm) {
      chr <- genome[[chrm]]  # brings chromosome into memory
      bind_rows(
        search_off_target_chr(guides, chr, '+', fixed_start, fixed_end, max_mismatch),
        search_off_target_chr(guides, chr, '-', fixed_start, fixed_end, max_mismatch)
      )
    }, cl = cores) %>%
    bind_rows(.id = 'chr') %>%
    mutate(match = Biostrings::getSeq(genome, chr, start, end, strand = strand) %>% as.character()) %>%
    select(guide, match, chr, strand, start, end)

  if (!is.null(secondary_filter)) result <- filter(result, stringr::str_detect(match, secondary_filter))
  return(result)
}

search_off_target_chr <- function(guides, chromosome, orientation,
                                  fixed_start, fixed_end,
                                  max_mismatch) {

  width <- unique(stringr::str_length(guides))
  assertthat::assert_that(
    length(width) == 1,
    msg = stringr::str_c(
      'Due to current implementation limitations, all sequences considered in ',
      'an off-target search must be the same width')
  )

  if (orientation == '+') {
    search <- guides
    start  <- fixed_start
    end    <- fixed_end
  } else {
    search <- Biostrings::reverseComplement(Biostrings::DNAStringSet(guides))
    start  <- width - (fixed_end - 1)
    end    <- width - (fixed_start - 1)
  }

  search %>%
    Biostrings::PDict(
      tb.start = start, # The trusted band (tb) does not tolerate any mismatches
      tb.end   = end    # The band boundaries are inclusive
    ) %>%
    Biostrings::matchPDict(
      subject      = chromosome,
      max.mismatch = max_mismatch,
      # ambiguities in the subject (chromosome) are counted as mismatch
      fixed = 'subject'
    ) %>%
    as.list %>%
    purrr::set_names(guides) %>%
    purrr::map_df(as_data_frame, .id = 'guide') %>%
    mutate(strand = orientation)
}
