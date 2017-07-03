#' @export
#' @md

search_off_target <- function(seqs, genome, chromosomes,
                              fixed_start = 9,
                              fixed_end = 20,
                              max_mismatch = 2,
                              cores = 1) {

    message('Counting possible off targets\n',
            '  Maximum mismatches: ', max_mismatch, '\n',
            '  Fixed band (no tolerated mismatch): ',
            fixed_band_start, '-', fixed_band_end)

  # Progress bar nonsense
  if (cores <= 1L) cores <- NULL
  pbo <- pbapply::pboptions(type = 'timer', char = '=')
  on.exit(pbapply::pboptions(pbo), add = TRUE)

  # Search each chromosome
  chromosomes %>%
    pbapply::pblapply(function(chrm) {
      chr <- genome[[chrm]]  # brings chromosome into memory
      result <-
        bind_rows(
          search_off_target_chr(seqs, chr, '+', fixed_start, fixed_end, max_mismatch),
          search_off_target_chr(seqs, chr, '-', fixed_start, fixed_end, max_mismatch)
        )
      result$chr <- chrm
      return(result)
    }, cl = cores) %>%
    bind_rows()
}

search_off_target_chr <- function(seqs, chromosome, orientation,
                                  fixed_start, fixed_end,
                                  max_mismatch) {

  width <- unique(stringr::str_length(seqs))
  assertthat::assert_that(
    length(width) == 1,
    msg = stringr::str_c(
      'Due to current implementation limitations, all sequences considered in ',
      'an off-target search must be the same width')
  )

  if (orientation == '+') {
    search <- seqs
    start  <- fixed_start
    end    <- fixed_end
  } else {
    search <- Biostrings::reverseComplement(Biostrings::DNAStringSet(seqs))
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
    purrr::set_names(seqs) %>%
    purrr::map_df(as_data_frame, .id = 'sequence') %>%
    mutate(strand = orientation)
}
