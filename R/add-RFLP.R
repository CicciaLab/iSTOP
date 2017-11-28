# ---- RFLP ----
#' Add RFLP enzymes to iSTOP data
#'
#' @param iSTOP A dataframe resulting from [locate_PAM].
#' @param recognizes The base recognized by the enzyme. One of "c" or "t".
#' In terms of base editing, "t" will return enzymes that only cut after editing,
#' whereas "c" will return enzymes that only cut before editing.
#' @param width An integer specifiying range of unique cutting.
#' @param enzymes A dataframe of enzymes considered. Defaults to NULL which will
#' use an enzyme dataset included in the iSTOP package.
#' @param cores Number of cores to use with [pbapply][pblapply]
#'
#' @details This function only works if changed base can be distinguished. Current strategy
#' is to mark the changed base by converting it to lower case, and constructing
#' an `enzymes` table that has each corresponding base chaned to lowercase.
#' See `enzymes` definition below for more information.
#'
#' @export
#' @importFrom tidyr gather unnest
#' @importFrom purrr pmap map_int map_chr map2 map
#' @md

add_RFLP <- function(iSTOP, recognizes = 'c', width = 150, enzymes = NULL, cores = 1) {

  assertthat::assert_that(assertthat::is.number(width), assertthat::is.number(cores), assertthat::is.string(recognizes))
  assertthat::assert_that(recognizes %in% c('c', 't'), msg = '`recognizes` argument for add_RFLP() must be one of "c", or "t".')

  message('Adding RFLP...')

  enzyme_combinations <- process_enzymes(enzymes, recognizes)

  if (cores <= 1L) cores <- NULL

  # Progress bar nonsense
  pbo <- pbapply::pboptions(type = 'timer', char = '=')
  on.exit(pbapply::pboptions(pbo), add = TRUE)

  iSTOP %>%
    split(.$chr) %>%   # split on chromosome for parallelization
    pbapply::pblapply(function(df) {

      # Make sure c is in the search
      search <- ifelse(str_detect(df$searched, '[ATGC]c'), df$searched, as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(df$searched))))
      if (any(!str_detect(search, 'c'))) warning('add_RFLP currently only supports searching of reference sequences that contain a "c" in either strand.')
      target <- str_locate(search, 'c')[,'start']
      search <- stringr::str_sub(search, start = target - width, end = target + width)

      if (recognizes == 't') {
        search <- stringr::str_replace(search, 'c', 't')
        column_name <- paste0('RFLP_T_', width)
      } else {
        column_name <- paste0('RFLP_C_', width)
      }

      df[[column_name]] <- identify_RFLP_enzymes(search, enzyme_combinations)
      return(df)
    }, cl = cores) %>%
    bind_rows()
}

# Converts all Cs to lowercase one at a time to generate all possible match patterns
# Supports ambiguities by converting compatible patterns to 'c' and excluding incombatible ambiguous matches
process_enzymes <- function(enzymes = NULL, recognizes) {

  if (is.null(enzymes)) {
    enzymes <-
      system.file('db/restriction-enzymes.csv', package = 'iSTOP') %>%
      read_csv(col_types = cols_only(enzyme  = 'c', forward = 'c', reverse = 'c', exclude = 'l')) %>%
      filter(!exclude)
  }

  switch(
    recognizes,
    c = process_enzymes_c(enzymes),
    t = process_enzymes_t(enzymes)
  )

}

# ---- Process enzymes for recognizing a particular base ----

process_enzymes_c <- function(enzymes) {
  enzymes %>%
    select(enzyme, forward, reverse) %>%
    # Be agnostic to strand orientation (i.e. keep both)
    tidyr::gather(orientation, pattern, forward, reverse) %>%
    # Change each occurence of 'C' to 'c' one at a time (since this is the intended base change)
    mutate(pattern = str_to_lower_each(pattern, 'C')) %>%
    tidyr::unnest() %>%
    # If there is an ambiguous pattern for 'c' that does not include 'T' then change
    # pattern to only match 'c' (since 'T' is what 'C' will be changed to)
    mutate(pattern = stringr::str_replace(pattern, '\\[Ac\\]|\\[cG\\]|\\[AcG\\]', 'c')) %>%
    # If 'c' is contained within an ambiguity that includes 'T' then remove pattern
    # (since we do not want to match after 'C' is chaned to 'T')
    filter(!stringr::str_detect(pattern, '\\[cT\\]|\\[AcT\\]|\\[cGT\\]')) %>%
    select(enzyme, pattern) %>%
    distinct %>%
    # Replace standard forward and reverse patterns
    left_join(enzymes, by = 'enzyme') %>%
    # Ensure that forward and reverse patterns will match even when overlapping
    mutate(
      forward = stringr::str_c('(?=', stringr::str_to_upper(forward), ')'),
      reverse = stringr::str_c('(?=', stringr::str_to_upper(reverse), ')'),
      fwd_rev = if_else(forward == reverse, forward, stringr::str_c(forward, '|', reverse))
    )
}

process_enzymes_t <- function(enzymes) {
  enzymes %>%
    select(enzyme, forward, reverse) %>%
    # Be agnostic to strand orientation (i.e. keep both)
    tidyr::gather(orientation, pattern, forward, reverse) %>%
    # Change each occurence of 'T' to 't' one at a time (since this is the intended base change)
    mutate(pattern = str_to_lower_each(pattern, 'T')) %>%
    tidyr::unnest() %>%
    # If there is an ambiguous pattern for 't' that does not include 'C' then change pattern to only match 't'
    mutate(pattern = stringr::str_replace(pattern, '\\[At\\]|\\[Gt\\]|\\[AGt\\]', 't')) %>%
    # If 't' is contained within an ambiguity that includes 'C' then remove pattern
    filter(!stringr::str_detect(pattern, '\\[Ct\\]|\\[ACt\\]|\\[CGt\\]')) %>%
    select(enzyme, pattern) %>%
    distinct %>%
    # Replace standard forward and reverse patterns
    left_join(enzymes, by = 'enzyme') %>%
    # Ensure that forward and reverse patterns will match even when overlapping
    mutate(
      forward = stringr::str_c('(?=', stringr::str_to_upper(forward), ')'),
      reverse = stringr::str_c('(?=', stringr::str_to_upper(reverse), ')'),
      fwd_rev = if_else(forward == reverse, forward, stringr::str_c(forward, '|', reverse))
    )
}

# ---- Utilities: add_RFLP ----

identify_RFLP_enzymes <- function(seqs, enzymes) {
  # Get all unique patterns
  basic <- enzymes %>% select(enzyme, fwd_rev) %>% distinct()
  purrr::map_chr(seqs, match_one_seq, basic, enzymes)
}


match_one_seq <- function(seq, basic, enzymes) {
  # The enzymes that match the base edit (a small number and faster pattern matching)
  match_base_edit <- enzymes$enzyme[which(stringr::str_detect(seq, enzymes$pattern))]

  # Get outta here if nothing matched the intended base edit
  if (length(match_base_edit) == 0) return(NA_character_)

  # Otherwise check how many times the enzyme cuts the sequence
  search <- basic[basic$enzyme %in% match_base_edit, ]
  nmatch <- stringr::str_locate_all(seq, regex(search$fwd_rev, ignore_case = T)) %>% purrr::map_int(nrow)
  uniquely_match_base_edit <- search$enzyme[which(nmatch == 1L)]

  # Get outta here if nothing uniquely matched the intended base edit
  if (length(uniquely_match_base_edit) == 0) return(NA_character_)

  # Otherwise collapse enzymes into a single string
  stringr::str_c(uniquely_match_base_edit, collapse = ' | ')
}

# ---- Misc utility functions ----

# Default value in case of NULL
#
# Replace length == 0 object (e.g. NULL) with some other default
#
# @param lhs Left hand side is the object
# @param rhs Right hand side is replacement in case of `length == 0`
#
# @example
#   object  <- 'lhs'
#   default <- 'rhs'
#   object %||% default
#
#   object <- NULL
#   object %||% default
#
#   object <- character(0)
#   object %||% default

'%||%' <- function(lhs, rhs) if (length(lhs)) lhs else rhs


# Make each occurence of a pattern in a sting lower case individually
str_to_lower_each <- function(string, pattern, locale = 'en') {
  loc <- stringr::str_locate_all(string, pattern)
  lhs <- purrr::map2(string, loc, ~stringr::str_sub(.x, start =  1L, end   = .y[, 'start'] - 1L))
  rhs <- purrr::map2(string, loc, ~stringr::str_sub(.x, end   = -1L, start = .y[, 'end']   + 1L))
  mid <- purrr::map2(string, loc, ~stringr::str_sub(.x, start = .y[, 'start'], end = .y[, 'end']) %>% stringr::str_to_lower(locale))

  purrr::pmap(list(lhs, mid, rhs), str_c)
}
