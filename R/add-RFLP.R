# ---- RFLP ----
#' Add RFLP enzymes to iSTOP data
#'
#' @param iSTOP A dataframe resulting from [locate_PAM].
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

add_RFLP <- function(iSTOP, width = 150, enzymes = NULL, cores = 1) {

  message('Adding RFLP...')

  enzyme_combinations <- process_enzymes(enzymes)

  if (cores <= 1L) cores <- NULL

  # Progress bar nonsense
  pbo <- pbapply::pboptions(type = 'timer', char = '=')
  on.exit(pbapply::pboptions(pbo), add = TRUE)

  iSTOP %>%
    split(.$chr) %>%
    pbapply::pblapply(function(df) {
      center <- ceiling(str_length(df$searched) / 2)
      search <- str_sub(df$searched, start = center - width, end = center + width)
      df[[paste0('RFLP_', width)]] <- identify_RFLP_enzymes(search, enzyme_combinations)
      return(df)
    }, cl = cores) %>%
    bind_rows
}

# Converts all Cs to lowercase one at a time to generate all possible match patterns
# Supports ambiguities by converting compatible patterns to 'c' and excluding incombatible ambiguous matches
process_enzymes <- function(enzymes = NULL) {

  if (is.null(enzymes)) {
    enzymes <-
      system.file('db/restriction-enzymes.csv', package = 'iSTOP') %>%
      read_csv(col_types = cols_only(enzyme  = 'c', forward = 'c', reverse = 'c', exclude = 'l')) %>%
      filter(!exclude)
  }

  enzymes %>%
    select(enzyme, forward, reverse) %>%
    # Be agnostic to strand orientation (i.e. keep both)
    tidyr::gather(orientation, pattern, forward, reverse) %>%
    # Change each occurence of 'C' to 'c' one at a time (since this is the intended base change)
    mutate(pattern = str_to_lower_each(pattern, 'C')) %>%
    tidyr::unnest() %>%
    # If there is an ambiguous pattern for 'c' that does not include 'T' then change
    # pattern to only match 'c' (since 'T' is what 'C' will be changed to)
    mutate(pattern = str_replace(pattern, '\\[Ac\\]|\\[cG\\]|\\[AcG\\]', 'c')) %>%
    # If 'c' is contained within an ambiguity that includes 'T' then remove pattern
    # (since we do not want to match after 'C' is chaned to 'T')
    filter(!str_detect(pattern, '\\[cT\\]|\\[AcT\\]|\\[cGT\\]')) %>%
    select(enzyme, pattern) %>%
    distinct %>%
    # Replace standard forward and reverse patterns
    left_join(enzymes, by = 'enzyme') %>%
    # Ensure that forward and reverse patterns will match even when overlapping
    mutate(
      forward = str_c('(?=', str_to_upper(forward), ')'),
      reverse = str_c('(?=', str_to_upper(reverse), ')'),
      fwd_rev = str_c(forward, '|', reverse)
    )
}

identify_RFLP_enzymes <- function(seqs, enzymes) {
  # Get all unique patterns
  basic <- enzymes %>% select(enzyme, fwd_rev) %>% distinct
  purrr::map_chr(seqs, match_one_seq, basic, enzymes)
}


match_one_seq <- function(seq, basic, enzymes) {
  # The enzymes that match the base edit (a small number and faster pattern matching)
  match_base_edit <- enzymes$enzyme[which(str_detect(seq, enzymes$pattern))]

  # Get outta here if nothing matched the intended base edit
  if (length(match_base_edit) == 0) return(NA_character_)

  # Otherwise check how many times the enzyme cuts the sequence
  search <- basic[basic$enzyme %in% match_base_edit, ]
  nmatch <- str_locate_all(seq, regex(search$fwd_rev, ignore_case = T)) %>% purrr::map_int(nrow)
  uniquely_match_base_edit <- search$enzyme[which(nmatch == 1L)]

  # Get outta here if nothing uniquely matched the intended base edit
  if (length(uniquely_match_base_edit) == 0) return(NA_character_)

  # Otherwise collapse enzymes into a single string
  str_c(uniquely_match_base_edit, collapse = ' | ')
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
  loc <- str_locate_all(string, pattern)
  lhs <- purrr::map2(string, loc, ~str_sub(.x, start =  1L, end   = .y[, 'start'] - 1L))
  rhs <- purrr::map2(string, loc, ~str_sub(.x, end   = -1L, start = .y[, 'end']   + 1L))
  mid <- purrr::map2(string, loc, ~str_sub(.x, start = .y[, 'start'], end = .y[, 'end']) %>% str_to_lower(locale))

  purrr::pmap(list(lhs, mid, rhs), str_c)
}
