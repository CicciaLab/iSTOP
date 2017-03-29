# ---- RFLP ----
#' Add RFLP enzymes to iSTOP data
#'
#'
#' @details This function only works if changed base can be distinguished. Current strategy
#' is to mark the changed base by converting it to lower case, and constructing
#' an `enzymes` table that has each corresponding base chaned to lowercase.
#' See `enzymes` definition below for more information.
#'
#' @export
#' @importFrom readr read_csv
#' @importFrom tidyr gather unnest
#' @importFrom purrr pmap map_int map_chr map2 map
#' @md

add_RFLP <- function(iSTOP, width = 150, enzymes = NULL, cores = 1) {

  if (is.null(enzymes)) {
  enzymes <-
    system.file('db/restriction-enzymes.csv', package = 'iSTOP') %>%
    read_csv(col_types = cols_only(enzyme  = 'c', forward = 'c', reverse = 'c', exclude = 'l')) %>%
    filter(!exclude)
  }

  enzyme_cominations <-
    enzymes %>%
    select(enzyme, forward, reverse) %>%
    # Be agnostic to strand orientation (i.e. keep both)
    gather(orientation, pattern, forward, reverse) %>%
    # Change each occurence of 'C' to 'c' one at a time (since this is the intended base change)
    mutate(pattern = str_to_lower_each(pattern, 'C')) %>%
    unnest %>%
    # If there is an ambiguous pattern for 'c' that does not include 'T' then change
    # pattern to only match 'c' (since 'T' is what 'C' will be changed to)
    mutate(pattern = str_replace(pattern, '\\[Ac\\]|\\[cG\\]|\\[AcG\\]', 'c')) %>%
    # If 'c' is contained within an ambiguity that includes 'T' then remove pattern
    # (since we do not want to match after 'C' is chaned to 'T')
    filter(!str_detect(pattern, '\\[cT\\]|\\[AcT\\]|\\[cGT\\]')) %>%
    select(enzyme, pattern) %>%
    distinct

  if (cores <= 1L) cores <- NULL

  # Progress bar nonsense
  pbo <- pbapply::pboptions(type = 'timer', char = '=')
  on.exit(pbapply::pboptions(pbo), add = TRUE)

  iSTOP %>%
    split(.$chr) %>%
    pbapply::pblapply(function(df) {
      center <- ceiling(str_length(df$searched) / 2)
      search <- str_sub(df$searched, start = center - width, end = center + width)
      df[[paste0('RFLP_', width)]] <- identify_RFLP_enzymes(search, enzyme_cominations)
      return(df)
    }, cl = cores) %>%
    bind_rows
}


identify_RFLP_enzymes <- function(seqs, enzymes) {
  basic <- enzymes %>% mutate(pattern = str_to_upper(pattern)) %>% distinct
  seqs %>% map_chr(match_one_seq, basic, enzymes)
}


match_one_seq <- function(seq, basic, enzymes) {
  # How many times does the pattern match in a case insensitive way
  nmatch <- seq %>%
    str_locate_all(regex(basic$pattern, ignore_case = T)) %>%
    map_int(nrow)

  # Do the final search with just those that match only once
  search <- enzymes %>%
    filter(enzyme %in% basic$enzyme[which(nmatch == 1L)])

  if (nrow(search) == 0) return('')

  seq %>%
    map(~search$enzyme[str_detect(., search$pattern)]) %>% # subset enzymes if pattern detected in seq
    map_chr(~str_c(., collapse = ' | ') %||% NA_character_)  # collapse matches into single string
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
  lhs <- map2(string, loc, ~str_sub(.x, start =  1L, end   = .y[, 'start'] - 1L))
  rhs <- map2(string, loc, ~str_sub(.x, end   = -1L, start = .y[, 'end']   + 1L))
  mid <- map2(string, loc, ~str_sub(.x, start = .y[, 'start'], end = .y[, 'end']) %>% str_to_lower(locale))

  pmap(list(lhs, mid, rhs), str_c)
}
