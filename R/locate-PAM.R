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
#' @param spacing The result of [PAM_spacing], which is a named list of two-element
#' numeric vectors. Each two-element vector defines a range of spacing between an
#' edited base and the PAM (i.e. the minimum and maximum allowable nucleotides between
#' the targeted base and the PAM). Each pair of ranges is checked in order with the
#' first range prioritized over the next. Defaults to `PAM_spacing()`.
#' @param PAM A named list of PAM patterns to be considered (see [PAM_pattern]). Names
#' correspond to resulting column names in the returned dataframe. Default PAM patterns
#' include `.GG`, `.GA`, `.GCG`, `.GAG`, `..G[AG][AG]T`, and `...[AG][AG]T`.
#' @param flanking An number specifying how much flanking genomic context to return
#' in the resulting dataframe. Defaults to 150.
#' @param keep_PAM Logical. Should the PAM sequence be included in the guides?
#' Defaults to `FALSE`.
#'
#' @export
#' @md

locate_PAM <- function(codons,
                       genome,
                       spacing = PAM_spacing(),
                       PAM = PAM_patterns_default(),
                       flanking = 150,
                       keep_PAM = FALSE) {

  assertthat::assert_that(
    length(flanking) == 1,
    is.numeric(flanking),
    msg = 'Please ensure that `flanking` is a single number. For more help, see ?locate_PAM'
  )

  assertthat::assert_that(
    !is.null(names(PAM)),
    !any(names(PAM) == ''),
    all(map_lgl(PAM, ~all(utils::hasName(., c('pattern', 'width'))))),
    all(map_lgl(PAM, ~is.numeric(.$width))),
    all(map_lgl(PAM, ~is.character(.$pattern))),
    msg = 'Please ensure that all `PAM` patterns are named and have a single pattern and single numeric width. For more help, see ?locate_PAM'
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
      searched =  stringr::str_c(
        stringr::str_sub(searched, end   = flanking),                                                     # LHS
        stringr::str_sub(searched, start = flanking + 1, end = flanking + 1) %>% stringr::str_to_lower(), # C -> c
        stringr::str_sub(searched, start = flanking + 2)                                                  # RHS
      ),
      no_upstream_G = !str_detect(searched, 'Gc')
    )

  # Add an sgSTOP column for each PAM
  for (i in 1:length(PAM)) {
    col_name_guide   <- names(PAM)[i]
    col_name_spacing <- names(PAM)[i] %>% stringr::str_c('_spacing')

    # Initialize columns
    sequences[[col_name_guide]]   <- NA_character_
    sequences[[col_name_spacing]] <- NA_character_
    for (j in 1:length(spacing)) {

      # Coalesce prioritizing previously detected guides
      guide_sequence <-
        coalesce(
          sequences[[col_name_guide]],
          sgSTOP(sequences$searched, pattern = PAM[[i]]$pattern, base_edit = 'c', spacing = spacing[[j]], width = 20 + PAM[[i]]$width)
        )

      guide_spacing <-
        coalesce(
          sequences[[col_name_spacing]],
          ifelse(!is.na(guide_sequence), names(spacing)[j], NA_character_)
        )

      # Optionally remove PAM sequence
      if (!keep_PAM) guide_sequence <- stringr::str_sub(guide_sequence, start = 1L, end = 20L)

      # Update
      sequences[[col_name_guide]]   <- guide_sequence
      sequences[[col_name_spacing]] <- guide_spacing
    }
  }

  # If not all of the PAM columns are NA then there is at least 1 match
  sequences$match_any <- !apply(apply(sequences[, names(PAM)], 2, is.na), 1, all)

  bind_rows(sequences, no_targetable_codons)
}

#' PAM patterns
#'
#' Specify PAM patterns for use with [locate_PAM]
#'
#' @param pattern A regular expression used to identify a PAM. See the [regex documentation][stringi::stringi-search-regex]
#' for extensive help on regular expressions.
#' @param width A single number indicating the expected width of the PAM sequence.
#'
#' @export
#' @md

PAM_pattern <- function(pattern, width) {
  assertthat::assert_that(
    is.character(pattern),
    is.numeric(width),
    msg = 'Please ensure each PAM has a pattern (e.g. ".GG") and numeric width (number of nucleotides matched by the pattern). For more help, see ?PAM_pattern'
  )
  tibble::lst(pattern, width)
}

#' @param sgNGG Pattern for NGG. Defaults to `PAM_pattern(".GG", 3)`
#' @param sgNGA Pattern for NGA. Defaults to `PAM_pattern(".GA", 3)`
#' @param sgNGCG Pattern for NGCG. Defaults to `PAM_pattern(".GCG", 4)`
#' @param sgNGAG Pattern for NGAG. Defaults to `PAM_pattern('.GAG', 4)`
#' @param sgNNGRRT Pattern for NNGRRT. Defaults to `PAM_pattern('..G[AG][AG]T', 6)`
#' @param sgNNNRRT Pattern for NNNRRT. Defaults to `PAM_pattern('...[AG][AG]T', 6)`
#' @param ... Additional PAM patterns to include in `PAM_patterns_default()`
#'
#' @rdname PAM_pattern
#' @export
#' @md

PAM_patterns_default <- function(sgNGG    = PAM_pattern('.GG', 3),
                                 sgNGA    = PAM_pattern('.GA', 3),
                                 sgNGCG   = PAM_pattern('.GCG', 4),
                                 sgNGAG   = PAM_pattern('.GAG', 4),
                                 sgNNGRRT = PAM_pattern('..G[AG][AG]T', 6),
                                 sgNNNRRT = PAM_pattern('...[AG][AG]T', 6),
                                 ...) {
  tibble::lst(sgNGG, sgNGA, sgNGCG, sgNGAG, sgNNGRRT, sgNNNRRT, ...)
}

#' PAM spacing specification
#'
#' Specify preferred PAM spacing when using [locate_PAM].
#'
#' @param optimal defaults to `c(14, 14)`
#' @param good defaults to `c(13, 15)`
#' @param ok defaults to `c(12, 16)`
#' @param ... Additional custom named ranges can be included.
#'
#' @details [locate_PAM] uses these spacing ranges to optimize PAM selection when designing
#' guide sequences. The priority for design is `optimal > good > ok > ...`. All arguments should
#' be two-element numeric vectors indicating a minimum and maximum of the range. The default
#' settings gives priority to PAMs that place the target base in the middle of the five base
#' range of high efficiency BE3-Cas9 editing.
#'
#' Note that a spacing of 14 means 14 nucleotides between the target base the first base in the PAM.
#' For example, `cNNNNNNNNNNNNNN.GG` has a spacing of 14 between the targeted `c` and the `.GG` PAM
#' pattern.
#'
#' @export
#' @md

PAM_spacing <- function(optimal = c(14, 14), good = c(13, 15), ok = c(12, 16), ...) {
  spacing <- list(optimal = optimal, good = good, ok = ok, ...)
  assertthat::assert_that(
    all(map_lgl(spacing, ~length(.) == 2)),
    all(map_lgl(spacing, is.numeric)),
    msg = 'Please ensure all PAM spacings are vectors of two numbers. For more help, see ?PAM_spacing'
  )
  return(spacing)
}

# Spacing is number of bases between edited base and PAM
sgSTOP <- function(sequence, pattern, base_edit, spacing, width) {
  regex  <- stringr::str_c(base_edit, '.{', spacing[1], ',', spacing[2], '}', pattern) # construct search pattern
  coords <- stringr::str_locate(sequence, regex)[,'end']
  return(stringr::str_sub(sequence, start = coords - (width - 1), end = coords))
}
