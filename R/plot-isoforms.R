#' Plot spliced isoforms
#'
#' @param gene A single gene name
#' @param coords A dataframe of [CDS coordinates][CDS] (will be filtered for `gene`)
#' @param colors A character vector of colors to be applied to each track
#' @param ... Tracks. Each track must be a dataframe with two columns - `gene` and `genome_coord`.
#'
#' @importFrom fuzzyjoin interval_left_join interval_inner_join
#' @importFrom purrr map_lgl
#'
#' @export
#' @md

plot_spliced_isoforms <- function(gene, coords, colors = NULL, ...) {

  tracks <- list(...)
  n_tracks <- length(tracks)

  if (length(colors) != n_tracks) {
    colors <- rainbow(n_tracks)
  }

  assertthat::assert_that(
    assertthat::is.scalar(gene),
    assertthat::has_name(coords, 'gene'),
    assertthat::has_name(coords, 'tx'),
    assertthat::has_name(coords, 'chr'),
    assertthat::has_name(coords, 'strand'),
    assertthat::has_name(coords, 'start'),
    assertthat::has_name(coords, 'end'),
    all(are_colors(colors)),
    all(map_lgl(tracks, ~assertthat::has_name(., 'gene'))),
    all(map_lgl(tracks, ~assertthat::has_name(., 'genome_coord')))
  )

  g <- gene
  exons <- filter(coords, gene == g)
  txs   <- exons %>% group_by(tx) %>% summarise(size = sum(abs(start - end))) %>% arrange(size)
  minus <- exons$strand[1] == '-'

  # Determine contiguous exon intervals
  intervals <-
    with(exons, IRanges::IRanges(start, end)) %>%
    IRanges::reduce() %>%
    IRanges::as.data.frame()

  if (minus) {
    intervals <- mutate(intervals, num = n():1)
    breaks <-
      lapply(1:nrow(intervals), function(i) {
        seq(intervals$end[i], intervals$start[i], by = -500)
      }) %>%
      unlist %>%
      c(min(intervals$start))
  } else {
    intervals <- mutate(intervals, num = 1:n())
    breaks <-
      lapply(1:nrow(intervals), function(i) {
        seq(intervals$start[i], intervals$end[i], by = 500)
      }) %>%
      unlist %>%
      c(max(intervals$end))
  }

  exons_enumerated <-
    select(exons, tx, start, end) %>%
    fuzzyjoin::interval_left_join(intervals, by = c('start', 'end')) %>%
    select(tx, num, start = start.x, end = end.x) %>%
    mutate(
      tx     = ordered(tx, levels = txs$tx),
      tx_num = as.numeric(tx)
    )

  figure <-
    ggplot(exons_enumerated, aes(y = tx_num, fill = tx)) +
    geom_rect(aes(xmin = start, xmax = end, ymin = tx_num - 0.4, ymax = tx_num + 0.4), alpha = 0.2) +
    facet_grid(~num, scales = 'free_x', space = 'free_x')

  if (minus) {
    figure <- figure + scale_x_reverse(breaks = breaks, expand = c(0, 0), labels = scales::comma)
  } else {
    figure <- figure + scale_x_continuous(breaks = breaks, expand = c(0, 0), labels = scales::comma)
  }

  figure <- figure +
    scale_y_continuous(breaks = 1:nrow(txs), labels = txs$tx) +
    guides(fill = 'none') +
    theme_minimal() +
    labs(y = gene, x = stringr::str_c(unique(exons$chr), ', ', unique(exons$strand), ' strand', collapse = '')) +
    theme(
      strip.background = element_blank(),
      #aspect.ratio = 500 * length(levels),
      panel.spacing.x = unit(0.001, 'npc'),
      axis.ticks = element_line(),
      axis.text.x = element_text(hjust = 1, angle = 45)
    )

  if (n_tracks) {
    ymaxs <- 0.3 - (0:(n_tracks - 1) * (0.6 / n_tracks))

    for (i in 1:n_tracks) {
      ymax <- ymaxs[i] - 0.05
      ymin <- ymax - (0.6 / n_tracks) + 0.05
      Track <- i

      df_ticks <-
        tracks[[i]] %>%
        select(gene, genome_coord) %>%
        filter(gene == g) %>%
        fuzzyjoin::interval_inner_join(exons_enumerated, by = c('genome_coord' = 'start', 'genome_coord' = 'end')) %>%
        mutate(ymin = ymin, ymax = ymax, Track = as.character(Track))

      figure <- figure +
        geom_rect(
          aes(xmin = genome_coord, xmax = genome_coord, ymin = tx_num + ymin, ymax = tx_num + ymax, color = Track),
          data = df_ticks)
    }
    figure <- figure +
      scale_color_manual(values = colors, labels = names(tracks)) +
      guides(color = guide_legend(override.aes = list(fill = colors)))
  }

  return(figure)
}


are_colors <- function(colors) {
   map_lgl(colors, ~tryCatch(is.matrix(col2rgb(.)), error = function(e) FALSE))
}


#' Has a value
#'
#' Does the element in the vector have a value? Equivalent to `!is.na(x)`
#'
#' @param x a vector
#'
#' @export
#' @md

has <- function(x) !is.na(x)

#plot_isoforms <- function(gene, coords) {
#  assertthat::assert_that(
#    assertthat::is.scalar(gene),
#    assertthat::has_name(coords, 'gene'),
#    assertthat::has_name(coords, 'tx'),
#    assertthat::has_name(coords, 'chr'),
#    assertthat::has_name(coords, 'strand'),
#    assertthat::has_name(coords, 'start'),
#    assertthat::has_name(coords, 'end')
#  )
#  g  <- gene
#  df <- filter(coords, gene == g)
#  minus_strand <- unique(df$strand) == '-'
#
#  if (length(minus_strand) != 1) stop('plot_isoforms does not currently support multiple strand orientations. Perhaps, try splitting the isoforms with different strand orientations and plotting separately.')
#
#  if (minus_strand) df <- mutate(df, en = end, end = start, start = en)
#
#  df <- df %>%
#    group_by(tx) %>%
#    arrange(exon) %>%
#    mutate(
#      seg1_start = end,
#      seg1_end   = end + ((lead(start) - end) * 0.5),
#      seg2_start = seg1_end,
#      seg2_end   = lead(start)
#    )
#
#  figure <-
#    ggplot(df, aes(color = tx)) +
#    geom_segment(aes(x = start,       xend = end,      y = tx, yend = tx), size = 10) +
#    geom_segment(aes(x = seg1_start,  xend = seg1_end, y = tx, yend = as.integer(factor(tx)) - 0.15), color = 'black', size = 0.5, data = na.omit(df)) +
#    geom_segment(aes(x = seg2_start,  xend = seg2_end, y = as.integer(factor(tx)) - 0.15, yend = tx), color = 'black', size = 0.5, data = na.omit(df)) +
#    facet_wrap(~ paste0(chr, ', ', strand, ' strand'), ncol = 1, scales = 'free_x') +
#    labs(y = g) +
#    guides(color = 'none') +
#    theme_bw() +
#    theme(aspect.ratio = .30, axis.title.x = element_blank(), strip.background = element_blank())
#
#  if (minus_strand) figure <- figure + scale_x_reverse()
#
#  return(figure)
#}
