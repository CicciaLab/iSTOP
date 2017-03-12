plot_isoforms <- function(gene, coords) {
  assertthat::assert_that(
    assertthat::is.scalar(gene),
    assertthat::has_name(coords, 'gene'),
    assertthat::has_name(coords, 'tx'),
    assertthat::has_name(coords, 'chr'),
    assertthat::has_name(coords, 'strand'),
    assertthat::has_name(coords, 'start'),
    assertthat::has_name(coords, 'end')
  )
  g  <- gene
  df <- filter(coords, gene == g)
  minus_strand <- unique(df$strand) == '-'

  if (length(minus_strand) != 1) stop('plot_isoforms does not currently support multiple strand orientations. Perhaps, try splitting the isoforms with different strand orientations and plotting separately.')

  if (minus_strand) df <- mutate(df, en = end, end = start, start = en)

  df <- df %>%
    group_by(tx) %>%
    arrange(exon) %>%
    mutate(
      seg1_start = end,
      seg1_end   = end + ((lead(start) - end) * 0.5),
      seg2_start = seg1_end,
      seg2_end   = lead(start)
    )

  figure <-
    ggplot(df, aes(color = tx)) +
    geom_segment(aes(x = start,       xend = end,      y = tx, yend = tx), size = 10) +
    geom_segment(aes(x = seg1_start,  xend = seg1_end, y = tx, yend = as.integer(factor(tx)) - 0.15), color = 'black', size = 0.5, data = na.omit(df)) +
    geom_segment(aes(x = seg2_start,  xend = seg2_end, y = as.integer(factor(tx)) - 0.15, yend = tx), color = 'black', size = 0.5, data = na.omit(df)) +
    facet_wrap(~ paste0(chr, ', ', strand, ' strand'), ncol = 1, scales = 'free_x') +
    labs(y = g) +
    guides(color = 'none') +
    theme_bw() +
    theme(aspect.ratio = .30, axis.title.x = element_blank(), strip.background = element_blank())

  if (minus_strand) figure <- figure + scale_x_reverse()

  return(figure)
}


plot_isoforms_spliced <- function(gene, coords, ticks_lower = NULL, ticks_upper = NULL) {
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
    interval_left_join(intervals, by = c('start', 'end')) %>%
    select(tx, num, start = start.x, end = end.x) %>%
    mutate(
      tx     = ordered(tx, levels = txs$tx),
      tx_num = as.numeric(tx)
    )

  figure <-
    ggplot(exons_enumerated, aes(y = tx_num, fill = tx)) +
    geom_rect(aes(xmin = start, xmax = end, ymin = tx_num - 0.4, ymax = tx_num + 0.4), alpha = 0.5) +
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
    labs(y = gene, x = str_c(unique(exons$chr), ', ', unique(exons$strand), ' strand', collapse = '')) +
    theme(
      strip.background = element_blank(),
      #aspect.ratio = 500 * length(levels),
      panel.spacing.x = unit(0.001, 'npc'),
      axis.ticks = element_line(),
      axis.text.x = element_text(hjust = 1, angle = 45)
    )

  if (!is.null(ticks_lower)) {
    df_ticks_lower <- data_frame(ticks_lower) %>%
      interval_inner_join(exons_enumerated, by = c('ticks_lower' = 'start', 'ticks_lower' = 'end'))
    figure <- figure +
      geom_rect(aes(xmin = ticks_lower, xmax = ticks_lower, ymin = tx_num - 0.3, ymax = tx_num - 0.05), color = 'black', data = df_ticks_lower)
  }

  if (!is.null(ticks_upper)) {
    df_ticks_upper <- data_frame(ticks_upper) %>%
      interval_inner_join(exons_enumerated, by = c('ticks_upper' = 'start', 'ticks_upper' = 'end'))
    figure <- figure +
      geom_rect(aes(xmin = ticks_upper, xmax = ticks_upper, ymin = tx_num + 0.05, ymax = tx_num + 0.3), color = 'black', data = df_ticks_upper)
  }

  return(figure)
}
