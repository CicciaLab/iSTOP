# Validate arguments for CDS function
Validate_CDS_args <- function(tx, gene, tx_cols, gene_cols, shift_start, shift_end) {

  invalid_args <- NA
  expected_tx_cols   <- c('tx', 'chr', 'strand', 'cds_start', 'cds_end', 'exon_start', 'exon_end')
  expected_gene_cols <- c('tx', 'gene')

  if (!is.string(tx)) {
    invalid_args[1] <- '`tx` must be a URL to a tab delimited trascript reference file'
  }
  if (!is.string(gene)) {
    invalid_args[2] <- '`gene` must be a URL to a tab delimited file mapping transcript IDs to gene names'
  }
  if (!all(expected_tx_cols %in% tx_cols)) {
    invalid_args[3] <- str_c('`tx_cols` requires the following missing values: "', str_c(setdiff(expected_tx_cols, tx_cols), collapse = '", "'), '"')
  }
  if (!all(expected_gene_cols %in% gene_cols)) {
    invalid_args[4] <- str_c('`gene_cols` requires the following missing values: "', str_c(setdiff(expected_gene_cols, gene_cols), collapse = '", "'), '"')
  }
  if (!is.scalar(shift_start) || !is.number(shift_start) || shift_start %% 1) {
    invalid_args[5] <- '`shift_start` must be a single integer'
  }
  if (!is.scalar(shift_end) || !is.number(shift_end) || shift_end   %% 1) {
    invalid_args[6] <- '`shift_end` must be a single integer'
  }

  invalid_args <- na.omit(invalid_args)
  if (length(invalid_args)) {
    stop('iSTOP::CDS() - Please resolve the following issues\n    ', str_c(invalid_args, collapse = '\n    '), call. = F)
  }
}

# Request user to install packages from Bioconductor
require_bioc_install <- function(require, recommend) {
  installed <- rownames(installed.packages())

  if (!all(recommend %in% installed)) {
    message(
      'Recommended package installation:\n',
      '  BiocInstaller::biocLite(c("', stringr::str_c(recommend, collapse = '", "'), '"))'
    )
  }

  if (!all(require %in% installed)) {
    stop(
      'Required package installation:\n',
      '  BiocInstaller::biocLite(c("', stringr::str_c(require, collapse = '", "'), '"))',
      call. = FALSE
    )
  }
}
