# ---- Locate iSTOP ----
#' @export

locate_iSTOP <- function(codons, genome) {

  if (nrow(codons) < 1) return(invisible(codons))

  codons <-
    codons %>%
    group_by(gene) %>%
    mutate(n_tx_in_gene = length(unique(tx)))

  no_targetable_codons <- filter(codons,  is.na(genome_coord))
  targetable_codons    <- filter(codons, !is.na(genome_coord))

  targetable_codons %>%
    group_by(gene, chr, genome_coord) %>%
    mutate(
      n_tx = n(),
      percent_tx = (n() / n_tx_in_gene) * 100
    ) %>%
    ungroup %>%
    mutate(
      searched = get_genomic_sequence(at = genome_coord, add_5prime = 150, add_3prime = 150, genome, chr, sg_strand),
      searched =  # Change target 'C' to 'c'
        str_c(
          str_sub(searched, end   = 150),                             # LHS
          str_sub(searched, start = 151, end = 151) %>% str_to_lower, # C -> c
          str_sub(searched, start = 152)                              # RHS
        ),
      #sgNG     = iSTOP(searched, PAM = 'G  ',          base_edit = 'c', spacing = c(12, 16), width = 22),
      sgNGG     = iSTOP(searched, PAM = '.GG',          base_edit = 'c', spacing = c(12, 16), width = 23),
      sgNGA     = iSTOP(searched, PAM = '.GA',          base_edit = 'c', spacing = c(12, 16), width = 23),
      #sgNGNG   = iSTOP(searched, PAM = '.G.G',         base_edit = 'c', spacing = c(12, 16), width = 24),
      sgNGCG    = iSTOP(searched, PAM = '.GCG',         base_edit = 'c', spacing = c(12, 16), width = 24),
      sgNGAG    = iSTOP(searched, PAM = '.GAG',         base_edit = 'c', spacing = c(12, 16), width = 24),
      sgNNGRRT  = iSTOP(searched, PAM = '..G[AG][AG]T', base_edit = 'c', spacing = c(12, 16), width = 26),
      sgNNNRRT  = iSTOP(searched, PAM = '...[AG][AG]T', base_edit = 'c', spacing = c(12, 16), width = 26),
      match_any =
        #!is.na(sgNG)    |
        !is.na(sgNGG)    |
        !is.na(sgNGA)    |
        #!is.na(sgNGNG)  |
        !is.na(sgNGCG)   |
        !is.na(sgNGAG)   |
        !is.na(sgNNGRRT) |
        !is.na(sgNNNRRT)
    ) %>%
    bind_rows(no_targetable_codons)
}

# Spacing is number of bases between edited base and PAM
iSTOP <- function(sequence, PAM, base_edit, spacing, width) {

  pattern <- str_c(base_edit, '.{', spacing[1], ',', spacing[2], '}', PAM)

  (str_locate(sequence, pattern)[,'end']) %>%
    str_sub(string = sequence, start = . - (width - 1), end = .)
}
