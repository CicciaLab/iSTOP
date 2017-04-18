# ---- Locate iSTOP ----
#' @export

locate_PAM <- function(codons, genome) {

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
      #sgNG     = PAM(searched, pattern = 'G  ',          base_edit = 'c', spacing = c(12, 16), width = 22),
      sgNGG     = PAM(searched, pattern = '.GG',          base_edit = 'c', spacing = c(12, 16), width = 23),
      sgNGA     = PAM(searched, pattern = '.GA',          base_edit = 'c', spacing = c(12, 16), width = 23),
      #sgNGNG   = PAM(searched, pattern = '.G.G',         base_edit = 'c', spacing = c(12, 16), width = 24),
      sgNGCG    = PAM(searched, pattern = '.GCG',         base_edit = 'c', spacing = c(12, 16), width = 24),
      sgNGAG    = PAM(searched, pattern = '.GAG',         base_edit = 'c', spacing = c(12, 16), width = 24),
      sgNNGRRT  = PAM(searched, pattern = '..G[AG][AG]T', base_edit = 'c', spacing = c(12, 16), width = 26),
      sgNNNRRT  = PAM(searched, pattern = '...[AG][AG]T', base_edit = 'c', spacing = c(12, 16), width = 26),
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
PAM <- function(sequence, pattern, base_edit, spacing, width) {

  pattern <- str_c(base_edit, '.{', spacing[1], ',', spacing[2], '}', pattern)

  (str_locate(sequence, pattern)[,'end']) %>%
    str_sub(string = sequence, start = . - (width - 1), end = .)
}
