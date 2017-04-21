essential_gene_library <- function() {
  CDS   <- read_csv(system.file('db/CDS-human.csv', package = 'iSTOP'), col_types = 'cciccii')
  genes <- read_csv(system.file('db/essential-genes.csv', package = 'iSTOP'), col_types = cols(gene = col_character()))
  genome <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
  genes$gene[which(!(genes$gene %in% CDS$gene))]
  codons <-
    filter(CDS, gene %in% genes$gene) %>%
    locate_codons(genome)

  Strict <-
    locate_PAM(codons, genome, spacing = c(13, 15), PAM = c(sgNGG_13_15 = '.GG'), PAM_widths = 3) %>%
    filter(match_any, NMD_pred, !stringr::str_detect(searched, 'Gc'), percent_tx >= 66)

  Semi_strict <-
    locate_PAM(codons, genome, spacing = c(12, 16), PAM = c(sgNGG_12_16 = '.GG'), PAM_widths = 3) %>%
    filter(match_any, NMD_pred, !stringr::str_detect(searched, 'Gc'), percent_tx >= 66)

  Relaxed <-
    locate_PAM(codons, genome, spacing = c(12, 16), PAM = c(sgNGG_12_16 = '.GG'), PAM_widths = 3) %>%
    filter(match_any, NMD_pred, !stringr::str_detect(searched, 'Gc'), percent_tx >= 50)

  Loose <-
    locate_PAM(codons, genome, spacing = c(12, 16), PAM = c(sgNGG_12_16 = '.GG'), PAM_widths = 3) %>%
    filter(match_any, NMD_pred, !stringr::str_detect(searched, 'Gc'))

  Very_loose <-
    locate_PAM(codons, genome, spacing = c(12, 16), PAM = c(sgNGG_12_16 = '.GG'), PAM_widths = 3) %>%
    filter(match_any, NMD_pred)

  Extremely_loose <-
    locate_PAM(codons, genome, spacing = c(12, 16), PAM = c(sgNGG_12_16 = '.GG'), PAM_widths = 3) %>%
    filter(match_any)
}

Strict_summary <-
  Strict %>%
  group_by(gene) %>%
  summarise(n = n(), max_tx = max(percent_tx), min_tx = min(percent_tx)) %>%
  filter(n >= 6)

Semi_strict_summary <-
  Semi_strict %>%
  group_by(gene) %>%
  summarise(n = n(), max_tx = max(percent_tx), min_tx = min(percent_tx)) %>%
  filter(!(gene %in% Strict_summary$gene), n >= 6)

Relaxed_summary <-
  Relaxed %>%
  group_by(gene) %>%
  summarise(n = n(), max_tx = max(percent_tx), min_tx = min(percent_tx)) %>%
  filter(
    !(gene %in% Strict_summary$gene),
    !(gene %in% Semi_strict_summary$gene),
    n >= 6
  )

Loose_summary <-
  Loose %>%
  group_by(gene) %>%
  summarise(n = n(), max_tx = max(percent_tx), min_tx = min(percent_tx)) %>%
  filter(
    !(gene %in% Strict_summary$gene),
    !(gene %in% Semi_strict_summary$gene),
    !(gene %in% Relaxed_summary$gene),
    n >= 6
  )

Very_loose_summary <-
  Very_loose %>%
  group_by(gene) %>%
  summarise(n = n(), max_tx = max(percent_tx), min_tx = min(percent_tx)) %>%
  filter(
    !(gene %in% Strict_summary$gene),
    !(gene %in% Semi_strict_summary$gene),
    !(gene %in% Relaxed_summary$gene),
    !(gene %in% Loose_summary$gene),
    n >= 6
  )

Extremely_loose_summary <-
  Extremely_loose %>%
  group_by(gene) %>%
  summarise(n = n(), max_tx = max(percent_tx), min_tx = min(percent_tx)) %>%
  filter(
    !(gene %in% Strict_summary$gene),
    !(gene %in% Semi_strict_summary$gene),
    !(gene %in% Relaxed_summary$gene),
    !(gene %in% Loose_summary$gene),
    !(gene %in% Very_loose_summary$gene),
    n >= 6
  )

Less_than_6 <-
  Extremely_loose %>%
  group_by(gene) %>%
  summarise(n = n(), max_tx = max(percent_tx), min_tx = min(percent_tx)) %>%
  filter(
    !(gene %in% Strict_summary$gene),
    !(gene %in% Semi_strict_summary$gene),
    !(gene %in% Relaxed_summary$gene),
    !(gene %in% Loose_summary$gene),
    !(gene %in% Very_loose_summary$gene),
    n < 6
  )

table(Less_than_6$n)

genes$gene[!(genes$gene %in% Extremely_loose$gene)]

genes$gene[!(genes$gene %in% codons$gene)]

length(unique(genes$gene)) - 1
length(which(summary$n >= 6))


abline(v = 6, col = 'red', lty = 'dotted')


