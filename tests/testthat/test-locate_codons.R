context("Codons")

test_that("All codons are found", {
  hg38 <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens

  CDS <-
    readr::read_csv(
      "tx,gene,exon,chr,strand,start,end
     ENST00000379143.10,PCNA,1,chr20,-,5119578,5119798
     ENST00000379143.10,PCNA,2,chr20,-,5118769,5118866
     ENST00000379143.10,PCNA,3,chr20,-,5118610,5118677
     ENST00000379143.10,PCNA,4,chr20,-,5117470,5117664
     ENST00000379143.10,PCNA,5,chr20,-,5115449,5115572
     ENST00000379143.10,PCNA,6,chr20,-,5115283,5115362",
      col_types = "cciccii"
    )

  # Attempt to look up an AAA codon
  result <-
    iSTOP::locate_codons(
      CDS,
      hg38,
      codons = "AAA",
      positions = 1L,
      switch_strand = TRUE
    )

  # Manually get the coding sequence for this transcript
  cds_dna <-
    BSgenome::getSeq(hg38, as(select(CDS, chr, strand, start, end), "GRanges"))
  cds_dna <- stringr::str_c(as.character(cds_dna), collapse = "")

  # Manually split into codons
  codon_starts <- seq(1, stringr::str_length(cds_dna), by = 3L)
  codon_ends   <- codon_starts + 2L
  codons <- stringr::str_sub(cds_dna, start = codon_starts, end = codon_ends)

  # We know that there should be "AAA" at these positions
  expect_equal(unique(codons[c(80, 164, 217, 248)]), "AAA")

  # All manually returned AAA codon coordinates should match the aa_coord in result
  expect_equal(which(codons == "AAA"), result$aa_coord)
})
