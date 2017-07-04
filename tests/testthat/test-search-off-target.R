context('Search off target')

test_that('Off target search gets same results regardless of strand', {

  seqs <-
    c('GGGATTCATGCAACGAGCGACGG', # NGG
      'GGGATTCATGCAACGAGCGACGA', # NGA 1 mismatch
      'GGGAGGCATGCAACGAGCGACGG', # 2 mismatch
      'GGGGGGCATGCAACGAGCGACGG', # 3 mismatch
      'GGGATTCAAAAAAAAAAAAAAAA', # Many mismatch
      'GGGATTCATGCAACGAGCGANGG', # 0 mismatch 1 ambiguity
      'GGGAGTCATGCAACGAGCGANGG', # 1 mismatch 1 ambiguity
      'GGGAGGCATGCAACGAGCGANGG', # 2 mismatch 1 ambiquity
      'GGGAGGGATGCAACGAGCGANGG', # 3 mismatch 1 ambiquity
      'GGGATTCATGCAAGGAGCGACGG'  # 1 mismatch in trusted band
    )

  expect_matching <- c(1, 2, 3, 6, 7, 8)

  plus <-
    iSTOP:::search_off_target_chr(
      seqs,
      Biostrings::DNAString('AAAATGGGATTCATGCAACGAGCGACGGAAAAT'),
      orientation  = '+',
      fixed_start  = 9,
      fixed_end    = 20,
      max_mismatch = 2
    )

  minus <-
    iSTOP:::search_off_target_chr(
      seqs,
      Biostrings::reverseComplement(Biostrings::DNAString('AAAATGGGATTCATGCAACGAGCGACGGAAAAT')),
      orientation  = '-',
      fixed_start  = 9,
      fixed_end    = 20,
      max_mismatch = 2
    )

  expect_equal(
    select(plus,  -strand),
    select(minus, -strand)
  )

  expect_equal(
    plus$sequence,
    minus$sequence
  )

  expect_equal(
    plus$sequence,
    seqs[expect_matching]
  )

})

