context('RFLP')

test_that('RFLP matching considers both strands', {

  enzymes <- iSTOP:::process_enzymes()
  enzymes <- enzymes[enzymes$enzyme == 'MnlI', ]

  # The basic pattern should match if edited base is present
  expect_equal(
    iSTOP:::identify_RFLP_enzymes('CCTc', enzymes),
    'MnlI'
  )
  # But should not match if no edited base is present
  expect_equal(
    iSTOP:::identify_RFLP_enzymes('CCTC', enzymes),
    NA_character_
  )
  # If match on opposite strand then should match nothing and return an empty string
  expect_equal(
    iSTOP:::identify_RFLP_enzymes('GAGGCCTc', enzymes),
    NA_character_
  )
  # If match two overlapping sites then should return nothing
  expect_equal(
    iSTOP:::identify_RFLP_enzymes('CCTCCTc', enzymes),
    NA_character_
  )
})
