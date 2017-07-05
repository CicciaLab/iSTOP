context('RFLP')

test_that('RFLP matching considers both strands', {

  enzymes <- iSTOP:::process_enzymes(recognizes = 'c')
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
  # If match two overlapping sites then should return nothing
  expect_equal(
    iSTOP:::identify_RFLP_enzymes('GGTcCAG', enzymes),
    NA_character_
  )
  # Overlapping on opposite forward and reverse patterns should match either
  expect_equal(
    stringr::str_locate_all('GGTcCAG', regex('(?=GGTc)|(?=cCAG)', ignore_case = T)) %>% purrr::map_int(nrow),
    2
  )
})


test_that('RFLP matching works for gained cutsite with t recognition', {

  enzymes <- iSTOP:::process_enzymes(recognizes = 't')
  enzymes <- enzymes[enzymes$enzyme == 'MnlI', ]

  # The basic pattern should match if edited base is present
  expect_equal(
    iSTOP:::identify_RFLP_enzymes('CCtC', enzymes),
    'MnlI'
  )
  # But should not match if no edited base is present
  expect_equal(
    iSTOP:::identify_RFLP_enzymes('CCTC', enzymes),
    NA_character_
  )
  # If match on opposite strand then should match nothing and return an empty string
  expect_equal(
    iSTOP:::identify_RFLP_enzymes('GAGGCCtC', enzymes),
    NA_character_
  )
  # If match two overlapping sites then should return nothing
  expect_equal(
    iSTOP:::identify_RFLP_enzymes('CCTCCtC', enzymes),
    NA_character_
  )
  # If match two overlapping sites then should return nothing
  expect_equal(
    iSTOP:::identify_RFLP_enzymes('GGTtCAG', enzymes),
    NA_character_
  )
  # Overlapping on opposite forward and reverse patterns should match either
  expect_equal(
    stringr::str_locate_all('GGTtCAG', regex('(?=GGTt)|(?=tCAG)', ignore_case = T)) %>% purrr::map_int(nrow),
    2
  )
})
