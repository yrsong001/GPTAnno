test_that("clean_and_match_annotation works with exact and fuzzy matches", {

  # Should map directly
  expect_equal(clean_and_match_annotation("t cell"), "t cell")

  # Should match ontology term
  expect_equal(clean_and_match_annotation("b cells"), "b cell")

  # Should handle trailing 'cells'
  expect_equal(clean_and_match_annotation("macrophage cells"), "macrophage")

  # Should return original with warning for unknown
  expect_match(clean_and_match_annotation("unknown"), "unknown")
})
