test_that("clean_and_match_annotation works with exact and fuzzy matches", {
  # Minimal mock mapping and ontology terms
  mapping_dict <- c("t cell" = "CL:0000084")
  ontology_terms <- c("t cell", "b cell", "macrophage")

  # Should map directly
  expect_equal(clean_and_match_annotation("t cell", mapping_dict, ontology_terms), "CL:0000084")

  # Should match ontology term
  expect_equal(clean_and_match_annotation("b cell", mapping_dict, ontology_terms), "b cell")

  # Should handle trailing 'cells'
  expect_equal(clean_and_match_annotation("macrophage cells", mapping_dict, ontology_terms), "macrophage")

  # Should return original with warning for unknown
  expect_match(clean_and_match_annotation("unknown type", mapping_dict, ontology_terms), "not in ontology")
})
