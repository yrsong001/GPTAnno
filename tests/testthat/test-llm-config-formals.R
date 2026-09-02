test_that("public annotation entry points expose llm_config", {
  expect_true("llm_config" %in% names(formals(OntoAnno::gptcelltype)))
  expect_true("llm_config" %in% names(formals(OntoAnno::summarize_gptcelltype)))
  expect_true("llm_config" %in% names(formals(OntoAnno::ontoanno)))
})
