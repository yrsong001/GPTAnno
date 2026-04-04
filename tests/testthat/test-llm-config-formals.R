test_that("public annotation entry points expose llm_config", {
  expect_true("llm_config" %in% names(formals(GPTAnno::gptcelltype)))
  expect_true("llm_config" %in% names(formals(GPTAnno::summarize_gptcelltype)))
  expect_true("llm_config" %in% names(formals(GPTAnno::gptanno)))
})
