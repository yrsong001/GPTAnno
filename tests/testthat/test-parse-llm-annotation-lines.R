test_that("LLM parser drops header lines but keeps valid predictions containing 'cell types'", {
  response <- paste(
    "Here are the cell types:",
    "1. T cell",
    "2. Mixed cell types",
    "3. NK cell",
    sep = "\n"
  )

  parsed <- GPTAnno:::.parse_llm_annotation_lines(response)

  expect_equal(parsed, c("T cell", "Mixed cell types", "NK cell"))
  expect_length(parsed, 3)
})

test_that("LLM parser handles common provider formatting without changing count", {
  response <- paste(
    "Cluster 1 - T cell",
    "Cluster 2 - B cell",
    "Cluster 3 - NK cell",
    sep = "\n"
  )

  parsed <- GPTAnno:::.parse_llm_annotation_lines(response)

  expect_equal(parsed, c("T cell", "B cell", "NK cell"))
  expect_length(parsed, 3)
})
