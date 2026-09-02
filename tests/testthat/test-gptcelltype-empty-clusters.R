test_that("gptcelltype preserves valid predictions when some clusters have no marker genes", {
  markers <- data.frame(
    cluster = factor(
      c(rep("0", 2), rep("1", 2)),
      levels = c("0", "1", "2")
    ),
    gene = c("COL1A1", "DCN", "CD3D", "TRBC1"),
    avg_log2FC = c(2, 1.5, 2.5, 1.2),
    stringsAsFactors = FALSE
  )

  prompt_seen <- NULL
  res <- testthat::with_mocked_bindings(
    call_llm = function(prompt, ...) {
      prompt_seen <<- prompt
      paste(c("fibroblast", "t cell"), collapse = "\n")
    },
    {
      OntoAnno::gptcelltype(
        input = markers,
        tissue_name = "test tissue",
        llm_config = list(provider = "openai", model = "gpt-5", api_key = "test-key")
      )
    }
  )

  expect_equal(unname(res[c("0", "1")]), c("fibroblast", "t cell"))
  expect_equal(unname(res["2"]), "unknown")
  expect_false(grepl("\\bNA\\b", prompt_seen))
})

test_that("gptcelltype returns unknown when no usable markers remain", {
  markers <- data.frame(
    cluster = factor(character(0), levels = c("0", "1")),
    gene = character(0),
    avg_log2FC = numeric(0),
    stringsAsFactors = FALSE
  )

  llm_called <- FALSE
  res <- testthat::with_mocked_bindings(
    call_llm = function(...) {
      llm_called <<- TRUE
      "unexpected"
    },
    {
      OntoAnno::gptcelltype(
        input = markers,
        tissue_name = "test tissue",
        llm_config = list(provider = "openai", model = "gpt-5", api_key = "test-key")
      )
    }
  )

  expect_false(llm_called)
  expect_equal(unname(res), c("unknown", "unknown"), ignore_attr = TRUE)
  expect_equal(names(res), c("0", "1"))
  expect_equal(attr(res, "run_note"), "no usable marker genes after filtering")
})
