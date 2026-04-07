test_that("summarize_gptcelltype excludes failed all-unknown runs and records simple run_summary notes", {
  markers <- data.frame(
    cluster = c("0", "0", "1", "1"),
    gene = c("COL1A1", "DCN", "CD3D", "TRBC1"),
    avg_log2FC = c(2, 1.5, 2.5, 1.2),
    stringsAsFactors = FALSE
  )

  call_idx <- 0L
  result <- NULL

  expect_message(
    {
      result <- testthat::with_mocked_bindings(
        gptcelltype = function(...) {
          call_idx <<- call_idx + 1L
          responses <- list(
            structure(
              c("0" = "unknown", "1" = "unknown"),
              run_note = "batch 1 response line count mismatch"
            ),
            c("0" = "fibroblast", "1" = "t cell"),
            c("0" = "fibroblast", "1" = "b cell")
          )
          responses[[call_idx]]
        },
        {
          GPTAnno::summarize_gptcelltype(markers, n_runs = 3)
        }
      )
    },
    "Excluding failed run 1: batch 1 response line count mismatch",
    fixed = TRUE
  )

  expect_equal(call_idx, 3L)
  expect_equal(result$combined_results$run, c("2", "3"))
  expect_equal(result$run_summary$run, "1")
  expect_equal(result$run_summary$note, "batch 1 response line count mismatch")

  cluster_0 <- result$final_summary[result$final_summary$cluster == "0", , drop = FALSE]
  cluster_1 <- result$final_summary[result$final_summary$cluster == "1", , drop = FALSE]

  expect_equal(cluster_0$most_frequent_annotation, "fibroblast")
  expect_equal(cluster_0$max_percentage, 100)
  expect_equal(cluster_1$max_percentage, 50)
})
