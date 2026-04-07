test_that("score_annotation_resolutions excludes clusters that are unknown in every run", {
  res_keep_unknown_vote <- list(
    combined_results = data.frame(
      run = c("1", "2"),
      a = c("unknown", "T cell"),
      b = c("unknown", "unknown"),
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    final_summary = data.frame(
      cluster = c("a", "b"),
      most_frequent_annotation = c("unknown", "unknown"),
      max_percentage = c(50, 100),
      other_annotations = c("t cell 50 %", ""),
      avg_distance = c(2, NA_real_),
      stringsAsFactors = FALSE
    )
  )

  res_all_scored <- list(
    combined_results = data.frame(
      run = c("1", "2"),
      a = c("T cell", "T cell"),
      b = c("B cell", "B cell"),
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    final_summary = data.frame(
      cluster = c("a", "b"),
      most_frequent_annotation = c("t cell", "b cell"),
      max_percentage = c(100, 100),
      other_annotations = c("", ""),
      avg_distance = c(1, 1),
      stringsAsFactors = FALSE
    )
  )

  scores <- score_annotation_resolutions(list(
    res_keep_unknown_vote = res_keep_unknown_vote,
    res_all_scored = res_all_scored
  ))

  keep_row <- scores[scores$resolution == "res_keep_unknown_vote", , drop = FALSE]
  scored_row <- scores[scores$resolution == "res_all_scored", , drop = FALSE]

  expect_equal(keep_row$avg_path_length, 2)
  expect_equal(keep_row$avg_max_percentage, 50)
  expect_equal(keep_row$min_max_percentage, 50)
  expect_equal(scored_row$avg_path_length, 1)
  expect_equal(scored_row$avg_max_percentage, 100)
  expect_equal(scored_row$min_max_percentage, 100)
})

test_that("score_annotation_resolutions returns NA when every cluster is unknown in every run", {
  res_failed <- list(
    combined_results = data.frame(
      run = c("1", "2"),
      a = c("unknown", "unknown"),
      b = c("unknown", "unknown"),
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    final_summary = data.frame(
      cluster = c("a", "b"),
      most_frequent_annotation = c("unknown", "unknown"),
      max_percentage = c(100, 100),
      other_annotations = c("", ""),
      avg_distance = c(NA_real_, NA_real_),
      stringsAsFactors = FALSE
    )
  )

  res_valid <- list(
    combined_results = data.frame(
      run = c("1", "2"),
      a = c("T cell", "T cell"),
      b = c("B cell", "B cell"),
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    final_summary = data.frame(
      cluster = c("a", "b"),
      most_frequent_annotation = c("t cell", "b cell"),
      max_percentage = c(100, 100),
      other_annotations = c("", ""),
      avg_distance = c(1, 2),
      stringsAsFactors = FALSE
    )
  )

  scores <- score_annotation_resolutions(list(
    res_failed = res_failed,
    res_valid = res_valid
  ))

  failed_row <- scores[scores$resolution == "res_failed", , drop = FALSE]
  valid_row <- scores[scores$resolution == "res_valid", , drop = FALSE]

  expect_true(is.na(failed_row$avg_path_length))
  expect_true(is.na(failed_row$avg_max_percentage))
  expect_true(is.na(failed_row$min_max_percentage))
  expect_true(is.na(failed_row$composite_score))
  expect_false(is.na(valid_row$composite_score))
  expect_equal(scores$resolution[1], "res_valid")
})
