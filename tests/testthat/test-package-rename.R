test_that("OntoAnno exports the canonical annotation entry point", {
  expect_identical(unname(getNamespaceName(asNamespace("OntoAnno"))), "OntoAnno")
  expect_true(is.function(OntoAnno::ontoanno))
  expect_identical(
    names(formals(OntoAnno::ontoanno)),
    names(formals(OntoAnno::gptanno))
  )
})

test_that("gptanno remains a working deprecated compatibility wrapper", {
  expect_warning(
    result <- OntoAnno::gptanno(
      seurat_obj = NULL,
      resolutions = numeric(0),
      cl = NULL,
      graph = NULL
    ),
    "ontoanno"
  )
  expect_identical(result, list())
  expect_identical(
    OntoAnno::ontoanno(
      seurat_obj = NULL,
      resolutions = numeric(0),
      cl = NULL,
      graph = NULL
    ),
    list()
  )
})
