test_that("predictdimension returns a single positive value", {
  skip_if_not_installed("Seurat")
  obj <- make_toy_seurat()
  pcs <- predictdimension(obj, plot = FALSE)
  expect_length(pcs, 1)
  expect_true(is.numeric(pcs))
  expect_gte(pcs, 1)
  expect_lte(pcs, ncol(Seurat::Embeddings(obj, "pca")))
})

test_that("predictdimension errors when PCA has not been run", {
  skip_if_not_installed("Seurat")
  bare <- make_bare_seurat()
  expect_error(predictdimension(bare, plot = FALSE), "pca", ignore.case = TRUE)
})
