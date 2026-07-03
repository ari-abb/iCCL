test_that("clusterselect validates its arguments", {
  skip_if_not_installed("Seurat")
  obj <- make_toy_seurat(npcs = 20)
  expect_error(clusterselect(obj, dims = 2, resolution = 0.5), "dims")     # dims < 3
  expect_error(clusterselect(obj, dims = 999, resolution = 0.5), "exceeds") # too many PCs
  bare <- make_bare_seurat()
  expect_error(clusterselect(bare, dims = 5, resolution = 0.5), "pca", ignore.case = TRUE)
})

test_that("clusterselect stores clusters and the chosen embeddings", {
  skip_if_not_installed("Seurat")
  obj <- make_toy_seurat(npcs = 20)
  out <- suppressWarnings(suppressMessages(
    clusterselect(obj, dims = 6, resolution = 0.5, reduction = "both")
  ))
  expect_s4_class(out, "Seurat")
  expect_true("seurat_clusters" %in% colnames(out[[]]))
  expect_true(all(c("umap", "tsne") %in% names(out@reductions)))
})
