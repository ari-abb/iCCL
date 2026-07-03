test_that("iCCL validates its arguments", {
  skip_if_not_installed("Seurat")
  obj <- make_toy_seurat(npcs = 20)
  expect_error(iCCL(obj, 2, 5), "min.dim")          # min.dim must be >= 3
  expect_error(iCCL(obj, 6, 4))                       # max.dim < min.dim
  expect_error(iCCL(obj, 5, 999), "exceeds")          # more dims than PCs
  bare <- make_bare_seurat()
  expect_error(iCCL(bare, 5, 6), "pca", ignore.case = TRUE)
})

test_that("iCCL writes one plot per dimension x resolution x reduction", {
  skip_if_not_installed("Seurat")
  obj <- make_toy_seurat(npcs = 20)

  tmp <- file.path(tempdir(), paste0("iccl_", as.integer(runif(1, 1, 1e7))))
  dir.create(tmp)
  old <- setwd(tmp)
  on.exit(setwd(old), add = TRUE)

  out <- suppressWarnings(suppressMessages(
    iCCL(obj, 5, 6, name = "unit", resolutions = c(0.5), reduction = "umap")
  ))

  expect_true(dir.exists(out))
  pngs <- list.files(out, pattern = "\\.png$")
  expect_length(pngs, 2)   # 2 dimensions x 1 resolution x 1 reduction
})
