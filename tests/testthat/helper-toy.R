# A tiny synthetic Seurat object carrying a PCA, cheap enough for unit tests.
make_toy_seurat <- function(n_cells = 120, n_genes = 400, npcs = 20, seed = 1) {
  set.seed(seed)
  counts <- matrix(rpois(n_genes * n_cells, lambda = 5), nrow = n_genes)
  rownames(counts) <- paste0("gene", seq_len(n_genes))
  colnames(counts) <- paste0("cell", seq_len(n_cells))
  suppressWarnings({
    obj <- Seurat::CreateSeuratObject(counts = counts)
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
    obj <- Seurat::FindVariableFeatures(obj, nfeatures = 150, verbose = FALSE)
    obj <- Seurat::ScaleData(obj, verbose = FALSE)
    obj <- Seurat::RunPCA(obj, npcs = npcs, verbose = FALSE)
  })
  obj
}

make_bare_seurat <- function(n_cells = 100, n_genes = 400, seed = 1) {
  set.seed(seed)
  counts <- matrix(rpois(n_genes * n_cells, lambda = 5), nrow = n_genes)
  rownames(counts) <- paste0("gene", seq_len(n_genes))
  colnames(counts) <- paste0("cell", seq_len(n_cells))
  suppressWarnings(Seurat::CreateSeuratObject(counts = counts))
}
