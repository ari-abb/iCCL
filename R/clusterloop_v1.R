#' iCCL: interactive Comparative Clustering Loop
#'
#' Runs Seurat clustering across a range of PCA dimensions and a set of
#' resolutions, saving one dimensionality-reduction plot per combination so the
#' results can be compared side by side. Intended to be used after
#' \code{\link{predictdimension}} has suggested a starting dimensionality.
#'
#' The embedding (UMAP and/or t-SNE) depends only on the number of dimensions,
#' so it is computed once per dimension and re-used across resolutions rather
#' than being recomputed for every resolution.
#'
#' @param SeuratObject A Seurat object that has already been through \code{RunPCA()}.
#' @param min.dim Minimum number of PCs to cluster on (must be >= 3).
#' @param max.dim Maximum number of PCs to cluster on (must be >= min.dim).
#' @param name Character label used for the output sub-directory and plot titles.
#' @param resolutions Numeric vector of clustering resolutions to sweep.
#' @param reduction Embedding(s) to compute for plotting: "umap", "tsne", or "both".
#'
#' @return (invisibly) the path to the directory containing the saved plots.
#' @export
#'
#' @examples
#' \dontrun{
#' iCCL(mySeuratObject, 20, 23, name = "Macs", reduction = "both")
#' }
iCCL <- function(SeuratObject, min.dim, max.dim, name = "iCCL",
                 resolutions = c(0.25, 0.50, 0.75, 1),
                 reduction = c("umap", "tsne", "both")) {

  reduction <- match.arg(reduction)
  stopifnot(min.dim >= 3, max.dim >= min.dim)

  reductions <- switch(reduction,
                       umap = "umap",
                       tsne = "tsne",
                       both = c("umap", "tsne"))

  # Results directory, created under the current working directory.
  out.dir <- file.path(getwd(), paste0("iCCL_results_", name))
  dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

  dim.list <- min.dim:max.dim

  for (i in dim.list) {
    neigh <- Seurat::FindNeighbors(SeuratObject, dims = 1:i, verbose = FALSE)

    # Embeddings depend only on dims -> compute once per dimension and reuse.
    embeds <- list()
    if ("umap" %in% reductions) {
      embeds$umap <- Seurat::RunUMAP(neigh, dims = 1:i, verbose = FALSE)
    }
    if ("tsne" %in% reductions) {
      embeds$tsne <- Seurat::RunTSNE(neigh, dims = 1:i)
    }

    for (r in resolutions) {
      clustered <- Seurat::FindClusters(neigh, resolution = r, verbose = FALSE)

      for (red in names(embeds)) {
        # Transfer the freshly computed clusters onto the embedded object
        # (identical cells in identical order share this lineage).
        obj <- embeds[[red]]
        obj$seurat_clusters <- clustered$seurat_clusters
        Seurat::Idents(obj) <- clustered$seurat_clusters

        nam <- paste0(name, "_Dim", i, "_Res", r, "_", red)
        message(nam)
        pl <- Seurat::DimPlot(obj, reduction = red, label = TRUE, repel = TRUE) +
          ggplot2::ggtitle(nam)
        ggplot2::ggsave(file.path(out.dir, paste0(nam, ".png")),
                        plot = pl, width = 7, height = 6, dpi = 150)
      }
    }
  }

  message("Saved plots to: ", out.dir)
  invisible(out.dir)
}
