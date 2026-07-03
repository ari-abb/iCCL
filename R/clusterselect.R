#' clusterselect: commit to a chosen dimensionality and resolution
#'
#' After inspecting the plots produced by \code{\link{iCCL}} and deciding which
#' combination of PC dimensionality and clustering resolution best fits the
#' biology, use \code{clusterselect} to apply exactly those settings to the
#' object in one call. It runs \code{FindNeighbors} -> \code{FindClusters} ->
#' \code{RunUMAP} and/or \code{RunTSNE} with the chosen values and returns the
#' updated Seurat object, so you do not have to re-type the individual Seurat
#' commands.
#'
#' @param SeuratObject A Seurat object that has already been through \code{RunPCA()}.
#' @param dims Number of principal components to use (clusters on 1:dims). Must be >= 3.
#' @param resolution The single clustering resolution you chose.
#' @param reduction Embedding(s) to compute: "umap", "tsne", or "both".
#'
#' @return The Seurat object with neighbours, clusters (\code{Idents} and
#'   \code{seurat_clusters}) and the chosen embedding(s) stored on it.
#' @export
#'
#' @examples
#' \dontrun{
#' # after comparing iCCL() output you settled on 9 PCs at resolution 0.5:
#' obj <- clusterselect(obj, dims = 9, resolution = 0.5, reduction = "both")
#' Seurat::DimPlot(obj, reduction = "umap", label = TRUE)
#' }
clusterselect <- function(SeuratObject, dims, resolution,
                          reduction = c("umap", "tsne", "both")) {

  reduction <- match.arg(reduction)
  stopifnot(length(dims) == 1, dims >= 3,
            length(resolution) == 1, resolution > 0)

  if (!"pca" %in% names(SeuratObject@reductions)) {
    stop("No 'pca' reduction found. Run Seurat::RunPCA() on the object first.")
  }
  npcs <- ncol(Seurat::Embeddings(SeuratObject, "pca"))
  if (dims > npcs) {
    stop(sprintf(paste0("dims (%d) exceeds the %d principal components available ",
                        "in the object. Lower dims or recompute PCA with more npcs."),
                 dims, npcs))
  }

  obj <- Seurat::FindNeighbors(SeuratObject, dims = 1:dims, verbose = FALSE)
  obj <- Seurat::FindClusters(obj, resolution = resolution, verbose = FALSE)
  if (reduction %in% c("umap", "both")) {
    obj <- Seurat::RunUMAP(obj, dims = 1:dims, verbose = FALSE)
  }
  if (reduction %in% c("tsne", "both")) {
    obj <- Seurat::RunTSNE(obj, dims = 1:dims)
  }

  message(sprintf("clusterselect: %d PCs, resolution %.2f, %d clusters, reduction '%s'.",
                  dims, resolution, length(levels(Seurat::Idents(obj))), reduction))
  obj
}
