#' predictdimension: quantitative PCA dimensionality estimate
#'
#' Estimates how many principal components to use for downstream clustering, as
#' an objective, reproducible alternative to reading a Seurat \code{ElbowPlot()}
#' by eye. Ports the heuristic from the Harvard Chan Bioinformatics Core: take
#' the minimum of (1) the first PC whose cumulative variation exceeds 90 percent
#' while contributing less than 5 percent individually, and (2) the last PC where
#' the drop in variation relative to the next PC still exceeds 0.1 percent.
#'
#' @param SeuratObject A Seurat object that has already been through \code{RunPCA()}.
#' @param plot Logical; if \code{TRUE} (default) draw the annotated diagnostic plot.
#'
#' @return (invisibly) the estimated number of PCs, as an integer.
#' @export
#'
#' @examples
#' \dontrun{
#' pcs <- predictdimension(Macs)
#' }
predictdimension <- function(SeuratObject, plot = TRUE) {

  if (!"pca" %in% names(SeuratObject@reductions)) {
    stop("No 'pca' reduction found. Run Seurat::RunPCA() on the object first.")
  }

  # Percent of variation associated with each PC
  pct <- SeuratObject[["pca"]]@stdev / sum(SeuratObject[["pca"]]@stdev) * 100

  # Cumulative percent variation for each PC
  cumu <- cumsum(pct)

  # (1) First PC with cumulative variation > 90% and individual contribution < 5%
  co1 <- which(cumu > 90 & pct < 5)[1]

  # (2) Last PC where the PC-to-PC drop in variation is still greater than 0.1%
  co2 <- sort(which((pct[1:(length(pct) - 1)] - pct[2:length(pct)]) > 0.1),
              decreasing = TRUE)[1] + 1

  # Proceed with the more conservative (smaller) of the two estimates
  pcs <- min(co1, co2, na.rm = TRUE)

  if (isTRUE(plot)) {
    plot_df <- data.frame(pct = pct, cumu = cumu, rank = seq_along(pct))
    prediction <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = cumu, y = pct, label = rank, color = rank > pcs)) +
      ggplot2::geom_text() +
      ggplot2::geom_vline(xintercept = 90, color = "grey") +
      ggplot2::geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
      ggplot2::theme_bw() +
      ggplot2::ggtitle(paste0("predictdimension estimate: ", pcs, " PCs"))
    print(prediction)
  }

  invisible(pcs)
}
