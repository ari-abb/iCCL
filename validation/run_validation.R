#!/usr/bin/env Rscript
# iCCL validation: does predictdimension() land on a PC count that recovers
# biology as well as (or better than) a subjectively-read elbow plot?
#
# Design (deliberately non-cherry-picked): we sweep a GRID of PC dimensions,
# compute clustering-vs-reference agreement (ARI/NMI) and rare-population
# recovery (best-matching-cluster Jaccard) at each, and mark where
# predictdimension() lands. The claim is not "our number beats the elbow" but
# "predict lands on the plateau where rare populations are resolved, removing
# the subjective low-dim choice that collapses them."

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratData)
  library(ggplot2)
  library(aricode)   # ARI(), NMI()
  library(patchwork)
})

set.seed(42)

# --- load the fixed iCCL functions directly (no build/install needed) ------
# Run this script from the validation/ directory.
pkg_dir <- "../iCCL1/R"
source(file.path(pkg_dir, "predict.R"))
source(file.path(pkg_dir, "clusterloop_v1.R"))

out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- helpers ----------------------------------------------------------------

standard_pipeline <- function(obj, npcs = 50) {
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst",
                              nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = npcs, verbose = FALSE)
  obj
}

cluster_labels_at <- function(obj, d, res = 0.8) {
  obj <- FindNeighbors(obj, dims = 1:d, verbose = FALSE)
  obj <- FindClusters(obj, resolution = res, verbose = FALSE)
  as.character(Idents(obj))
}

# aricode coerces character labels with as.integer() (-> NA), so hand it
# integer-coded factors.
ari_nmi <- function(a, b) {
  a <- as.integer(factor(a)); b <- as.integer(factor(b))
  c(ARI = aricode::ARI(a, b), NMI = aricode::NMI(a, b))
}

# Best-matching-cluster Jaccard for one reference cell type: high when the type
# gets its own tight cluster, low when it is merged into a larger cluster.
type_jaccard <- function(clusters, ref, type) {
  in_type <- ref == type
  sapply_max <- 0
  for (cl in unique(clusters)) {
    in_cl <- clusters == cl
    inter <- sum(in_cl & in_type)
    union <- sum(in_cl | in_type)
    j <- if (union > 0) inter / union else 0
    if (j > sapply_max) sapply_max <- j
  }
  sapply_max
}

evaluate_dataset <- function(obj, ref_col, dataset_name,
                             dims_grid = c(3, 5, 7, 10, 15, 20, 30, 50),
                             res = 0.8, n_rare = 3) {
  message("\n==== ", dataset_name, " ====")
  obj <- standard_pipeline(obj)

  # predictdimension's objective pick (also save its diagnostic plot)
  png(file.path(out_dir, paste0(dataset_name, "_predictdimension.png")),
      width = 800, height = 600)
  pcs <- predictdimension(obj, plot = TRUE)
  dev.off()
  message(dataset_name, ": predictdimension picked ", pcs, " PCs")

  ref <- as.character(obj[[ref_col]][, 1])
  keep <- !is.na(ref) & ref != ""

  # pick the rarest annotated types (>= 10 cells) to track for recovery
  tab <- sort(table(ref[keep]))
  rare_types <- names(tab[tab >= 10])[seq_len(min(n_rare, sum(tab >= 10)))]
  message(dataset_name, ": tracking rare types -> ",
          paste(rare_types, collapse = ", "))

  npcs <- ncol(Embeddings(obj, "pca"))
  dims_grid <- sort(unique(c(dims_grid[dims_grid <= npcs], pcs)))

  rows <- list()
  for (d in dims_grid) {
    cl <- cluster_labels_at(obj, d, res)
    clk <- cl[keep]; refk <- ref[keep]
    rec <- vapply(rare_types, function(t) type_jaccard(clk, refk, t), numeric(1))
    an <- ari_nmi(clk, refk)
    rows[[as.character(d)]] <- data.frame(
      dataset    = dataset_name,
      dims       = d,
      is_predict = (d == pcs),
      n_clusters = length(unique(cl)),
      ARI        = an[["ARI"]],
      NMI        = an[["NMI"]],
      rare_mean_jaccard = mean(rec),
      t(setNames(rec, paste0("jac_", make.names(rare_types)))),
      check.names = FALSE
    )
    message(sprintf("  dims=%2d  clusters=%2d  ARI=%.3f  NMI=%.3f  rareJ=%.3f",
                    d, rows[[as.character(d)]]$n_clusters,
                    rows[[as.character(d)]]$ARI, rows[[as.character(d)]]$NMI,
                    rows[[as.character(d)]]$rare_mean_jaccard))
  }
  res_df <- do.call(rbind, rows)
  attr(res_df, "pcs") <- pcs
  res_df
}

plot_dataset <- function(res_df, dataset_name) {
  pcs <- attr(res_df, "pcs")
  long <- rbind(
    data.frame(dims = res_df$dims, metric = "ARI",              value = res_df$ARI),
    data.frame(dims = res_df$dims, metric = "NMI",              value = res_df$NMI),
    data.frame(dims = res_df$dims, metric = "rare mean Jaccard", value = res_df$rare_mean_jaccard)
  )
  p <- ggplot(long, aes(dims, value, colour = metric)) +
    geom_line() + geom_point() +
    geom_vline(xintercept = pcs, linetype = "dashed") +
    annotate("text", x = pcs, y = 0.02, label = paste0("predict = ", pcs),
             angle = 90, vjust = -0.4, hjust = 0, size = 3) +
    ylim(0, 1) +
    labs(title = paste0(dataset_name, ": clustering quality vs PC dimensionality"),
         subtitle = "dashed line = predictdimension() pick",
         x = "number of PCs used", y = "score") +
    theme_bw()
  ggsave(file.path(out_dir, paste0(dataset_name, "_dims_sweep.png")),
         p, width = 8, height = 5, dpi = 150)
  p
}

# --- PBMC3k (anchor) --------------------------------------------------------
suppressWarnings(InstallData("pbmc3k"))
pbmc <- LoadData("pbmc3k")
pbmc <- UpdateSeuratObject(pbmc)
res_pbmc <- evaluate_dataset(pbmc, "seurat_annotations", "pbmc3k")
plot_dataset(res_pbmc, "pbmc3k")
write.csv(res_pbmc, file.path(out_dir, "pbmc3k_metrics.csv"), row.names = FALSE)

# --- bmcite (harder / payoff) ----------------------------------------------
suppressWarnings(InstallData("bmcite"))
bm <- LoadData("bmcite")
bm <- UpdateSeuratObject(bm)
# downsample for a laptop-friendly dim sweep
set.seed(42)
if (ncol(bm) > 10000) bm <- bm[, sample(colnames(bm), 10000)]
ref_col_bm <- if ("celltype.l2" %in% colnames(bm[[]])) "celltype.l2" else "celltype.l1"
res_bm <- evaluate_dataset(bm, ref_col_bm, "bmcite")
plot_dataset(res_bm, "bmcite")
write.csv(res_bm, file.path(out_dir, "bmcite_metrics.csv"), row.names = FALSE)

message("\nAll done. Outputs in ", normalizePath(out_dir))
