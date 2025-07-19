#' Preprocess a Seurat Object for Downstream Analysis
#'
#' This function normalizes, identifies variable features, scales, and performs PCA on a Seurat object.
#' Optionally saves the processed object to disk.
#'
#' @param seurat_obj A Seurat object.
#' @param assay Character. The assay to use (default: "RNA").
#' @param nfeatures Integer. Number of highly variable features to keep (default: 3000).
#' @param scale_factor Numeric. Scale factor for normalization (default: 10000).
#' @param npcs Integer. Number of principal components to compute (default: 30).
#' @param save_path Character. Path to save the processed Seurat object as RDS (default: "preprocessed_seurat.rds").
#'
#' @return The processed Seurat object.
#' @importFrom Seurat DefaultAssay NormalizeData FindVariableFeatures ScaleData RunPCA VariableFeatures
#' @export
#' @examples
#' # seurat_obj <- preprocess_seurat_object(seurat_obj)
preprocess_seurat_object <- function(seurat_obj,
                                     assay = "RNA",
                                     nfeatures = 3000,
                                     scale_factor = 10000,
                                     npcs = 30,
                                     save_path = "preprocessed_seurat.rds") {
  Seurat::DefaultAssay(seurat_obj) <- assay
  seurat_obj <- Seurat::NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = scale_factor)
  seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = nfeatures)
  seurat_obj <- Seurat::ScaleData(seurat_obj, features = Seurat::VariableFeatures(seurat_obj))
  seurat_obj <- Seurat::RunPCA(seurat_obj, features = Seurat::VariableFeatures(seurat_obj), npcs = npcs)
  saveRDS(seurat_obj, file = save_path)
  return(seurat_obj)
}

#' Run Multi-Resolution Clustering and Marker Detection
#'
#' This function runs graph-based clustering at multiple resolutions and finds marker genes for each, saving results to disk.
#'
#' @param seurat_obj A Seurat object.
#' @param resolutions Numeric vector of resolution parameters (e.g., c(0.1, 0.2, 0.3)).
#' @param result_dir Character. Directory to save marker RDS files (default: "results/markers").
#' @param dims Numeric vector. Which dimensions to use for clustering (default: 1:30).
#' @param assay Character. Assay to use, or NULL for current (default: NULL).
#' @param group.by Character. Metadata column to use for group assignment (default: "seurat_clusters").
#' @param reduction Character. Reduction to use for graph construction (default: "pca").
#' @param use_existing_neighbors Logical. Use an existing neighbor graph if present (default: TRUE).
#'
#' @return A list with: updated Seurat object, marker lists, and cluster assignments for each resolution.
#' @importFrom Seurat DefaultAssay Assays FindNeighbors FindClusters VariableFeatures Idents FindAllMarkers
#' @export
#' @examples
#' # result <- run_multi_resolution_clustering(seurat_obj, c(0.1, 0.2, 0.3))
run_multi_resolution_clustering <- function(seurat_obj,
                                            resolutions,
                                            result_dir = "output/markers",
                                            dims = 1:30,
                                            assay = NULL,
                                            group.by = "seurat_clusters",
                                            reduction = "pca",
                                            use_existing_neighbors = TRUE) {
  if (!dir.exists(result_dir)) {
    dir.create(result_dir, recursive = TRUE)
  }
  if (!is.null(assay)) {
    if (!(assay %in% Seurat::Assays(seurat_obj))) {
      stop(paste("Assay", assay, "not found. Available assays:", paste(Seurat::Assays(seurat_obj), collapse = ", ")))
    }
    Seurat::DefaultAssay(seurat_obj) <- assay
  }
  graph.name <- NULL
  if (use_existing_neighbors && length(names(seurat_obj@graphs)) > 0) {
    graph.name <- names(seurat_obj@graphs)[1]
    message("Using existing graph for clustering: ", graph.name)
  } else {
    message("Running FindNeighbors to generate neighbor graph...")
    seurat_obj <- Seurat::FindNeighbors(seurat_obj, reduction = reduction, dims = dims, assay = assay, verbose = TRUE)
    graph.name <- names(seurat_obj@graphs)[1]
  }
  all_markers_list <- list()
  cluster_assignments <- list()
  for (res in resolutions) {
    cat("Running clustering at resolution:", res, "\n")
    seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = res, graph.name = graph.name)
    cluster_col <- paste0("cluster_res.", res)
    seurat_obj@meta.data[[cluster_col]] <- seurat_obj@meta.data$seurat_clusters
    cluster_assignments[[as.character(res)]] <- seurat_obj@meta.data[[cluster_col]]
    Seurat::Idents(seurat_obj) <- seurat_obj@meta.data[[cluster_col]]
    markers <- Seurat::FindAllMarkers(seurat_obj, assay = assay, group.by = group.by)
    all_markers_list[[paste0("res_", res)]] <- markers
    saveRDS(markers, file = file.path(result_dir, paste0("markers_res_", res, ".rds")))
  }
  return(list(
    seurat_obj = seurat_obj,
    markers = all_markers_list,
    clusters = cluster_assignments
  ))
}
