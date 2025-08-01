% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess.R
\name{run_multi_resolution_clustering}
\alias{run_multi_resolution_clustering}
\title{Run Multi-Resolution Clustering and Marker Detection}
\usage{
run_multi_resolution_clustering(
  seurat_obj,
  resolutions,
  result_dir = "output/markers",
  dims = 1:30,
  assay = NULL,
  group.by = "seurat_clusters",
  reduction = "pca",
  use_existing_neighbors = TRUE
)
}
\arguments{
\item{seurat_obj}{A Seurat object.}

\item{resolutions}{Numeric vector of resolution parameters (e.g., c(0.1, 0.2, 0.3)).}

\item{result_dir}{Character. Directory to save marker RDS files (default: "results/markers").}

\item{dims}{Numeric vector. Which dimensions to use for clustering (default: 1:30).}

\item{assay}{Character. Assay to use, or NULL for current (default: NULL).}

\item{group.by}{Character. Metadata column to use for group assignment (default: "seurat_clusters").}

\item{reduction}{Character. Reduction to use for graph construction (default: "pca").}

\item{use_existing_neighbors}{Logical. Use an existing neighbor graph if present (default: TRUE).}
}
\value{
A list with: updated Seurat object, marker lists, and cluster assignments for each resolution.
}
\description{
This function runs graph-based clustering at multiple resolutions and finds marker genes for each, saving results to disk.
}
\examples{
# result <- run_multi_resolution_clustering(seurat_obj, c(0.1, 0.2, 0.3))
}
