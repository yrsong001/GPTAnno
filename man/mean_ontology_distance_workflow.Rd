% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/score.R
\name{mean_ontology_distance_workflow}
\alias{mean_ontology_distance_workflow}
\title{Mean Ontology Distance Workflow}
\usage{
mean_ontology_distance_workflow(
  seurat_obj,
  col1,
  col2,
  cl_term_map,
  graph,
  verbose = TRUE
)
}
\arguments{
\item{seurat_obj}{Seurat object.}

\item{col1}{Column name of first cell type annotation.}

\item{col2}{Column name of second cell type annotation.}

\item{cl_term_map}{Name/synonym-to-ID table from \code{build_cl_term_map}.}

\item{graph}{igraph CL ontology object.}

\item{verbose}{Print mapping/progress.}
}
\value{
List with mean distance and mapping dataframe.
}
\description{
Runs mapping and computes mean ontology distance between two cell annotation columns in a Seurat object.
}
