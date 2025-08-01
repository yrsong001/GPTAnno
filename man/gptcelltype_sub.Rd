% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subcluster_anno.R
\name{gptcelltype_sub}
\alias{gptcelltype_sub}
\title{GPT-based Cell Type Annotation for Subclusters}
\usage{
gptcelltype_sub(
  input,
  tissue_name = NULL,
  model = "gpt-4o",
  topgenenumber = 10,
  restrict_to = NULL
)
}
\arguments{
\item{input}{Named character vector or marker gene data.frame/list per subcluster.}

\item{tissue_name}{Character. Tissue or context for prompt (default: NULL).}

\item{model}{Character. GPT model to use (default: 'gpt-4o').}

\item{topgenenumber}{Number of top marker genes to use (default: 10).}

\item{restrict_to}{Character vector of allowed cell types for GPT output (default: NULL).}
}
\value{
Named character vector of predicted cell types for each subcluster.
}
\description{
Calls GPT model for each subcluster to predict cell types using marker genes.
}
