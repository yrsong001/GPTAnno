% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/score.R
\name{extract_synonyms}
\alias{extract_synonyms}
\title{Extract Synonyms from OBO-Style Strings}
\usage{
extract_synonyms(x)
}
\arguments{
\item{x}{Character vector of OBO-style synonym fields.}
}
\value{
Character vector of synonyms.
}
\description{
Helper function to extract synonyms from CL ontology synonym strings.
}
\examples{
extract_synonyms('"fibroblast"')
}
