#' Prebuilt Cell Ontology Name/Synonym-to-ID Map
#'
#' This data object is a precomputed lookup table mapping cell type names and synonyms
#' to Cell Ontology (CL) term IDs and labels, used for annotation and ontology distance calculations.
#' Generated with `build_cl_term_map()` on the official CL ontology (see source).
#'
#' @format A data.frame with columns:
#'   \describe{
#'     \item{key}{Character. Lowercase cell type name or synonym.}
#'     \item{clid}{Character. Cell Ontology term ID (e.g., "CL:0000127").}
#'     \item{cl_label}{Character. Official Cell Ontology label.}
#'   }
#' @source \url{http://purl.obolibrary.org/obo/cl.obo}
#' @seealso \code{\link{build_cl_term_map}}, \code{\link{map_celltypes_to_cl}}
#' @usage data(cl_term_map)
#' @examples
#' data(cl_term_map)
#' head(cl_term_map)
"cl_term_map"

#' Reference Mapping: GPTCelltyp_mapping
#'
#' A curated data frame mapping GPT-4-predicted cell type annotations to standardized Cell Ontology (CL) names,
#' based on the supplementary material from paper: Hou, Wenpin, and Zhicheng Ji. "Assessing GPT-4 for cell type annotation in single-cell RNA-seq analysis." Nature methods 21.8 (2024): 1462-1465.
#'
#' @format A named character vector (or data.frame if you want). Names are GPT-4 annotations, values are CL names.
#' @source \url{https://doi.org/10.1038/s41592-024-02235-4}
#' @seealso \code{\link{clean_and_match_annotation}}, \code{\link{map_celltypes_to_cl}}
#' @usage data(GPTCelltyp_mapping)
#' @examples
#' data(GPTCelltyp_mapping)
#' head(GPTCelltyp_mapping)
"GPTCelltyp_mapping"
