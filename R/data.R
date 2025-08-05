#' Prebuilt Cell Ontology Name/Synonym-to-ID Map
#'
#' A precomputed lookup table mapping cell type names and synonyms
#' to Cell Ontology (CL) term IDs and labels. Used for cell type annotation
#' and ontology distance calculations. Generated via \code{\link{build_cl_term_map}}
#' on the official Cell Ontology.
#'
#' @docType data
#' @name cl_term_map
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

#' GPTCelltyp_mapping: Reference Mapping of GPT-4 Annotations to Cell Ontology
#'
#' A curated mapping from GPT-4-predicted cell type names to standardized Cell Ontology (CL) terms,
#' based on the supplementary material from:
#' Hou, Wenpin, and Zhicheng Ji. "Assessing GPT-4 for cell type annotation in single-cell RNA-seq analysis." \emph{Nature Methods} 21.8 (2024): 1462-1465.
#'
#' @docType data
#' @name GPTCelltyp_mapping
#' @format A named character vector. Names are GPT-4 annotation results, values are CL names.
#' @source \url{https://doi.org/10.1038/s41592-024-02235-4}
#' @seealso \code{\link{clean_and_match_annotation}}, \code{\link{map_celltypes_to_cl}}
#' @usage data(GPTCelltyp_mapping)
#' @examples
#' data(GPTCelltyp_mapping)
#' head(GPTCelltyp_mapping)
"GPTCelltyp_mapping"
