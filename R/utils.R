#' Plot Cell Type Comparison
#'
#' Creates a side-by-side DimPlot for original clusters and annotated cell types in a Seurat object.
#'
#' @param seurat_obj A Seurat object.
#' @param original_col Character. Metadata column for original clusters (optional; if NULL, uses Idents).
#' @param annotation_col Character. Metadata column for predicted/annotated cell types (default: "annotated_celltype").
#' @param label Logical. Whether to label clusters (default: TRUE).
#' @param pt.size Numeric. Point size (default: 0.3).
#'
#' @return A patchwork ggplot object with two UMAP/tSNE plots.
#' @importFrom Seurat DimPlot Idents
#' @importFrom ggplot2 ggtitle
#' @importFrom patchwork `+`
#' @export
#' @examples
#' # plot_celltype_comparison(seurat_obj, "seurat_clusters", "annotated_celltype")
plot_celltype_comparison <- function(seurat_obj, original_col = NULL, annotation_col = "annotated_celltype",
                                     label = TRUE, pt.size = 0.3) {
  if (!annotation_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Annotation column", annotation_col, "not found in Seurat object metadata."))
  }
  if (is.null(original_col)) {
    p1 <- Seurat::DimPlot(seurat_obj, label = label, pt.size = pt.size) +
      ggplot2::ggtitle("Original Clusters (Idents)") + Seurat::NoLegend()
  } else {
    if (!original_col %in% colnames(seurat_obj@meta.data)) {
      stop(paste("Original group column", original_col, "not found in Seurat object metadata."))
    }
    p1 <- Seurat::DimPlot(seurat_obj, group.by = original_col, label = label, pt.size = pt.size) +
      ggplot2::ggtitle(paste("Original Clusters (", original_col, ")", sep = "")) + Seurat::NoLegend()
  }
  p2 <- Seurat::DimPlot(seurat_obj, group.by = annotation_col, label = label, pt.size = pt.size) +
    ggplot2::ggtitle(paste("Predicted Cell Types (", annotation_col, ")", sep = ""))
  p1 + p2
}

#' Clean and Standardize Cell Type Annotation
#'
#' Cleans, standardizes, and maps raw cell type predictions to Cell Ontology names using direct and fuzzy matching.
#'
#' @param annotation Character. Raw annotation string (possibly mixed or ambiguous).
#' @param mapping_dict Named character vector or mapping table from GPT annotation to CL names.
#' @param ontology_terms Character vector of valid ontology term names.
#'
#' @return Character. Cleaned and mapped annotation string, or the original with "(not in ontology)" if not matched.
#' @importFrom stringr str_to_lower str_remove str_trim str_replace str_detect str_split str_replace_all
#' @export
#' @examples
#' mapping_dict <- c("t cell" = "CL:0000084")
#' ontology_terms <- c("t cell", "b cell", "macrophage")
#' clean_and_match_annotation("t cell", mapping_dict, ontology_terms)
clean_and_match_annotation <- function(annotation, mapping_dict, ontology_terms) {
  annotation <- stringr::str_to_lower(annotation)
  annotation <- stringr::str_remove(annotation, "^\\d+\\.\\s*")
  annotation <- stringr::str_trim(annotation)

  clean_single <- function(candidate) {
    candidate <- stringr::str_trim(candidate)
    if (candidate %in% names(mapping_dict)) return(as.character(mapping_dict[[candidate]]))
    match_ignore_case <- names(mapping_dict)[tolower(names(mapping_dict)) == tolower(candidate)]
    if (length(match_ignore_case) > 0) return(as.character(mapping_dict[[match_ignore_case[1]]]))
    if (candidate %in% ontology_terms) return(candidate)
    candidate_cleaned <- stringr::str_remove(candidate, "\\s*cells?$")
    if (candidate_cleaned %in% ontology_terms) return(candidate_cleaned)
    match_cleaned <- names(mapping_dict)[tolower(names(mapping_dict)) == tolower(candidate_cleaned)]
    if (length(match_cleaned) > 0) return(as.character(mapping_dict[[match_cleaned[1]]]))
    candidate_singular <- stringr::str_replace(candidate, "(?<!e)s$", "")
    candidate_singular_es <- stringr::str_replace(candidate, "es$", "e")
    if (candidate_singular %in% ontology_terms) return(candidate_singular)
    if (candidate_singular_es %in% ontology_terms) return(candidate_singular_es)
    match_singular <- names(mapping_dict)[tolower(names(mapping_dict)) == tolower(candidate_singular)]
    if (length(match_singular) > 0) return(as.character(mapping_dict[[match_singular[1]]]))
    match_singular_es <- names(mapping_dict)[tolower(names(mapping_dict)) == tolower(candidate_singular_es)]
    if (length(match_singular_es) > 0) return(as.character(mapping_dict[[match_singular_es[1]]]))
    candidate_appended <- paste0(candidate, " cell")
    if (candidate_appended %in% ontology_terms) return(candidate_appended)
    match_appended <- names(mapping_dict)[tolower(names(mapping_dict)) == tolower(candidate_appended)]
    if (length(match_appended) > 0) return(as.character(mapping_dict[[match_appended[1]]]))
    suggestions <- tryCatch({ search_ols(candidate) }, error = function(e) NULL)
    if (!is.null(suggestions) && "label" %in% colnames(suggestions)) {
      suggestions <- suggestions[grepl("^CL:", suggestions$obo_id), ]
      if (nrow(suggestions) > 0) {
        return(paste0(suggestions$label[1], collapse = " | "))
      }
    }
    return(paste0(candidate, " (not in ontology)"))
  }
  if (stringr::str_detect(annotation, "/|\\s+or\\s+|\\s+and\\s+|\\(or\\s+[^)]+\\)")) {
    annotation <- stringr::str_replace_all(annotation, "\\(or\\s+([^\\)]+)\\)", "or \\1")
    parts <- unlist(stringr::str_split(annotation, "\\s*(/|or|and)\\s*"))
    parts <- trimws(parts)
    cleaned_parts <- unique(sapply(parts, clean_single, USE.NAMES = FALSE))
    return(paste(cleaned_parts, collapse = " | "))
  }
  clean_single(annotation)
}


#' Build igraph Representation of Ontology from Parent Relationships
#'
#' Given an ontologyIndex object (e.g., as returned by `ontologyIndex::get_ontology()`),
#' constructs an igraph directed graph using 'is_a' parent relationships.
#'
#' @param ontology An ontologyIndex object with a $parents list (e.g., Cell Ontology).
#'
#' @return An igraph object, nodes named by ontology term IDs.
#' @importFrom igraph graph_from_data_frame
#' @export
#' @examples
#' # cl <- ontologyIndex::get_ontology("http://purl.obolibrary.org/obo/cl.obo", extract_tags = "everything")
#' # graph <- build_ontology_graph(cl)
build_ontology_graph <- function(ontology) {
  # Extract edges from parents
  edges <- data.frame(
    from = unlist(ontology$parents),
    to = rep(names(ontology$parents), sapply(ontology$parents, length)),
    stringsAsFactors = FALSE
  )
  edges <- na.omit(edges)
  graph <- igraph::graph_from_data_frame(edges, directed = TRUE)
  return(graph)
}

#' Build Ancestor Map for Cell Ontology
#'
#' For each CL ID, store all parent/ancestor CL IDs.
#' @param cl An ontologyIndex object for Cell Ontology.
#' @return Named list: clid → character vector of ancestor clids (including self)
#' @export
build_ancestor_type_map <- function(cl) {
  stopifnot(requireNamespace("ontologyIndex"))
  # ontologyIndex::get_ancestors(cl, term) returns all ancestors for a term
  ancestor_map <- setNames(
    lapply(names(cl$name), function(id) ontologyIndex::get_ancestors(cl, id)),
    names(cl$name)
  )
  return(ancestor_map)
}

#' Search Cell Ontology Labels and Synonyms via OLS API
#'
#' Queries the EBI Ontology Lookup Service (OLS) API for Cell Ontology (CL) terms matching a given search string.
#' Returns the top results (by label and synonym) as a data.frame with CL labels and OBO IDs.
#'
#' @param query Character. The search query string (cell type name or synonym).
#' @param size Integer. Number of top results to return (default: 3).
#'
#' @return Data.frame of search results with columns: \code{label} and \code{obo_id} (if found); \code{NULL} if no result or API failure.
#'
#' @details
#' This function uses the OLS public API endpoint \url{https://www.ebi.ac.uk/ols/api/search}
#' to search within the Cell Ontology (\code{ontology = "cl"}), using both label and synonym fields.
#' It is primarily intended to provide fuzzy or suggested CL terms when an exact match is not found in a local ontology.
#'
#' @examples
#' \dontrun{
#' # Search for a cell type name
#' search_ols("fibroblast")
#'
#' # Get more suggestions for ambiguous cell types
#' search_ols("helper T", size = 5)
#' }
#'
#' @importFrom httr GET status_code content
#' @importFrom jsonlite fromJSON
#' @export
search_ols <- function(query, size = 3) {
  base_url <- "https://www.ebi.ac.uk/ols/api/search"
  response <- httr::GET(base_url, query = list(
    q = query,
    ontology = "cl",
    size = size,
    queryFields = "label,synonym"
  ))
  if (httr::status_code(response) == 200) {
    results_text <- httr::content(response, as = "text", encoding = "UTF-8")
    results <- jsonlite::fromJSON(results_text)
    docs <- results$response$docs
    if (!is.null(docs) && length(docs) > 0) {
      docs <- as.data.frame(docs)
      docs <- docs[!duplicated(docs$label), ]
      docs <- docs[1:min(nrow(docs), size), ]
      return(docs[, intersect(c("label", "obo_id"), colnames(docs)), drop = FALSE])
    }
  }
  return(NULL)
}
