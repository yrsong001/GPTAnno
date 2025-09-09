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
#' @import patchwork
#' @export
#' @examples
#' # plot_celltype_comparison(seurat_obj, "seurat_clusters", "annotated_celltype")
plot_celltype_comparison <- function(seurat_obj, original_col = NULL, annotation_col = "annotated_celltype",
                                     label = TRUE, pt.size = 0.3) {
  if (!annotation_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Annotation column", annotation_col, "not found in Seurat object metadata."))
  }
  reduction_to_use <- if ("umap" %in% Reductions(seurat_obj)) "umap" else NULL

  if (is.null(original_col)) {
    p1 <- Seurat::DimPlot(seurat_obj, label = label, reduction = reduction_to_use, pt.size = pt.size) +
      ggplot2::ggtitle("Original Clusters (Idents)") + Seurat::NoLegend()
  } else {
    if (!original_col %in% colnames(seurat_obj@meta.data)) {
      stop(paste("Original group column", original_col, "not found in Seurat object metadata."))
    }
    p1 <- Seurat::DimPlot(seurat_obj, group.by = original_col, label = label, reduction = reduction_to_use, pt.size = pt.size) +
      ggplot2::ggtitle(paste("Original Clusters (", original_col, ")", sep = "")) + Seurat::NoLegend()
  }
  p2 <- Seurat::DimPlot(seurat_obj, group.by = annotation_col, label = label, reduction = reduction_to_use, pt.size = pt.size) +
    ggplot2::ggtitle(paste("Predicted Cell Types (", annotation_col, ")", sep = "")) + Seurat::NoLegend()
  p1 + p2
}

#' Clean and Standardize Cell Type Annotation
#'
#' Cleans, standardizes, and maps raw cell type predictions to Cell Ontology names using:
#' 1. GPTCelltype_mapping (from GPT-4 annotations to CL names)
#' 2. cl_term_map (synonyms/names to CL IDs and labels)
#'
#' @param annotation Character. Raw annotation string (possibly mixed or ambiguous).
#' @param use_ols Logical. Whether to use OLS search as fallback (default: FALSE). Please use this with Caution.
#'
#' @return Character. Cleaned and mapped annotation string, or the original in lower singular form if not matched.
#' @importFrom stringr str_to_lower str_remove str_trim str_replace str_detect str_split str_replace_all
#' @export
#' @examples
#' clean_and_match_annotation("t cell")
#' clean_and_match_annotation("Helper T cells")
#' clean_and_match_annotation("unknown cell")
clean_and_match_annotation <- function(annotation, use_ols = FALSE) {
  # Access package data directly
  # GPTCelltype_mapping and cl_term_map are available as package data

  annotation <- stringr::str_to_lower(annotation)
  annotation <- stringr::str_remove(annotation, "^\\d+\\.\\s*")
  annotation <- stringr::str_trim(annotation)

  clean_single <- function(candidate) {
    candidate <- stringr::str_trim(candidate)

    # Special case: substitute cardiomyocyte(s) with cardiac muscle cell
    if (grepl("cardiomyocytes?", candidate, ignore.case = TRUE)) {
      candidate <- gsub("cardiomyocytes?", "cardiac muscle cell", candidate, ignore.case = TRUE)
      # message("Substituted cardiomyocyte(s) with 'cardiac muscle cell' in: ", candidate)
    }

    candidate_lower <- tolower(candidate)

    # Step 1: Check GPTCelltype_mapping (GPT-4 annotation -> CL name)
    gpt_match <- names(GPTCelltype_mapping)[tolower(names(GPTCelltype_mapping)) == candidate_lower]
    if (length(gpt_match) > 0) {
      return(as.character(GPTCelltype_mapping[[gpt_match[1]]]))
    }

    # Step 2: Check cl_term_map for exact matches with key (case-insensitive)
    key_match_idx <- which(tolower(cl_term_map$key) == candidate_lower & !is.na(cl_term_map$key))
    if (length(key_match_idx) > 0) {
      return(cl_term_map$key[key_match_idx[1]])
    }

    # Step 3: Check cl_term_map for exact matches with cl_label (case-insensitive)
    label_match_idx <- which(tolower(cl_term_map$cl_label) == candidate_lower & !is.na(cl_term_map$cl_label))
    if (length(label_match_idx) > 0) {
      return(cl_term_map$cl_label[label_match_idx[1]])
    }

    # Step 4: Try cleaning variations
    # Remove "cells?" suffix and try again
    candidate_cleaned <- stringr::str_remove(candidate, "\\s*cells?$")
    if (candidate_cleaned != candidate) {
      candidate_cleaned_lower <- tolower(candidate_cleaned)

      # Check GPTCelltype_mapping with cleaned version
      gpt_match_cleaned <- names(GPTCelltype_mapping)[tolower(names(GPTCelltype_mapping)) == candidate_cleaned_lower]
      if (length(gpt_match_cleaned) > 0) {
        return(as.character(GPTCelltype_mapping[[gpt_match_cleaned[1]]]))
      }

      # Check cl_term_map key with cleaned version
      key_match_cleaned_idx <- which(tolower(cl_term_map$key) == candidate_cleaned_lower & !is.na(cl_term_map$key))
      if (length(key_match_cleaned_idx) > 0) {
        return(cl_term_map$key[key_match_cleaned_idx[1]])
      }

      # Check cl_term_map cl_label with cleaned version
      label_match_cleaned_idx <- which(tolower(cl_term_map$cl_label) == candidate_cleaned_lower & !is.na(cl_term_map$cl_label))
      if (length(label_match_cleaned_idx) > 0) {
        return(cl_term_map$cl_label[label_match_cleaned_idx[1]])
      }
    }

    # Step 5: Try singular forms
    candidate_singular <- stringr::str_replace(candidate, "(?<!e)s$", "")
    candidate_singular_es <- stringr::str_replace(candidate, "es$", "e")

    for (sing_candidate in unique(c(candidate_singular, candidate_singular_es))) {
      if (sing_candidate != candidate) {
        sing_candidate_lower <- tolower(sing_candidate)

        # Check GPTCelltype_mapping with singular
        gpt_match_sing <- names(GPTCelltype_mapping)[tolower(names(GPTCelltype_mapping)) == sing_candidate_lower]
        if (length(gpt_match_sing) > 0) {
          return(as.character(GPTCelltype_mapping[[gpt_match_sing[1]]]))
        }

        # Check cl_term_map key with singular
        key_match_sing_idx <- which(tolower(cl_term_map$key) == sing_candidate_lower & !is.na(cl_term_map$key))
        if (length(key_match_sing_idx) > 0) {
          return(cl_term_map$key[key_match_sing_idx[1]])
        }

        # Check cl_term_map cl_label with singular
        label_match_sing_idx <- which(tolower(cl_term_map$cl_label) == sing_candidate_lower & !is.na(cl_term_map$cl_label))
        if (length(label_match_sing_idx) > 0) {
          return(cl_term_map$cl_label[label_match_sing_idx[1]])
        }
      }
    }

    # Step 6: Try appending " cell" to original candidate
    candidate_appended <- paste0(candidate, " cell")
    candidate_appended_lower <- tolower(candidate_appended)

    # Check GPTCelltype_mapping with appended version
    gpt_match_app <- names(GPTCelltype_mapping)[tolower(names(GPTCelltype_mapping)) == candidate_appended_lower]
    if (length(gpt_match_app) > 0) {
      return(as.character(GPTCelltype_mapping[[gpt_match_app[1]]]))
    }

    # Check cl_term_map key with appended version
    key_match_app_idx <- which(tolower(cl_term_map$key) == candidate_appended_lower & !is.na(cl_term_map$key))
    if (length(key_match_app_idx) > 0) {
      return(cl_term_map$key[key_match_app_idx[1]])
    }

    # Check cl_term_map cl_label with appended version
    label_match_app_idx <- which(tolower(cl_term_map$cl_label) == candidate_appended_lower & !is.na(cl_term_map$cl_label))
    if (length(label_match_app_idx) > 0) {
      return(cl_term_map$cl_label[label_match_app_idx[1]])
    }

    # Step 6.5: NEW - Try appending " cell" to SINGULAR forms
    # This handles cases like "myeloids" -> "myeloid" + " cell" = "myeloid cell"
    for (sing_candidate in unique(c(candidate_singular, candidate_singular_es))) {
      if (sing_candidate != candidate) {
        sing_candidate_with_cell <- paste0(sing_candidate, " cell")
        sing_candidate_with_cell_lower <- tolower(sing_candidate_with_cell)

        # Check GPTCelltype_mapping with singular + " cell"
        gpt_match_sing_cell <- names(GPTCelltype_mapping)[tolower(names(GPTCelltype_mapping)) == sing_candidate_with_cell_lower]
        if (length(gpt_match_sing_cell) > 0) {
          return(as.character(GPTCelltype_mapping[[gpt_match_sing_cell[1]]]))
        }

        # Check cl_term_map key with singular + " cell"
        key_match_sing_cell_idx <- which(tolower(cl_term_map$key) == sing_candidate_with_cell_lower & !is.na(cl_term_map$key))
        if (length(key_match_sing_cell_idx) > 0) {
          return(cl_term_map$key[key_match_sing_cell_idx[1]])
        }

        # Check cl_term_map cl_label with singular + " cell"
        label_match_sing_cell_idx <- which(tolower(cl_term_map$cl_label) == sing_candidate_with_cell_lower & !is.na(cl_term_map$cl_label))
        if (length(label_match_sing_cell_idx) > 0) {
          return(cl_term_map$cl_label[label_match_sing_cell_idx[1]])
        }
      }
    }

    # Step 7: Fallback to OLS search if requested and available
    if (use_ols) {
      suggestions <- tryCatch({
        search_ols(candidate)
      }, error = function(e) NULL)

      if (!is.null(suggestions) && "label" %in% colnames(suggestions)) {
        suggestions <- suggestions[grepl("^CL:", suggestions$obo_id), ]
        if (nrow(suggestions) > 0) {
          message("Found OLS suggestion for '", candidate, "': ", suggestions$label[1])
          return(suggestions$label[1])
        }
      }
    }

    # Step 8: Return candidate in lower singular form with log message
    # Convert to singular form
    candidate_singular <- stringr::str_replace(candidate, "(?<!e)s$", "")
    candidate_singular <- stringr::str_replace(candidate_singular, "es$", "e")
    candidate_final <- tolower(candidate_singular)

    # Log the no-match case
    # message("No ontology match found for '", candidate, "'. Returning in lower singular form: '", candidate_final, "'")

    return(candidate_final)
  }

  # Handle multiple annotations separated by delimiters
  if (stringr::str_detect(annotation, "/|\\s+or\\s+|\\s+and\\s+|\\(or\\s+[^)]+\\)")) {
    # Clean up parenthetical "or" expressions first
    annotation <- stringr::str_replace_all(annotation, "\\(or\\s+([^\\)]+)\\)", "or \\1")

    # Split more carefully, preserving word boundaries
    # Split on: / or "or" or "and" (with word boundaries)
    parts <- unlist(stringr::str_split(annotation, "\\s*/\\s*|\\b\\s+or\\s+\\b|\\b\\s+and\\s+\\b"))
    parts <- trimws(parts)
    parts <- parts[parts != ""]  # Remove empty parts

    cleaned_parts <- unique(sapply(parts, function(part) clean_single(part), USE.NAMES = FALSE))
    return(paste(cleaned_parts, collapse = " | "))
  }

  # Single annotation
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

#' Search Cell Ontology Labels and Synonyms via OLS API with Enhanced Matching
#'
#' Queries the EBI Ontology Lookup Service (OLS) API for Cell Ontology (CL) terms matching a given search string.
#' Uses multiple search strategies and similarity scoring to find the most relevant CL terms, with automatic
#' synonym mapping (e.g., "cardiomyocyte" to "cardiac muscle cell"). Returns top results ranked by relevance.
#'
#' @param query Character. The search query string (cell type name or synonym).
#' @param size Integer. Number of top results to return (default: 3).
#'
#' @return Data.frame of search results with columns: \code{label} and \code{obo_id}, ranked by relevance; \code{NULL} if no result or API failure.
#'
#' @details
#' This function uses the OLS public API endpoint \url{https://www.ebi.ac.uk/ols/api/search}
#' to search within the Cell Ontology (\code{ontology = "cl"}), using both label and synonym fields.
#'
#' The enhanced matching includes:
#' \itemize{
#'   \item Multiple search strategies with automatic synonym replacement
#'   \item Similarity scoring based on word overlap and semantic relevance
#'   \item Intelligent ranking to prioritize the most relevant CL terms
#'   \item Filtering to ensure only Cell Ontology (CL:) terms are returned
#' }
#'
#' It is primarily intended to provide fuzzy or suggested CL terms when an exact match is not found in a local ontology,
#' with improved accuracy for common cell type synonyms.
#'
#' @examples
#' \dontrun{
#' # Search for a cell type name
#' search_ols("fibroblast")
#'
#' # Get more suggestions for ambiguous cell types
#' search_ols("helper T", size = 5)
#'
#' # Automatic synonym matching
#' search_ols("ventricular cardiomyocyte")  # Returns "ventricular cardiac muscle cell"
#' }
#'
#' @importFrom httr GET status_code content
#' @importFrom jsonlite fromJSON
#' @export
search_ols <- function(query, size = 3) {
  base_url <- "https://www.ebi.ac.uk/ols/api/search"

  # Function to calculate similarity score
  calculate_similarity <- function(query_terms, label) {
    query_lower <- tolower(query_terms)
    label_lower <- tolower(label)

    # Split into words
    query_words <- unlist(strsplit(query_lower, "\\s+"))
    label_words <- unlist(strsplit(label_lower, "\\s+"))

    # Calculate different similarity metrics
    exact_match <- identical(query_lower, label_lower)
    word_overlap <- length(intersect(query_words, label_words)) / length(union(query_words, label_words))
    contains_all_query_words <- all(sapply(query_words, function(w) grepl(w, label_lower, fixed = TRUE)))
    partial_matches <- sum(sapply(query_words, function(w) any(grepl(w, label_words, fixed = TRUE))))

    # Weighted score
    score <- 0
    if (exact_match) score <- score + 100
    score <- score + (word_overlap * 50)
    if (contains_all_query_words) score <- score + 25
    score <- score + (partial_matches * 5)

    return(score)
  }

  # Try multiple search strategies
  search_queries <- c(
    query,  # Original query
    gsub("cardiomyocyte", "cardiac muscle cell", query),  # Replace cardiomyocyte with cardiac muscle cell
    gsub("myocyte", "muscle cell", query)  # More general replacement
  )

  all_results <- list()

  for (search_query in unique(search_queries)) {
    response <- httr::GET(base_url, query = list(
      q = search_query,
      ontology = "cl",
      size = 20,  # Get more results to find better matches
      queryFields = "label,synonym"
    ))

    if (httr::status_code(response) == 200) {
      results_text <- httr::content(response, as = "text", encoding = "UTF-8")
      results <- jsonlite::fromJSON(results_text)
      docs <- results$response$docs

      if (!is.null(docs) && length(docs) > 0) {
        docs <- as.data.frame(docs)

        # Filter to only include CL (Cell Ontology) IDs
        if ("obo_id" %in% colnames(docs)) {
          docs <- docs[grepl("^CL:", docs$obo_id), ]
        }

        if (nrow(docs) > 0) {
          all_results[[length(all_results) + 1]] <- docs
        }
      }
    }
  }

  if (length(all_results) > 0) {
    # Combine all results
    combined_docs <- do.call(rbind, all_results)

    # Remove duplicates based on obo_id
    combined_docs <- combined_docs[!duplicated(combined_docs$obo_id), ]

    # Calculate similarity scores
    if ("label" %in% colnames(combined_docs)) {
      combined_docs$similarity_score <- sapply(combined_docs$label,
                                               function(label) calculate_similarity(query, label))

      # Sort by similarity score (descending)
      combined_docs <- combined_docs[order(-combined_docs$similarity_score), ]

      # Select top results
      combined_docs <- combined_docs[1:min(nrow(combined_docs), size), ]

      # Return only the requested columns
      result_cols <- intersect(c("label", "obo_id"), colnames(combined_docs))
      return(combined_docs[, result_cols, drop = FALSE])
    }
  }

  return(NULL)
}

#' Get Common Ancestors of Two Cell Types
#'
#' Finds all common ancestors shared between two cell types, returning both
#' cell type names and their corresponding CL IDs.
#'
#' @param celltype1 Character. First cell type name (e.g., "helper T cell")
#' @param celltype2 Character. Second cell type name (e.g., "cytotoxic T cell")
#' @param cl_term_map Data.frame with mapping from cell type names to CL IDs.
#'   Should have columns: key, clid, cl_label. Defaults to package built-in map.
#' @param cl_ontology The CL ontology object from \code{ontologyIndex::get_ontology()}.
#' @param ancestor_type_map Named list: CL ID to character vector of ancestors (including self).
#'        \strong{Run \code{ancestor_type_map <- build_ancestor_type_map(cl)}} to create the ancestor map.
#' @param verbose Logical. Print mapping information (default: TRUE).
#'
#' @return Data.frame with columns:
#' \describe{
#'   \item{ancestor_name}{Character. Cell type name of common ancestor}
#'   \item{ancestor_clid}{Character. CL ID of common ancestor}
#'   \item{specificity_rank}{Integer. Rank by specificity (1 = most specific common ancestor)}
#' }
#' Returns NULL if either cell type cannot be mapped or no common ancestors found.
#'
#' @examples
#' \dontrun{
#' # Setup required objects
#' cl <- ontologyIndex::get_ontology("http://purl.obolibrary.org/obo/cl.obo",
#'                                   extract_tags = "everything")
#' cl_term_map <- build_cl_term_map(cl)
#' ancestor_type_map <- build_ancestor_type_map(cl)
#'
#' # Get common ancestors
#' common <- get_common_ancestors("helper T cell", "cytotoxic T cell",
#'                               cl_term_map, cl, ancestor_type_map)
#' print(common)
#'
#' # Check B cell vs plasma cell
#' common2 <- get_common_ancestors("B cell", "plasma cell",
#'                                cl_term_map, cl, ancestor_type_map)
#' print(common2)
#' }
#'
#' @export
get_common_ancestors <- function(celltype1,
                                 celltype2,
                                 cl_term_map = GPTAnno::cl_term_map,
                                 cl_ontology,
                                 ancestor_type_map,
                                 verbose = TRUE) {

  # Clean and map cell type names to CL IDs
  clean_celltype1 <- clean_and_match_annotation(celltype1)
  clean_celltype2 <- clean_and_match_annotation(celltype2)

  if (verbose) {
    message("Mapping cell types to CL IDs...")
    message("  '", celltype1, "' -> '", clean_celltype1, "'")
    message("  '", celltype2, "' -> '", clean_celltype2, "'")
  }

  # Map to CL IDs
  map1 <- map_celltypes_to_cl(clean_celltype1, cl_term_map, verbose = verbose)
  map2 <- map_celltypes_to_cl(clean_celltype2, cl_term_map, verbose = verbose)

  clid1 <- map1$clid[1]
  clid2 <- map2$clid[1]

  # Check if mapping was successful
  if (is.na(clid1)) {
    warning("Could not map '", celltype1, "' to a CL ID")
    return(NULL)
  }
  if (is.na(clid2)) {
    warning("Could not map '", celltype2, "' to a CL ID")
    return(NULL)
  }

  if (verbose) {
    message("Mapped CL IDs: ", clid1, " and ", clid2)
  }

  # Check if both cell types exist in ancestor map
  if (!clid1 %in% names(ancestor_type_map)) {
    warning("CL ID '", clid1, "' not found in ancestor type map")
    return(NULL)
  }
  if (!clid2 %in% names(ancestor_type_map)) {
    warning("CL ID '", clid2, "' not found in ancestor type map")
    return(NULL)
  }

  # Get all ancestors for each cell type
  ancestors1 <- ancestor_type_map[[clid1]]
  ancestors2 <- ancestor_type_map[[clid2]]

  if (verbose) {
    message("Cell type 1 has ", length(ancestors1), " ancestors")
    message("Cell type 2 has ", length(ancestors2), " ancestors")
  }

  # Find common ancestors
  common_ancestor_ids <- intersect(ancestors1, ancestors2)

  if (length(common_ancestor_ids) == 0) {
    if (verbose) message("No common ancestors found")
    return(NULL)
  }

  # Filter to keep only CL: terms
  cl_common_ancestors <- common_ancestor_ids[grepl("^CL:", common_ancestor_ids)]

  if (length(cl_common_ancestors) == 0) {
    if (verbose) message("No CL: common ancestors found after filtering")
    return(NULL)
  }

  # Get names for common ancestors
  common_ancestor_names <- character(length(cl_common_ancestors))
  for (i in seq_along(cl_common_ancestors)) {
    clid <- cl_common_ancestors[i]
    if (clid %in% names(cl_ontology$name)) {
      common_ancestor_names[i] <- cl_ontology$name[[clid]]
    } else {
      common_ancestor_names[i] <- paste0("Unknown (", clid, ")")
    }
  }

  # Create result dataframe
  result <- data.frame(
    ancestor_name = common_ancestor_names,
    ancestor_clid = cl_common_ancestors,
    stringsAsFactors = FALSE
  )

  # Add specificity ranking (heuristic: shorter CL ID number = more general)
  # Extract numeric part of CL ID for rough specificity ordering
  numeric_ids <- as.numeric(gsub("CL:", "", result$ancestor_clid))
  result$specificity_rank <- rank(numeric_ids, ties.method = "first")

  # Sort by specificity (most specific first)
  result <- result[order(result$specificity_rank), ]

  message("Found ", nrow(result), " CL: common ancestors")

  return(result)
}
