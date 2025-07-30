#' Score GPT Annotation Results Across Multiple Resolutions
#'
#' Calculates composite scores for annotation results at different resolutions,
#' considering ontology path length and annotation confidence.
#'
#' @param annotation_result_list A named list of annotation results (e.g., output from `run_annotation_all_resolutions()`).
#' @param output_csv Optional. File path to write the summary table as CSV.
#'
#' @return A data.frame summarizing each resolution's metrics and a composite score (1 is ideal: short ontology distance, high percentage).
#' @examples
#' # result_list <- list(res_01 = ..., res_02 = ...)  # Your annotation results
#' # score_annotation_resolutions(result_list)
#' @export
score_annotation_resolutions <- function(annotation_result_list, output_csv = NULL) {
  if (!is.list(annotation_result_list)) stop("Input must be a named list of results.")
  rng <- function(x) diff(range(x, na.rm = TRUE))
  norm_vec <- function(x) {
    if (rng(x) == 0) rep(1, length(x)) else (x - min(x, na.rm = TRUE)) / rng(x)
  }
  sum_path_length <- sapply(annotation_result_list, function(res) sum(res$final_summary$avg_distance, na.rm = TRUE))
  avg_max_perc    <- sapply(annotation_result_list, function(res) mean(res$final_summary$max_percentage, na.rm = TRUE))
  min_max_perc    <- sapply(annotation_result_list, function(res) min(res$final_summary$max_percentage, na.rm = TRUE))

  # Composite score: 1 is ideal (shortest ontology distance, max % = 100)
  norm_sum <- 1 - norm_vec(sum_path_length) # smaller is better
  norm_avg <-      norm_vec(avg_max_perc)   # larger is better
  norm_min <-      norm_vec(min_max_perc)   # larger is better
  composite_score <- (norm_sum + norm_avg + norm_min) / 3

  summary_table <- data.frame(
    resolution          = names(sum_path_length),
    sum_path_length     = sum_path_length,
    avg_max_percentage  = avg_max_perc,
    min_max_percentage  = min_max_perc,
    composite_score     = composite_score,
    row.names           = NULL
  ) |>
    dplyr::arrange(dplyr::desc(composite_score))
  if (!is.null(output_csv)) {
    utils::write.csv(summary_table, output_csv, row.names = FALSE)
    message("✓  Score summary written to: ", output_csv)
  }
  return(summary_table)
}

#' Extract Synonyms from OBO-Style Strings
#'
#' Helper function to extract synonyms from CL ontology synonym strings.
#' @param x Character vector of OBO-style synonym fields.
#' @return Character vector of synonyms.
#' @examples
#' extract_synonyms('"fibroblast"')
#' @export
extract_synonyms <- function(x) {
  if (length(x) == 0) return(character(0))
  syns <- stringr::str_match(x, '^"([^"]+)"')[,2]
  syns[!is.na(syns)]
}

#' Build Name/Synonym-to-CL-ID Map from Ontology Object
#'
#' Constructs a lookup table mapping CL names and synonyms to CL IDs.
#'
#' @param cl An ontologyIndex CL ontology object.
#' @param verbose Logical; print mapping stats.
#' @return Data.frame with columns: key, clid, cl_label.
#' @examples
#' # cl <- ontologyIndex::get_ontology("http://purl.obolibrary.org/obo/cl.obo", extract_tags = "everything")
#' # cl_term_map <- build_cl_term_map(cl)
#' @export
build_cl_term_map <- function(cl, verbose = TRUE) {
  all_ids <- names(cl$name)
  all_labels <- cl$name
  all_labels_lower <- tolower(all_labels)
  df <- data.frame(key = all_labels_lower, clid = all_ids, cl_label = all_labels, stringsAsFactors = FALSE)
  for (clid in names(cl$synonym)) {
    syn_vec <- extract_synonyms(cl$synonym[[clid]])
    if (length(syn_vec)) {
      syn_vec <- tolower(syn_vec)
      syn_vec <- setdiff(syn_vec, df$key)
      if (length(syn_vec)) {
        df <- rbind(
          df,
          data.frame(key = syn_vec, clid = clid, cl_label = cl$name[[clid]], stringsAsFactors = FALSE)
        )
      }
    }
  }
  if (verbose) message("Mapping includes ", nrow(df), " name/synonym entries.")
  df <- unique(df)
  rownames(df) <- NULL
  return(df)
}

#' Map Cell Type Names to CL IDs
#'
#' Maps character vector of cell type names to CL ontology IDs using a prebuilt mapping table.
#'
#' @param terms Character vector of cell type names.
#' @param cl_term_map Data.frame; defaults to package built-in map, or returned by `build_cl_term_map()`.
#' @param verbose Print mapping results.
#' @return Data.frame with columns: term, clid, cl_label.
#' @examples
#' # cell_names <- c("fibroblast", "endothelial cell")
#' # map_celltypes_to_cl(cell_names, cl_term_map)
#' @export
map_celltypes_to_cl <- function(terms, cl_term_map = GPTAnno::cl_term_map, verbose = TRUE) {
  terms <- as.character(terms)
  keys <- tolower(terms)
  match_idx <- match(keys, cl_term_map$key)
  clid <- cl_term_map$clid[match_idx]
  cl_label <- cl_term_map$cl_label[match_idx]
  if (verbose) {
    for (i in seq_along(terms)) {
      if (!is.na(clid[i])) {
        message("Mapped: '", terms[i], "' to ", clid[i], " (", cl_label[i], ")")
      } else {
        message("No match: '", terms[i], "'")
      }
    }
  }
  return(data.frame(term = terms, clid = clid, cl_label = cl_label, stringsAsFactors = FALSE))
}

#' Calculate Mean Shortest Path Distance in CL Ontology
#'
#' For each row of a CL ID dataframe, calculates the ontology graph shortest path distance between pairs.
#'
#' @param clid_df Data.frame with columns clid1 and clid2 (CL IDs).
#' @param graph igraph object representing CL ontology DAG.
#' @param verbose Print progress and stats.
#' @return List: mean_distance, inverse_mean, dist_vector.
#' @importFrom igraph V shortest_paths
#' @importFrom utils flush.console
#' @examples
#' # res <- calculate_mean_ontology_distance(clid_df, graph)
#' @export
calculate_mean_ontology_distance <- function(clid_df, graph, verbose = TRUE) {
  stopifnot(requireNamespace("igraph", quietly = TRUE))
  n <- nrow(clid_df)
  dist_vector <- rep(NA_real_, n)
  last_pct <- -1
  for (i in seq_len(n)) {
    id1 <- clid_df$clid1[i]
    id2 <- clid_df$clid2[i]
    dist <- tryCatch({
      if (is.na(id1) || is.na(id2)) NA_real_
      else if (id1 %in% igraph::V(graph)$name && id2 %in% igraph::V(graph)$name) {
        sp <- igraph::shortest_paths(graph, from = id1, to = id2, mode = "all")$vpath[[1]]
        if (length(sp) > 0) length(sp) - 1 else NA_real_
      } else NA_real_
    }, error = function(e) NA_real_)
    dist_vector[i] <- dist
    if (verbose && n > 1) {
      pct <- floor(100 * i / n)
      if (pct != last_pct) {
        cat(sprintf("\rProgress: %d%%", pct))
        utils::flush.console()
        last_pct <- pct
      }
    }
  }
  if (verbose && n > 1) cat("\n")
  mean_distance <- mean(dist_vector, na.rm = TRUE)
  inverse_mean <- 1 / mean_distance
  if (verbose) {
    message("Mean ontology distance: ", mean_distance)
    message("Inverse (1 / mean): ", inverse_mean)
  }
  return(list(mean_distance = mean_distance, inverse_mean = inverse_mean, dist_vector = dist_vector))
}

#' Mean Ontology Distance Workflow
#'
#' Runs mapping and computes mean ontology distance between two cell annotation columns in a Seurat object.
#'
#' @param seurat_obj Seurat object.
#' @param col1 Column name of first cell type annotation.
#' @param col2 Column name of second cell type annotation.
#' @param cl_term_map Name/synonym-to-ID table from `build_cl_term_map`.
#' @param graph igraph CL ontology object.
#' @param verbose Print mapping/progress.
#' @return List with mean distance and mapping dataframe.
#' @export
mean_ontology_distance_workflow <- function(seurat_obj, col1, col2, cl_term_map, graph, verbose = TRUE) {
  meta <- seurat_obj@meta.data
  meta_subset <- meta[!is.na(meta[[col1]]) & !is.na(meta[[col2]]), ]
  terms1 <- as.character(meta_subset[[col1]])
  terms2 <- as.character(meta_subset[[col2]])
  map1 <- map_celltypes_to_cl(terms1, cl_term_map, verbose = FALSE)
  map2 <- map_celltypes_to_cl(terms2, cl_term_map, verbose = FALSE)
  clid_df <- data.frame(
    cell = rownames(meta_subset),
    term1 = terms1, term2 = terms2,
    clid1 = map1$clid, clid2 = map2$clid,
    label1 = map1$cl_label, label2 = map2$cl_label,
    stringsAsFactors = FALSE
  )
  dist_res <- calculate_mean_ontology_distance(clid_df, graph, verbose = verbose)
  dist_res$clid_df <- clid_df
  return(dist_res)
}


#' Score Agreement Based on Cell Ontology Ancestry
#'
#' Fully match: same CL ID. Partially match: manual is ancestor of predicted or vice versa.
#'
#' @param seurat_obj A Seurat object with both annotation columns.
#' @param manual_col Character. Metadata column for manual annotation.
#' @param predicted_col Character. Metadata column for predicted annotation.
#' @param cl_term_map Data.frame. Term name/synonym to CL ID.
#' @param ancestor_type_map Named list: CL ID to character vector of ancestors (including self).
#'        \strong{Run \code{ancestor_type_map <- build_ancestor_type_map(cl)}} to create the ancestor map before using this function.
#' @param output_csv Optional. Save the detailed score of each cells. Path for CSV export. if NULL, no
#' @return List with per-cell scores and summary.
#' @details
#' \strong{NOTE:} You need to first run \code{ancestor_type_map <- build_ancestor_type_map(cl)}
#' to create the ancestor mapping for your ontology.
#' @export
score_annotation_agreement_ontology <- function(
    seurat_obj,
    manual_col      = "manual_celltype",
    predicted_col   = "predicted_celltype",
    cl_term_map,
    ancestor_type_map,
    output_csv      = NULL
) {
  meta <- seurat_obj@meta.data
  meta_sub <- meta[!is.na(meta[[manual_col]]) & !is.na(meta[[predicted_col]]), ]
  if (nrow(meta_sub) == 0) stop("No overlapping annotated cells found.")

  map1 <- map_celltypes_to_cl(meta_sub[[manual_col]], cl_term_map, verbose = FALSE)
  map2 <- map_celltypes_to_cl(meta_sub[[predicted_col]], cl_term_map, verbose = FALSE)

  # Helper: check ancestry
  is_ancestor <- function(a, b) {
    if (is.na(a) || is.na(b)) return(FALSE)
    a %in% ancestor_type_map[[b]]
  }

  # Agreement score: 1 (exact), 0.5 (ancestor/descendant), 0 (unrelated)
  scores <- mapply(function(c1, c2) {
    if (is.na(c1) || is.na(c2)) return(NA_real_)
    if (c1 == c2) return(1)
    if (is_ancestor(c1, c2) || is_ancestor(c2, c1)) return(0.5)
    return(0)
  }, map1$clid, map2$clid)

  result_df <- data.frame(
    cell = rownames(meta_sub),
    manual_label = meta_sub[[manual_col]],
    predicted_label = meta_sub[[predicted_col]],
    manual_clid = map1$clid,
    predicted_clid = map2$clid,
    agreement_score = scores,
    stringsAsFactors = FALSE
  )

  summary <- list(
    mean_score    = mean(scores, na.rm = TRUE),
    n_cells       = sum(!is.na(scores)),
    fully_matched = mean(scores == 1, na.rm = TRUE),
    partial_match = mean(scores == 0.5, na.rm = TRUE),
    mismatch      = mean(scores == 0, na.rm = TRUE)
  )

  if (!is.null(output_csv)) {
    write.csv(result_df, file = output_csv, row.names = FALSE)
    message("Agreement scores written to: ", output_csv)
  }

  return(list(scores = result_df, summary = summary))
}

