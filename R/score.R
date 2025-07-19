
# Score the GPT annotation
score_annotation_resolutions <- function(annotation_result_list, output_csv = NULL) {
  if (!is.list(annotation_result_list)) stop("Input must be a named list of results.")
  rng <- function(x) diff(range(x, na.rm = TRUE))
  norm_vec <- function(x) {
    if (rng(x) == 0) rep(1, length(x)) else (x - min(x, na.rm = TRUE)) / rng(x)
  }
  sum_path_length <- sapply(annotation_result_list, function(res) sum(res$final_summary$avg_distance, na.rm = TRUE))
  avg_max_perc <- sapply(annotation_result_list, function(res) mean(res$final_summary$max_percentage, na.rm = TRUE))
  min_max_perc <- sapply(annotation_result_list, function(res) min(res$final_summary$max_percentage, na.rm = TRUE))
  norm_sum <- 1 - norm_vec(sum_path_length)
  norm_avg <- norm_vec(avg_max_perc)
  norm_min <- norm_vec(min_max_perc)
  composite_score <- norm_sum / 2 + norm_avg / 4 + norm_min / 4
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
    write.csv(summary_table, output_csv, row.names = FALSE)
    message("✓  Score summary written to: ", output_csv)
  }
  return(summary_table)
}


extract_synonyms <- function(x) {
  if (length(x) == 0) return(character(0))
  syns <- stringr::str_match(x, '^"([^"]+)"')[,2]
  syns[!is.na(syns)]
}

build_cl_term_map <- function(cl, verbose = TRUE) {
  # stopifnot(requireNamespace("stringr", quietly = TRUE))
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

map_celltypes_to_cl <- function(terms, cl_term_map, verbose = TRUE) {
  terms <- as.character(terms)
  keys <- tolower(terms)
  match_idx <- match(keys, cl_term_map$key)
  clid <- cl_term_map$clid[match_idx]
  cl_label <- cl_term_map$cl_label[match_idx]
  if (verbose) {
    for (i in seq_along(terms)) {
      if (!is.na(clid[i])) {
        message("Mapped: '", terms[i], "' → ", clid[i], " (", cl_label[i], ")")
      } else {
        message("No match: '", terms[i], "'")
      }
    }
  }
  return(data.frame(term = terms, clid = clid, cl_label = cl_label, stringsAsFactors = FALSE))
}

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
        flush.console()
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

# Score the two annotation similarity
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
