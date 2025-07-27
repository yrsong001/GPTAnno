
### Utility function ###
get_all_descendant_names <- function(ontology, term_name) {
  clean_name <- gsub("_", " ", term_name)
  term_ids <- names(ontology$name[ontology$name == clean_name])
  if (length(term_ids) == 0) stop("Term name not found in the ontology.")
  all_descendant_names <- character()
  for (term_id in term_ids) {
    descendant_ids <- get_descendants(ontology, roots = term_id, exclude_roots = TRUE)
    descendant_names <- ontology$name[descendant_ids]
    all_descendant_names <- c(all_descendant_names, descendant_names)
  }
  unique(unname(all_descendant_names))
}


subcluster_and_find_markers <- function(seurat_obj,
                                        predicted_celltype_column = "new_celltype",
                                        cluster_col = "seurat_clusters",
                                        output_dir = "results/subcluster_markers",
                                        resolutions = c(0.1, 0.2, 0.3, 0.4, 0.5),
                                        dims = 1:30,
                                        assay = NULL,
                                        min_cell_count = 10000) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  meta <- seurat_obj@meta.data
  celltype_stats <- meta %>%
    dplyr::group_by(!!sym(predicted_celltype_column)) %>%
    dplyr::summarize(
      cluster_count = dplyr::n_distinct(!!sym(cluster_col)),
      cell_count = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::filter(cluster_count > 1 | cell_count > min_cell_count)
  message("Cell types to be subclustered: ", paste(celltype_stats[[predicted_celltype_column]], collapse = ", "))
  results <- list()
  for (ct in celltype_stats[[predicted_celltype_column]]) {
    message("\nProcessing: ", ct)
    ct_safe <- gsub("[^[:alnum:]_]+", "_", ct)
    ct_dir <- file.path(output_dir, ct_safe)
    dir.create(ct_dir, showWarnings = FALSE)
    subset_path <- file.path(ct_dir, "seurat_subset.rds")
    if (file.exists(subset_path)) {
      message(" - Loading existing subset: ", subset_path)
      subset_obj <- readRDS(subset_path)
    } else {
      message(" - Creating subset...")
      subset_obj <- subset(seurat_obj, subset = !!sym(predicted_celltype_column) == ct)
    }
    existing_marker_files <- list.files(ct_dir, pattern = "^markers_res_.*\\.rds$", full.names = TRUE)
    existing_resolutions <- gsub("markers_res_(.*)\\.rds", "\\1", basename(existing_marker_files))
    resolutions_to_run <- resolutions[!as.character(resolutions) %in% existing_resolutions]
    if (length(resolutions_to_run) == 0) {
      message(" - Marker files exist for all resolutions. Will update metadata only.")
    } else {
      message(" - Running marker detection for resolutions: ", paste(resolutions_to_run, collapse = ", "))
    }
    subcluster_results <- run_multi_resolution_clustering(
      seurat_obj = subset_obj,
      resolutions = resolutions_to_run,
      result_dir = ct_dir,
      dims = dims,
      assay = assay,
      group.by = "seurat_clusters"
    )
    for (res_name in names(subcluster_results$markers)) {
      marker_df <- subcluster_results$markers[[res_name]]
      # csv_path <- file.path(ct_dir, paste0("markers_", res_name, ".csv"))
      # write.csv(marker_df, file = csv_path, row.names = FALSE)
      rds_path <- file.path(ct_dir, paste0("markers_res_", res_name, ".rds"))
      saveRDS(marker_df, file = rds_path)
    }
    all_resolutions <- union(names(subcluster_results$clusters), existing_resolutions)
    for (res in all_resolutions) {
      colname <- paste0("subcluster_res.", res)
      if (res %in% names(subcluster_results$clusters)) {
        subset_obj@meta.data[[colname]] <- subcluster_results$clusters[[res]]
      } else if (!colname %in% colnames(subset_obj@meta.data)) {
        warning(" - Resolution ", res, " exists as a marker file but not in subset metadata. Skipping column.")
      }
    }
    saveRDS(subset_obj, file = subset_path)
    results[[ct]] <- subcluster_results
  }
  return(results)
}

summarize_gptcelltype_sub <- function(markers, model = 'gpt-4o', tissue_name = "", n_runs = 2, restrict_to = NULL) {
  results_list <- vector("list", n_runs)
  for (i in seq_len(n_runs)) {
    res <- gptcelltype_sub(markers, model = model, tissue_name = tissue_name, restrict_to = restrict_to )
    results_list[[i]] <- res
  }
  combined_results <- dplyr::bind_rows(results_list, .id = "run")
  split_results <- combined_results %>%
    tidyr::pivot_longer(cols = -run, names_to = "cluster", values_to = "annotation") %>%
    dplyr::mutate(
      annotation = stringr::str_to_lower(annotation),
      annotation = stringr::str_remove_all(annotation, "^\\s*(-|\\d+\\.)\\s*"),
      annotation = stringr::str_replace(annotation, "^[-\\s]+", ""),
      annotation = stringr::str_trim(annotation, side = "right"),
      annotation = stringr::str_extract(annotation, "[a-zA-Z].*"),
      standardized_annotation = sapply(annotation, clean_and_match_annotation,
                                       mapping_dict = gpt4_to_clname_mapping,
                                       ontology_terms = unique(cl$name))
    ) %>%
    dplyr::mutate(
      annotation_split = strsplit(as.character(standardized_annotation), "\\s*\\|\\s*")
    ) %>%
    tidyr::unnest(annotation_split) %>%
    dplyr::group_by(cluster, run) %>%
    dplyr::mutate(weight = 1 / length(annotation_split)) %>%
    dplyr::ungroup()
  summary <- split_results %>%
    dplyr::group_by(cluster, annotation_split) %>%
    dplyr::summarize(total_weight = sum(weight), .groups = 'drop') %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(total = sum(total_weight),
                  percentage = (total_weight / total) * 100) %>%
    dplyr::arrange(cluster, dplyr::desc(percentage))
  final_summary <- summary %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarize(
      most_frequent_annotation = annotation_split[1],
      max_percentage = round(percentage[1], 2),
      other_annotations = if (n() > 1) {
        paste(annotation_split[-1], round(percentage[-1], 2), "%", collapse = ", ")
      } else { "" },
      .groups = 'drop'
    )
  return(list(
    combined_results = combined_results,
    summary = summary,
    final_summary = final_summary
  ))
}

gptcelltype_sub <- function(input, tissue_name = NULL, model = 'gpt-4o', topgenenumber = 20, restrict_to = NULL) {
  OPENAI_API_KEY <- Sys.getenv("OPENAI_API_KEY")
  if (OPENAI_API_KEY == "") {
    print("Note: OpenAI API key not found: returning the prompt itself.")
    API.flag <- 0
  } else {
    API.flag <- 1
  }
  if (class(input) == 'list') {
    input <- sapply(input, paste, collapse = ',')
  } else {
    input <- input[input$avg_log2FC > 0, , drop = FALSE]
    input <- tapply(input$gene, list(input$cluster), function(i) paste0(i[1:topgenenumber], collapse = ','))
  }
  if (!API.flag) {
    message <- paste0(
      'Identify cell types of ', tissue_name, ' using the following markers separately for each\n row. ',
      'Only provide the cell type name. Do not show numbers before the name.\n',
      'Some can be a mixture of multiple cell types.\n'
    )
    if (!is.null(restrict_to)) {
      message <- paste0(
        message,
        'Restrict your predictions to the following cell types:\n',
        paste(restrict_to, collapse = ', '), '\n'
      )
    }
    message <- paste0(message, paste0(names(input), ':', unlist(input), collapse = "\n"))
    return(message)
  } else {
    print("Note: OpenAI API key found: returning the cell type annotations.")
    cutnum <- ceiling(length(input) / 30)
    if (cutnum > 1) {
      cid <- as.numeric(cut(1:length(input), cutnum))
    } else {
      cid <- rep(1, length(input))
    }
    allres <- sapply(1:cutnum, function(i) {
      id <- which(cid == i)
      flag <- 0
      while (flag == 0) {
        prompt <- paste0(
          'Identify cell types of ', tissue_name, ' using the following markers separately for each\n row. ',
          'Only provide the cell type name. Do not show numbers before the name.\n',
          'Some can be a mixture of multiple cell types.\n'
        )
        if (!is.null(restrict_to)) {
          prompt <- paste0(
            prompt,
            'Restrict your predictions to the following cell types:\n',
            paste(restrict_to, collapse = ', '), '\n'
          )
        }
        prompt <- paste0(prompt, paste(input[id], collapse = '\n'))
        k <- openai::create_chat_completion(
          model = model,
          message = list(list("role" = "user", "content" = prompt))
        )
        res <- strsplit(k$choices[,'message.content'], '\n')[[1]]
        if (length(res) == length(id))
          flag <- 1
      }
      names(res) <- names(input)[id]
      res
    }, simplify = F)
    return(gsub(',$', '', unlist(allres)))
  }
}

# prev: run_annotation_on_subclusters
anno_subcluster_ontology <- function(base_dir = "output/subclusters",
                                          cl,
                                          mapping_dict,
                                          model = 'gpt-4o',
                                          tissue_name = NULL,
                                          n_runs = 2,
                                          resolutions = c(0.1, 0.3, 0.5),
                                          ontology_graph = NULL,
                                          subcluster_prefix = "subcluster_res.") {
  subtypes <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
  all_results <- list()
  for (ct_dir in subtypes) {
    ct_name <- basename(ct_dir)
    message("\nProcessing annotations for: ", ct_name)
    seurat_path <- file.path(ct_dir, "seurat_subset.rds")
    if (!file.exists(seurat_path)) {
      warning("Subset Seurat object missing for: ", ct_name)
      next
    }
    seurat_obj <- readRDS(seurat_path)
    tissue_name <- paste(gsub("_", " ", ct_name) , "in", tissue_name )
    children_names <- get_all_descendant_names(cl, ct_name)
    annotation_results <- gptanno_sub(
      seurat_obj = seurat_obj,
      resolutions = resolutions,
      cl = cl,
      mapping_dict = mapping_dict,
      marker_path = ct_dir,
      model = model,
      tissue_name = tissue_name,
      n_runs = n_runs,
      save_dir = file.path(ct_dir, "annotation"),
      save_plots = TRUE,
      save_objects = FALSE,
      ontology_graph = ontology_graph,
      subcluster_prefix = subcluster_prefix,
      children_names = children_names
    )
    all_results[[ct_name]] <- annotation_results
  }
  return(all_results)
}
# prev: annotate_subclusters_with_inherited_parent_markers
anno_subcluster_inherit <- function(
    sub_seurat,
    parent_marker_file,
    subcluster_res = "subcluster_res.0.1",
    original_cluster_col = "cluster_res.0.1",
    celltype_col = "celltype_res.0.1",
    top_n = 10,
    model = "gpt-4o",
    tissue_name = "Your tissue name",
    n_runs = 1,
    marker_dir = "output/subclusters",
    ontology_graph = NULL,
    ontology = NULL,
    save_dir = NULL,
    save_plots = TRUE,
    save_objects = FALSE
) {
  res_value <- gsub("subcluster_res\\.", "", subcluster_res)
  marker_filename <- paste0("markers_res_", res_value, ".rds")
  sub_markers_path <- file.path(marker_dir, marker_filename)
  if (!file.exists(sub_markers_path)) stop("Missing subcluster marker file: ", sub_markers_path)
  sub_markers <- readRDS(sub_markers_path)
  if (!file.exists(parent_marker_file)) stop("Missing parent marker file: ", parent_marker_file)
  parent_markers <- readRDS(parent_marker_file)
  Idents(sub_seurat) <- subcluster_res
  sub_ids <- unique(Idents(sub_seurat))
  base_celltype <- gsub(" \\(not in GO\\)", "", unique(sub_seurat[[celltype_col]][,1]))
  gpt_input <- list()
  for (sub_id in sub_ids) {
    sub_cells <- WhichCells(sub_seurat, idents = sub_id)
    parent_cluster_mode <- sub_seurat@meta.data[sub_cells, original_cluster_col] |>
      as.character() |> as.numeric() |> na.omit() |> as.integer() |>
      as.data.frame() |> setNames("cluster") |>
      dplyr::count(cluster) |> dplyr::arrange(desc(n)) |>
      dplyr::slice(1) |> dplyr::pull(cluster)
    parent_top_genes <- parent_markers |>
      dplyr::filter(cluster == parent_cluster_mode, avg_log2FC > 0) |>
      dplyr::arrange(desc(avg_log2FC)) |>
      dplyr::slice(1:top_n) |>
      dplyr::pull(gene)
    sub_top_genes <- sub_markers |>
      dplyr::filter(cluster == sub_id, avg_log2FC > 0) |>
      dplyr::arrange(desc(avg_log2FC)) |>
      dplyr::slice(1:top_n) |>
      dplyr::pull(gene)
    combined_genes <- unique(c(parent_top_genes, sub_top_genes))[1:min(20, length(c(parent_top_genes, sub_top_genes)))]
    gpt_input[[as.character(sub_id)]] <- paste(combined_genes, collapse = ", ")
  }
  annotation_summary <- summarize_gptcelltype_sub(
    markers = gpt_input,
    model = model,
    tissue_name = paste0(base_celltype, " in ", tissue_name),
    n_runs = n_runs
  )
  if (!is.null(ontology_graph)) {
    annotation_summary <- calculate_ontology_distance(
      annotation_summary,
      ontology_graph = ontology_graph,
      cl_term_map
    )
  }
  annotation_summary$cluster <- as.character(annotation_summary$cluster)
  Idents(sub_seurat) <- as.character(Idents(sub_seurat))
  annotated_seurat <- assign_celltype(
    sub_seurat,
    annotation_summary,
    new_celltype = paste0("annotated_sub_", res_value)
  )
  if (!is.null(save_dir)) {
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    if (save_plots) {
      pdf(file.path(save_dir, paste0("annotation_plot_res_", res_value, ".pdf")), width = 18, height = 10)
      print(plot_celltype_comparison(
        annotated_seurat,
        original_col = subcluster_res,
        annotation_col = paste0("annotated_sub_", res_value)
      ))
      dev.off()
    }
    if (save_objects) {
      saveRDS(annotated_seurat, file.path(save_dir, paste0("annotated_seurat_res_", res_value, ".rds")))
    }
    saveRDS(annotation_summary, file.path(save_dir, paste0("annotation_summary_res_", res_value, ".rds")))
  }
  return(list(
    summary = annotation_summary,
    seurat = annotated_seurat
  ))
}

# prev: run_annotation_on_subclusters_v2
anno_subcluster <- function( base_dir = "output/subclusters", cl, mapping_dict, model = "gpt-4o",
    tissue_name = NULL,
    n_runs = 2,
    resolutions = c(0.1, 0.3, 0.5),
    ontology_graph        = NULL,
    subcluster_prefix     = "subcluster_res.",
    strategy              = c("ontology", "marker_inheritance"),
    parent_marker_root    = NULL,
    parent_res_val        = "0.3",
    parent_celltype_col   = "celltype_parent"
) {
  strategy  <- match.arg(strategy)
  subtypes  <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
  all_out   <- list()
  for (ct_dir in subtypes) {
    ct_name  <- basename(ct_dir)
    message("\nProcessing ", ct_name, " (", strategy, ")")
    seurat_path <- file.path(ct_dir, "seurat_subset.rds")
    if (!file.exists(seurat_path)) {
      warning("No Seurat subset in ", ct_dir)
      next
    }
    sub_obj <- readRDS(seurat_path)
    tissue_string <- paste(gsub("_", " ", ct_name), "in", tissue_name)
    if (strategy == "ontology") {
      children <- get_all_descendant_names(cl, ct_name)
      annot <- gptanno_sub(
        seurat_obj        = sub_obj,
        resolutions       = resolutions,
        cl                = cl,
        mapping_dict      = mapping_dict,
        marker_path       = ct_dir,
        model             = model,
        tissue_name        = tissue_string,
        n_runs            = n_runs,
        save_dir          = file.path(ct_dir, "annotation_ontology"),
        save_plots        = TRUE,
        save_objects      = FALSE,
        ontology_graph    = ontology_graph,
        subcluster_prefix = subcluster_prefix,
        children_names    = children
      )
      all_out[[ct_name]] <- annot
      next
    }
    original_cluster_col <- paste0("cluster_res.", parent_res_val)
    parent_marker_file <- file.path(parent_marker_root, paste0("markers_res_", parent_res_val, ".rds"))
    message(paste("Access Parent marker in", parent_marker_file))
    annot_list <- list()
    for (res in resolutions) {
      subcluster_res <- paste0(subcluster_prefix, format(res, nsmall = 1))
      annot_save_dir <- file.path(ct_dir, paste0("annotation_marker_inheritance_res_", res))
      annot <- anno_subcluster_inherit(
        sub_seurat           = sub_obj,
        parent_marker_file   = parent_marker_file,
        subcluster_res       = subcluster_res,
        original_cluster_col = original_cluster_col,
        celltype_col         = parent_celltype_col,
        top_n                = 10,
        model                = model,
        tissue_name           = tissue_string,
        n_runs               = n_runs,
        ontology_graph       = ontology_graph,
        ontology             = cl,
        save_dir             = annot_save_dir,
        save_plots           = TRUE,
        save_objects         = FALSE,
        marker_dir           = ct_dir
      )
      annot_list[[paste0("res_", res)]] <- annot
    }
    score_tbl <- summarize_annotation_scores_subclusters(list(tmp = annot_list))
    best_row  <- score_tbl |>
      dplyr::slice_max(composite_score, n = 1, with_ties = FALSE)
    best_res  <- best_row$resolution
    message("   best subcluster resolution: ", best_res,
            " (score ", round(best_row$composite_score, 3), ")")
    all_out[[ct_name]] <- annot_list[best_res]
  }
  return(all_out)
}

gptanno_sub <- function(seurat_obj, resolutions, cl,  mapping_dict,
                        marker_path = "output/marker_genes", model = 'gpt-4o',
                        tissue_name = NULL, n_runs = 2,
                        save_dir = NULL, save_plots = TRUE,
                        save_objects = FALSE, ontology_graph = NULL,
                        subcluster_prefix = "subcluster_res.", children_names = NULL) {
  results_list <- list()
  for (res in resolutions) {
    message("\nRunning annotation for resolution: ", res)
    col_name <- paste0(subcluster_prefix, format(res, nsmall = 1))
    if (!col_name %in% colnames(seurat_obj@meta.data)) {
      warning("Column ", col_name, " not found in metadata. Skipping.")
      next
    } else {
      message("Setting Idents to ", col_name)
      Idents(seurat_obj) <- col_name
    }
    marker_file <- file.path(marker_path, paste0("markers_res_", res, ".rds"))
    if (!file.exists(marker_file)) {
      warning("Marker file not found for resolution ", res, ". Skipping.")
      next
    }
    markers <- readRDS(marker_file)
    annotation_summary <- summarize_gptcelltype_sub(markers, model = model, tissue_name = tissue_name, n_runs = n_runs, restrict_to = children_names)
    if (!is.null(ontology_graph)) {
      annotation_summary <- calculate_ontology_distance(annotation_summary, ontology_graph = ontology_graph, cl_term_map)
    }
    annotation_summary$cluster <- as.character(annotation_summary$cluster)
    Idents(seurat_obj) <- as.character(Idents(seurat_obj))
    annotated_seurat <- assign_celltype(seurat_obj, annotation_summary, new_celltype = "annotated_sub_celltype")
    if (!is.null(save_dir)) {
      dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
      if (save_plots) {
        pdf(file.path(save_dir, paste0("annotation_plot_res_", res, ".pdf")), width = 18, height = 10)
        print(plot_celltype_comparison(annotated_seurat, original_col = paste0("subcluster_res.", res), annotation_col = "annotated_sub_celltype"))
        dev.off()
      }
      if (save_objects) {
        saveRDS(annotated_seurat, file.path(save_dir, paste0("annotated_seurat_res_", res, ".rds")))
      }
    }
    results_list[[paste0("res_", res)]] <- list(
      summary = annotation_summary,
      seurat = annotated_seurat
    )
  }
  return(results_list)
}

## need to check in new seaneario
summarize_annotation_scores_subclusters <- function(annotation_sub) {
  all_scores <- list()
  for (celltype in names(annotation_sub)) {
    celltype_entry <- annotation_sub[[celltype]]
    res_names <- names(celltype_entry)
    sum_path_length <- numeric(0)
    avg_max_perc <- numeric(0)
    valid_res <- character(0)
    for (res in res_names) {
      res_obj <- celltype_entry[[res]]
      fs <- tryCatch(res_obj$summary$final_summary, error = function(e) NULL)
      if (is.null(fs) || !is.data.frame(fs) || nrow(fs) == 0) next
      if (!"max_percentage" %in% names(fs) || all(is.na(fs$max_percentage))) next
      sum_path_length[res] <- sum(fs$avg_distance, na.rm = TRUE)
      avg_max_perc[res]    <- mean(fs$max_percentage, na.rm = TRUE)
      valid_res <- c(valid_res, res)
    }
    if (length(valid_res) == 0) {
      dummy_res <- if (length(res_names) > 0) res_names[1] else NA_character_
      score_df <- data.frame(
        celltype = celltype,
        resolution = dummy_res,
        sum_path_length = NA_real_,
        avg_max_percentage = NA_real_,
        composite_score = 1,
        stringsAsFactors = FALSE
      )
      all_scores[[celltype]] <- score_df
      next
    }
    # --- New: robust normalization for path (smaller is better), percent (bigger is better) ---
    max_possible_path <- max(sum_path_length, na.rm = TRUE)
    # If max_possible_path == 0 (ideal), set all path_scores = 1
    if (max_possible_path == 0) {
      path_score <- rep(1, length(sum_path_length))
    } else {
      path_score <- 1 - (sum_path_length / max_possible_path)
    }
    percent_score <- avg_max_perc / 100
    composite <- (path_score + percent_score) / 2
    # If both ideal, force composite to 1
    is_ideal <- (sum_path_length == 0) & (avg_max_perc == 100)
    if (any(is_ideal)) composite[is_ideal] <- 1

    score_df <- data.frame(
      celltype = rep(celltype, length(valid_res)),
      resolution = valid_res,
      sum_path_length = sum_path_length,
      avg_max_percentage = avg_max_perc,
      composite_score = composite,
      stringsAsFactors = FALSE
    )
    all_scores[[celltype]] <- score_df
  }
  if (length(all_scores) == 0) return(data.frame())
  final_table <- do.call(rbind, all_scores)
  final_table <- dplyr::arrange(final_table, celltype, dplyr::desc(composite_score))
  return(final_table)
}


assign_best_subcluster_annotations <- function(seurat_obj,
                                               summary_scores_sub,
                                               annotation_sub,
                                               parent_column = "celltype_parent",
                                               subcluster_prefix = "subcluster_res.",
                                               final_colname = "celltype_final") {
  seurat_obj@meta.data[[final_colname]] <- seurat_obj@meta.data[[parent_column]]
  best_res_df <- summary_scores_sub |>
    dplyr::group_by(celltype) |>
    dplyr::slice_max(composite_score, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
  for (i in seq_len(nrow(best_res_df))) {
    celltype <- best_res_df$celltype[i]
    best_res <- best_res_df$resolution[i]
    best_obj <- annotation_sub[[celltype]][[best_res]]
    if (is.null(best_obj)) {
      warning("Skipping ", celltype, " no annotation found for ", best_res)
      next
    }
    annotated_seurat <- best_obj$seurat
    sub_annot_col <- grep("^annotated_sub_", colnames(annotated_seurat@meta.data), value = TRUE)
    if (length(sub_annot_col) != 1) {
      warning("Skipping ", celltype, " could not uniquely identify subcluster annotation column.")
      next
    }
    annotations <- annotated_seurat@meta.data[[sub_annot_col]]
    cell_names <- colnames(annotated_seurat)
    seurat_obj@meta.data[cell_names, final_colname] <- annotations
  }
  return(seurat_obj)
}

assign_inherited_subcluster_annotations <- function(
    full_seurat,
    annotation_sub_inherited,
    subcluster_annotation_col_prefix = "annotated_sub_",
    final_colname = "celltype_subcluster",
    parent_colname = "celltype_parent",
    clean_labels = TRUE
) {
  annotations <- rep(NA, ncol(full_seurat))
  names(annotations) <- colnames(full_seurat)
  for (ct in names(annotation_sub_inherited)) {
    entry <- annotation_sub_inherited[[ct]]
    if (!is.list(entry) || length(entry) == 0) {
      warning("Skipping ", ct, ": invalid annotation entry.")
      next
    }
    best_res <- names(entry)[[1]]
    res_data <- entry[[best_res]]
    if (is.null(res_data$seurat)) {
      warning("Skipping ", ct, ": missing Seurat object in ", best_res)
      next
    }
    seurat_obj <- res_data$seurat
    annot_col <- grep(paste0("^", subcluster_annotation_col_prefix), colnames(seurat_obj@meta.data), value = TRUE)
    if (length(annot_col) != 1) {
      warning("Skipping ", ct, ": expected one annotation column, found ", length(annot_col), " in ", best_res)
      next
    }
    cells <- colnames(seurat_obj)
    annots <- as.character(seurat_obj[[annot_col]][, 1])
    if (clean_labels) {
      annots <- gsub(" \\(not in GO\\)", "", annots)
    }
    annotations[cells] <- annots
  }
  if (!parent_colname %in% colnames(full_seurat@meta.data)) {
    stop("`parent_colname` (", parent_colname, ") not found in Seurat metadata.")
  }
  parent_annots <- as.character(full_seurat[[parent_colname]][, 1])
  if (clean_labels) {
    parent_annots <- gsub(" \\(not in GO\\)", "", parent_annots)
  }
  missing_idx <- which(is.na(annotations))
  annotations[missing_idx] <- parent_annots[missing_idx]
  full_seurat[[final_colname]] <- annotations
  return(full_seurat)
}
