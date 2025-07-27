#' Query GPT-4 for Cell Type Annotation Using Marker Genes
#'
#' Calls OpenAI's GPT model to predict cell types given marker gene sets for each cluster.
#'
#' @param input A named character vector or data frame: cluster-wise marker genes (e.g., c("Cluster1"="Cd3e, Cd4, Cd8a")) or marker gene data frame.
#' @param tissue_name Character. Name of tissue for context (optional).
#' @param model Character. GPT model to use (default: 'gpt-4o').
#' @param topgenenumber Integer. Number of top marker genes per cluster (default: 10).
#'
#' @return A named character vector: predicted cell types for each cluster.
#' @export
gptcelltype <- function(input, tissue_name=NULL, model='gpt-4', topgenenumber = 10) {
  OPENAI_API_KEY <- Sys.getenv("OPENAI_API_KEY")
  if (OPENAI_API_KEY == "") {
    stop("Error: OpenAI API key not found. Please set the OPENAI_API_KEY environment variable.")
  } else {
    API.flag <- 1
  }

  if (class(input)=='list') {
    input <- sapply(input,paste,collapse=',')
  } else {
    input <- input[input$avg_log2FC > 0,,drop=FALSE]
    input <- tapply(input$gene,list(input$cluster),function(i) paste0(i[1:topgenenumber],collapse=','))
  }

  if (API.flag){
    # print("Note: OpenAI API key found: returning the cell type annotations.")
    cutnum <- ceiling(length(input)/30)
    if (cutnum > 1) {
      cid <- as.numeric(cut(1:length(input),cutnum))
    } else {
      cid <- rep(1,length(input))
    }

    allres <- sapply(1:cutnum,function(i) {
      id <- which(cid==i)
      flag <- 0
      while (flag == 0) {
        k <- openai::create_chat_completion(
          model = model,
          message = list(list("role" = "user", "content" = paste0('Identify cell types of ',tissue_name,' cells using the following markers separately for each\n row. Only provide the cell type name. Do not show numbers before the name.\n Some can be a mixture of multiple cell types.\n',paste(input[id],collapse = '\n'))))
        )
        res <- strsplit(k$choices[,'message.content'],'\n')[[1]]
        if (length(res)==length(id))
          flag <- 1
      }
      names(res) <- names(input)[id]
      res
    },simplify = F)
    #print('Note: It is always recommended to check the results returned by GPT-4 in case of\n AI hallucination, before going to down-stream analysis.')
    return(gsub(',$','',unlist(allres)))
  }

}

#' Summarize GPT-4 Cell Type Annotations Across Multiple Runs
#'
#' Calls GPT model multiple times for cell type annotation and summarizes results.
#'
#' @param markers Named character vector or list of marker genes per cluster.
#' @param model Character. Model to use (default: 'gpt-4o').
#' @param tissue_name Character. Optional context for prompt.
#' @param n_runs Integer. Number of GPT calls to aggregate (default: 2).
#'
#' @return List with: \itemize{
#'   \item combined_results: raw predictions,
#'   \item summary: weighted summary table,
#'   \item final_summary: most frequent annotation per cluster.
#' }
#' @importFrom dplyr bind_rows mutate group_by ungroup summarize arrange
#' @importFrom tidyr pivot_longer unnest
#' @importFrom stringr str_to_lower str_remove_all str_replace str_trim str_extract
#' @export
summarize_gptcelltype <- function(markers, model = 'gpt-4o', tissue_name = "", n_runs = 2) {
  results_list <- vector("list", n_runs)
  for (i in seq_len(n_runs)) {
    res <- gptcelltype(markers, model = model, tissue_name = tissue_name)
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
      standardized_annotation = sapply(annotation, function(x) {
        if (is.na(x)) {
          return(NA_character_)
        } else {
          clean_and_match_annotation(x, mapping_dict = GPTAnno::GPTCelltyp_mapping, ontology_terms = unique(cl$name))
        }
      })
    ) %>%
    dplyr::mutate(annotation_split = strsplit(as.character(standardized_annotation), "\\s*\\|\\s*")) %>%
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


#' Run Annotation Workflow for All Cluster Resolutions
#'
#' Applies GPT annotation workflow for all specified clustering resolutions.
#' For each, retrieves markers, runs annotation, calculates ontology distances, and stores output.
#'
#' @param seurat_obj A Seurat object.
#' @param resolutions Numeric vector of resolutions to annotate.
#' @param cl The ontology object.
#' @param mapping_dict Named character vector or data.frame; mapping from GPT-predicted to CL names. Defaults to package data \code{GPTCelltyp_mapping}.
#' @param model Character. Model to use (default: 'gpt-4o').
#' @param tissue_name Character. Optional context for prompt.
#' @param n_runs Integer. Number of GPT calls to aggregate (default: 2).
#'
#' @return A named list of annotation summary objects for each resolution.
#' @importFrom dplyr arrange
#' @export
gptanno <- function(seurat_obj, resolutions, cl, mapping_dict = GPTAnno::GPTCelltyp_mapping, model = 'gpt-4o',
                    tissue_name = NULL, n_runs = 2) {
  results_list <- list()
  prediction_dir <- "./output/prediction"
  if (!dir.exists(prediction_dir)) {
    message("Prediction directory not found. Creating one at: ", prediction_dir)
    dir.create(prediction_dir, recursive = TRUE)
  }
  for (res in resolutions) {
    message("\nRunning annotation for resolution: ", res)
    col_name <- paste0("cluster_res.", res)
    if (!col_name %in% colnames(seurat_obj@meta.data)) {
      warning("Column ", col_name, " not found in metadata. Skipping.")
      next
    }
    Seurat::Idents(seurat_obj) <- col_name
    marker_file <- paste0("output/marker_genes/markers_res_", res, ".rds")
    if (!file.exists(marker_file)) {
      warning("Marker file not found for resolution ", res, ". Skipping.")
      next
    }
    markers <- readRDS(marker_file)
    annotation_summary <- summarize_gptcelltype(markers, model = model, tissue_name = tissue_name, n_runs = n_runs)
    all_clusters <- unique(seurat_obj@meta.data[[col_name]])
    annotated_clusters <- unique(annotation_summary$summary$cluster)
    missing_clusters <- setdiff(all_clusters, annotated_clusters)
    if (length(missing_clusters) > 0) {
      warning("The following clusters lack annotations and will be labeled 'unannotated': ", paste(missing_clusters, collapse = ", "))
    }
    annotation_summary <- calculate_ontology_distance(annotation_summary, ontology_graph = graph, cl_term_map)
    annotated_seurat <- assign_celltype(seurat_obj, annotation_summary)
    results_list[[paste0("res_", res)]] <- annotation_summary
    pdf(file = file.path(prediction_dir, paste0("res_", res, ".pdf")), width = 30, height = 10)
    print(plot_celltype_comparison(annotated_seurat), group = paste0("res_", res))
    dev.off()
  }
  return(results_list)
}

#' Assign Annotated Cell Types to Seurat Object Metadata
#'
#' Adds a new metadata column to a Seurat object with annotated cell types for each cluster.
#'
#' @param seurat_obj A Seurat object.
#' @param annotation_summary Annotation summary as returned by `summarize_gptcelltype`.
#' @param cluster_col Character. Column to use as cluster identity (default: NULL, use Idents).
#' @param new_celltype Character. Name of new metadata column (default: "annotated_celltype").
#'
#' @return The updated Seurat object with new cell type annotations.
#' @importFrom Seurat Idents
#' @export
assign_celltype <- function(seurat_obj, annotation_summary, cluster_col = NULL, new_celltype = "annotated_celltype") {
  if (!is.null(cluster_col)) {
    if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
      stop("assign_celltype(): column '", cluster_col, "' not found in metadata.")
    }
    Seurat::Idents(seurat_obj) <- cluster_col
  }
  current_idents <- as.character(Seurat::Idents(seurat_obj))
  final_summary <- annotation_summary$final_summary
  if (is.null(final_summary) || !"cluster" %in% colnames(final_summary))
    stop("assign_celltype(): annotation_summary$final_summary is missing or malformed.")
  cluster_annotations <- setNames(final_summary$most_frequent_annotation, as.character(final_summary$cluster))
  all_clusters    <- unique(current_idents)
  missing_clusters <- setdiff(all_clusters, names(cluster_annotations))
  if (length(missing_clusters) > 0) {
    warning("Unannotated clusters: ", paste(missing_clusters, collapse = ", "),
            " labeled 'unannotated'.")
    cluster_annotations[missing_clusters] <- "unannotated"
  }
  seurat_obj@meta.data[[new_celltype]] <- cluster_annotations[current_idents]
  return(seurat_obj)
}

#' Calculate Average Ontology Distance Between Cluster Annotations (using CL term map)
#'
#' Maps annotation names (including synonyms) to CL IDs using `cl_term_map` and computes mean ontology distances for each cluster.
#'
#' @param result_summary List as returned by `summarize_gptcelltype` (must include `summary` and `final_summary`).
#' @param ontology_graph An igraph object representing the Cell Ontology DAG.
#' @param cl_term_map Data frame produced by `build_cl_term_map()`, with columns 'key', 'clid', 'cl_label'.
#'
#' @return The input `result_summary` list, with `final_summary` augmented by `avg_distance` per cluster.
#' @importFrom dplyr left_join mutate group_by ungroup summarise filter
#' @importFrom igraph V shortest_paths
#' @export
calculate_ontology_distance <- function(result_summary, ontology_graph, cl_term_map) {
  # Defensive checks
  if (!inherits(ontology_graph, "igraph")) {
    warning("Invalid ontology graph provided. Skipping distance calculation.")
    result_summary$final_summary$avg_distance <- NA
    return(result_summary)
  }
  if (is.null(result_summary$summary) || !"annotation_split" %in% colnames(result_summary$summary)) {
    stop("result_summary$summary is missing or does not contain 'annotation_split'.")
  }
  if (!"cluster" %in% colnames(result_summary$summary)) {
    stop("result_summary$summary does not contain a 'cluster' column.")
  }
  # Helper: map annotation (synonyms or names) to CL IDs using the term map
  name_to_id <- function(name) {
    i <- match(tolower(name), cl_term_map$key)
    if (is.na(i)) return(NA_character_)
    cl_term_map$clid[[i]]
  }
  mapped_summary <- dplyr::mutate(
    result_summary$summary,
    cl_term = sapply(annotation_split, name_to_id)
  )
  # Compute average ontology distance for each cluster
  distances <- mapped_summary %>%
    dplyr::filter(!is.na(cl_term)) %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(
      avg_distance = {
        terms <- unique(cl_term)
        if (length(terms) < 2) {
          NA_real_
        } else {
          dist_values <- utils::combn(terms, 2, function(x) {
            if (all(x %in% igraph::V(ontology_graph)$name)) {
              path_result <- tryCatch({
                igraph::shortest_paths(ontology_graph, from = x[1], to = x[2], mode = "all")$vpath[[1]]
              }, error = function(e) NULL)
              if (is.null(path_result) || length(path_result) == 0) return(Inf)
              return(length(path_result) - 1)
            } else {
              return(NA_real_)
            }
          }, simplify = TRUE)
          finite_vals <- dist_values[is.finite(dist_values)]
          if (length(finite_vals) == 0) NA_real_ else mean(finite_vals, na.rm = TRUE)
        }
      },
      .groups = "drop"
    )
  result_summary$final_summary <- dplyr::left_join(
    result_summary$final_summary,
    distances,
    by = "cluster"
  )
  return(result_summary)
}

