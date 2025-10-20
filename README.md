# GPTAnno

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%3E%3D3.5-blue.svg)](https://www.r-project.org/)

> An Ontology-Guided, LLM-Driven Framework for Robust and Automated Cell Type Annotation in Single-Cell RNA-seq Data

## Overview

**GPTAnno** is the first ontology-referenced LLM pipeline for automated scRNA-seq cell type annotation. Starting from a gene expression matrix, GPTAnno outputs annotations with reproducibility measures, without requiring any external reference dataset. GPTAnno provides:

- **Adaptive clustering with LLM reasoning**: Evaluates multiple clustering resolutions and Selects the optimal resolutions through repeated GPT queries.
- **Cell Ontology (CL) guidance**: Standardizes annotations and enables ontology-aware evaluation.
- **Composite scoring framework**: Combines consistency (agreement across runs), robustness (cluster stability), and reliability (ontology distance) to select optimal resolution.
- **Hierarchical subclustering**: Refines broad categories into functional or disease-relevant subtypes.
- **Flexible subcluster annotation strategies**:
  - **CL-restricted prompting**: Restricts predictions to ontology-defined child terms. This makes sure the predictions has standard annotated cell type names.
  - **Parent marker inheritance**: Combines parent and subcluster markers for prompting, utilizing the flexibility of LLMs to discover emerging cell types beyond current CL coverage.

## Installation

Install from GitHub using `devtools`:

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install GPTAnno
devtools::install_github("yrsong001/GPTAnno")
```

## Prerequisites

1. **OpenAI API Key**: Set your API key as an environment variable:

```r
Sys.setenv(OPENAI_API_KEY = "your-api-key-here")
```

2. **Cell Ontology**: Load the Cell Ontology object and build the ontology graph:

```r
library(ontologyIndex)
cl <- get_ontology("http://purl.obolibrary.org/obo/cl.obo", extract_tags = "everything")
graph <- build_ontology_graph(cl)
```

## Quick Start

Here's a minimal example for parent-level annotation:

```r
library(GPTAnno)
library(Seurat)
library(ontologyIndex)

# Setup
Sys.setenv(OPENAI_API_KEY = "your-api-key")
cl <- get_ontology("http://purl.obolibrary.org/obo/cl.obo", extract_tags = "everything")
graph <- build_ontology_graph(cl)

# Run clustering and find markers
cluster_results <- run_multi_resolution_clustering(
  seurat_obj = your_seurat_object,
  resolutions = c(0.1, 0.3, 0.5),
  result_dir = "output/marker_genes"
)
your_seurat_object <- cluster_results$seurat_obj

# Run repeated GPT queries at multiple resolutions
results <- gptanno(
  seurat_obj = your_seurat_object,
  resolutions = c(0.1, 0.3, 0.5),
  cl = cl,
  graph = graph,
  tissue_name = "mouse heart",
  model = "gpt-4",
  n_runs = 2  # Number of independent queries for reproducibility assessment
)

# Evaluate resolutions with composite scoring
scores <- score_annotation_resolutions(results)
print(scores)  # Highest composite score is optimal resolution
```

## Complete Workflow

### Step 1: Optimal Clustering and Annotation by Multi-Resolution Clustering

Identify major cell types through adaptive clustering at multiple resolutions with automated optimal resolution selection:

```r
# Setup ontology
cl <- get_ontology("http://purl.obolibrary.org/obo/cl.obo", extract_tags = "everything")
graph <- build_ontology_graph(cl)

# Step 1a: Run clustering and find markers
cluster_results <- run_multi_resolution_clustering(
  seurat_obj = seurat_obj,
  resolutions = c(0.1, 0.3, 0.5, 0.7),
  result_dir = "output/marker_genes",
  dims = 1:30
)
seurat_obj <- cluster_results$seurat_obj

# Step 1b: Perform repeated GPT queries for each resolution
parent_results <- gptanno(
  seurat_obj = seurat_obj,
  resolutions = c(0.1, 0.3, 0.5, 0.7),
  cl = cl,
  graph = graph,
  tissue_name = "mouse adult heart",
  model = "gpt-4",
  n_runs = 10,  # Repeated queries for reproducibility measures
  topgenenumber = 10
)

# Step 1c: Score resolutions using composite framework
scores <- score_annotation_resolutions(parent_results)
print(scores)

# Step 1d: Select optimal resolution and assign annotations
best_res_name <- scores$resolution[1]  # e.g., "res_0.3"
best_res_obj <- parent_results[[best_res_name]]
res_val <- sub("^res_", "", best_res_name)  # Extract "0.3"
cluster_col <- paste0("cluster_res.", res_val)  # "cluster_res.0.3"

seurat_obj <- assign_celltype(
  seurat_obj = seurat_obj,
  annotation_summary = best_res_obj,
  cluster_col = cluster_col,
  new_celltype = "celltype_parent"
)
```

### Step 2: Subclustering, Marker Discovery and Annotation with Strategy Selection

Depends on the research needs, subclustering annotated cell type (user speciy or by default criteria) and find marker genes:

```r
subcluster_results <- subcluster_and_find_markers(
  seurat_obj = seurat_obj,
  cl = cl,
  predicted_celltype_column = "celltype_parent",
  output_dir = "output/subclusters",
  resolutions = c(0.1, 0.2, 0.3),
  celltypes_to_subcluster = c("T cell", "B cell", "endothelial cell"), 
  dims = 1:30
)
```
Tip: If celltypes_to_subcluster is not specified, GPTAnno automatically selects cell types that meet minimum cell count criteria and have descendants in the Cell Ontology.


Annotate subclusters using CL-restricted prompting or parent marker inheritance prompting:

- Strategy A: CL-Restricted (Ontology-Constrained), Use for: Well-characterized cell types where biological validity is paramount
- Strategy B: Marker Inheritance (Discovery-Oriented), Use for: Discovering subtypes beyond CL coverage like disease-specific cell states

```r
# Option A: CL-restricted prompting (constrains to Cell Ontology descendants)
annotation_results <- run_subcluster_annotation_workflow(
  base_dir = "output/subclusters",
  strategy = "ontology",
  cl = cl,
  tissue_name = "mouse adult heart",
  resolutions = c(0.1, 0.2, 0.3),
  model = "gpt-4",
  n_runs = 10,
  select_best = TRUE,
  user_restrict_to = list(
    "T cell" = c("regulatory T cell", "effector T cell", "memory T cell"),
    "B cell" = c("memory B cell", "naive B cell", "plasma cell")
  ),
  combine_restrictions = TRUE  # Combine user restrictions with ontology descendants
)

# Option B: Parent marker inheritance (discovers novel subtypes beyond CL)
annotation_results <- run_subcluster_annotation_workflow(
  base_dir = "output/subclusters",
  strategy = "marker_inheritance",
  cl = cl,
  parent_marker_root = "output/markers",
  parent_res = "0.3",
  parent_cluster_col = "cluster_res.0.3",
  tissue_name = "mouse adult heart",
  model = "gpt-4",
  n_runs = 10
)
```

Assign the best subcluster annotations back to your full Seurat object:

```r
seurat_annotated <- assign_subcluster_annotations_to_full(
  full_seurat = seurat_obj,
  workflow_results = annotation_results,
  use_best = TRUE,
  parent_column = "celltype_parent"
)

# View results
table(seurat_annotated$celltype_subcluster)
```

## Key Features

### Automatic Resolution Selection

GPTAnno scores clustering resolutions using a composite framework integrating three metrics:
- **Reliability (Ontology distance)**: Measures semantic coherence through Cell Ontology distance calculations
- **Consistency**: Quantifies agreement across repeated GPT queries to assess reproducibility
- **Robustness (Cluster quality)**: Evaluates stability of cell type predictions across cluster distributions

### Per-Celltype Restrictions

Specify different candidate cell types for each parent celltype using named lists:

```r
user_restrict_to = list(
  "T cell" = c("regulatory T cell", "effector T cell", "CD8+ T cell"),
  "B cell" = c("memory B cell", "naive B cell", "plasma cell"),
  "endothelial cell" = c("vascular endothelial cell", "lymphatic endothelial cell")
)
```

### Automated Name Standardization

All cell type predictions are automatically cleaned and mapped to standardized Cell Ontology terms:
- "Helper T cells" → "T-helper cell"
- "Treg" → "regulatory T cell"
- "cardiomyocyte" → "cardiac muscle cell"

This ensures ontology-grounded evaluation and reproducible nomenclature.

### Two Annotation Strategies

Choose between two complementary approaches for subcluster annotation:

| Strategy | When to Use | Advantages |
|----------|-------------|------------|
| **CL-restricted prompting** | General purpose, well-characterized cell types | Constrains predictions to Cell Ontology descendants, ensures biological validity through ontology structure |
| **Parent marker inheritance** | Cell types beyond CL coverage | Discovers cell types by combining parent and subcluster marker profiles |

## Example Dataset

The package has been tested on the **Aging mouse heart dataset** (non-cardiomyocytes):

```r
# Download dataset (example)
# aging <- readRDS("path/to/Aging_preprocessed.rds")

# Setup
cl <- get_ontology("http://purl.obolibrary.org/obo/cl.obo", extract_tags = "everything")
graph <- build_ontology_graph(cl)

# Run clustering
cluster_results <- run_multi_resolution_clustering(
  seurat_obj = aging,
  resolutions = c(0.1, 0.3, 0.5, 0.7),
  result_dir = "output/marker_genes"
)
aging <- cluster_results$seurat_obj

# Perform repeated GPT queries for reference-free annotation
results <- gptanno(
  seurat_obj = aging,
  resolutions = c(0.1, 0.3, 0.5, 0.7),
  cl = cl,
  graph = graph,
  tissue_name = "mouse adult heart non-cardiomyocytes",
  model = "gpt-4",
  n_runs = 10  # Multiple queries enable reproducibility assessment
)
```

## Documentation

- **Getting Started**: See `vignette("getting-started", package = "GPTAnno")`
- **Subcluster Workflow**: See `vignette("subcluster-annotation", package = "GPTAnno")`
- **Function Reference**: See [package documentation](man/)
- **Complete Example**: See [Annotation_example.R](Annotation_example.R)

## Citation

If you find GPTAnno useful in your research, please cite:

```
Song, Y. (2024). GPTAnno: Automated Cell Type Annotation for Single-Cell
RNA-seq Analysis. R package version 0.0.0.9000.
https://github.com/yrsong001/GPTAnno
```

## Related Work

This package builds upon GPT-4 cell type annotation methods described in:

> Hou, W., & Ji, Z. (2024). Assessing GPT-4 for cell type annotation in
> single-cell RNA-seq analysis. *Nature Methods*, 21(8), 1462-1465.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

## Issues and Support

- **Bug reports**: [GitHub Issues](https://github.com/yrsong001/GPTAnno/issues)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
