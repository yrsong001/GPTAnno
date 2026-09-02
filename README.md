# GPTAnno

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%3E%3D3.5-blue.svg)](https://www.r-project.org/)

> An Ontology-Guided, LLM-Driven Framework for Robust and Automated Cell Type Annotation in Single-Cell RNA-seq Data

## Overview

**GPTAnno** is the first ontology-referenced LLM pipeline for automated scRNA-seq cell type annotation. Starting from a gene expression matrix, GPTAnno outputs annotations with reproducibility measures, without requiring any external reference dataset. GPTAnno provides:

- **Adaptive clustering with LLM reasoning**: Evaluates multiple clustering resolutions and selects the optimal resolutions through repeated GPT queries.
- **Cell Ontology (CL) guidance**: Standardizes annotations and enables ontology-aware evaluation.
- **Composite scoring framework**: Combines consistency (agreement across runs), robustness (cluster stability), and reliability (ontology distance) to select optimal resolution.
- **Hierarchical subclustering**: Refines broad categories into functional or disease-relevant subtypes.
- **Flexible subcluster annotation strategies**:
  - **CL-restricted prompting**: Restricts predictions to ontology-defined child terms so that predictions have standard annotated cell type names.
  - **Parent marker inheritance**: Combines parent and subcluster markers for prompting, utilizing the flexibility of LLMs to discover emerging cell types beyond current CL coverage.
- **Multiple LLM backends**: OpenAI, Anthropic, Google Gemini, and local models (Ollama, vLLM) via the [ellmer](https://github.com/ropensci/ellmer) package.
- **Python Tools: PDF to Markers Extraction**: Standalone pipeline for extracting cell type markers from scientific papers (see [`pdf2markers/`](pdf2markers/))

## Installation

Install from GitHub using `devtools`. LLM calls require the **ellmer** package, which will be installed as a dependency.

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install GPTAnno
devtools::install_github("yrsong001/GPTAnno")
```

Optional speedup for large datasets:

```r
install.packages("presto")
```

GPTAnno works without `presto`. If it is installed, Seurat can use it automatically to speed up Wilcoxon-based marker detection in `FindAllMarkers()`.

## Prerequisites

1. **LLM access** (choose one or more):
   - **OpenAI**: Set `OPENAI_API_KEY` (e.g. `Sys.setenv(OPENAI_API_KEY = "your-api-key-here")`).
   - **Anthropic**: Set `ANTHROPIC_API_KEY`.
   - **Google Gemini**: Set `GOOGLE_API_KEY` or `GEMINI_API_KEY`.
   - **Local models (Ollama or vLLM)**: No API key needed; use `llm_config` with `provider = "ollama"` or `provider = "vllm"` (see [Using local LLMs (Ollama)](#using-local-llms-ollama) later in this README).

2. **Cell Ontology**: Load the Cell Ontology object and build the ontology graph:

```r
library(ontologyIndex)
cl <- get_ontology("http://purl.obolibrary.org/obo/cl.obo", extract_tags = "everything")
graph <- build_ontology_graph(cl)
```

## Vignettes

For a fuller walkthrough, start with:

- [`GPTAnno Pipeline`](vignettes/gptanno-pipeline.Rmd): main parent annotation and optional subclustering workflow.
- [`Helper Utilities in GPTAnno`](vignettes/helper-utilities.Rmd): ontology lookup, relationship checks, agreement scoring, and other helper functions.

If you install GPTAnno with vignettes built, you can also open them from R:

```r
devtools::install_github("yrsong001/GPTAnno", build_vignettes = TRUE)
browseVignettes("GPTAnno")
```

## Quick Start

Here's a minimal example for general annotation. Use the same directory for `result_dir` (clustering) and `marker_dir` (gptanno) so that marker files are found.

```r
library(GPTAnno)
library(Seurat)
library(ontologyIndex)
library(dplyr)

# Setup
Sys.setenv(OPENAI_API_KEY = "your-api-key")
cl <- get_ontology("http://purl.obolibrary.org/obo/cl.obo", extract_tags = "everything")
graph <- build_ontology_graph(cl)

# Run clustering and find markers (writes markers to result_dir)
marker_dir <- "output/marker_genes"
cluster_results <- run_multi_resolution_clustering(
  seurat_obj = your_seurat_object,
  resolutions = c(0.1, 0.3, 0.5),
  result_dir = marker_dir
)
your_seurat_object <- cluster_results$seurat_obj

# Run repeated GPT queries at multiple resolutions (reads markers from marker_dir)
results <- gptanno(
  seurat_obj = your_seurat_object,
  resolutions = c(0.1, 0.3, 0.5),
  cl = cl,
  graph = graph,
  tissue_name = "mouse heart",
  llm_config = list(
    provider = "openai",
    model = "gpt-5" # Optionally can change params, eg. params = ellmer::params(temperature = 1)
  ),
  marker_dir = marker_dir,
  n_runs = 2  # Number of independent queries for reproducibility assessment
)

# Evaluate resolutions with composite scoring
scores <- score_annotation_resolutions(results)
print(scores)  # Highest composite score is optimal resolution
```

**Using another LLM provider:** Use the same `llm_config` pattern and set the appropriate provider and model (and API key). For example, with Anthropic:

```r
results <- gptanno(
  seurat_obj = your_seurat_object,
  resolutions = c(0.1, 0.3, 0.5),
  cl = cl,
  graph = graph,
  tissue_name = "mouse heart",
  llm_config = list(
    provider = "anthropic",
    model = "claude-sonnet-4-20250514",
    params = ellmer::params(temperature = 0.2, top_p = 0.9)
  ),
  marker_dir = marker_dir,
  n_runs = 2
)
```

`llm_config$params` is passed directly to `ellmer` chat clients. This is the recommended way to control decoding settings (for example `temperature`, `top_p`, `top_k`, and other provider-supported parameters).

### Parameter compatibility by model (important)

Not every model supports every decoding parameter. For example, some OpenAI GPT models may reject specific fields (like `temperature`) and return HTTP 400 with a message such as:

`Unsupported parameter: 'temperature' is not supported with this model.`

GPTAnno now surfaces this provider message in batch logs (including provider/model/params) when `gptanno()` or `gptcelltype()` calls fail.

How to check support for a model:

1. List models available to your API account:

```r
ellmer::models_openai()
```

2. Probe your target model with a minimal call and the exact params you plan to use:

```r
call_llm(
  prompt = "Reply with OK",
  provider = "openai",
  model = "gpt-5",
  params = ellmer::params(temperature = 0)
)
```

3. If it fails, remove unsupported fields from `llm_config$params` or switch to a compatible model.

There is no single cross-provider parameter-compatibility table in GPTAnno; the authoritative source is each provider's model documentation plus a quick probe call like above.

Set the corresponding API key (e.g. `Sys.setenv(ANTHROPIC_API_KEY = "your-key")`). For local Ollama, use `llm_config = list(provider = "ollama", model = "llama2")` (see [Using local LLMs (Ollama)](#using-local-llms-ollama)).

## Annotation Workflow

### Preprocessing (optional)

If your Seurat object is not yet normalized or has no PCA/UMAP, run `preprocess_seurat_object()` before clustering:

```r
seurat_obj <- preprocess_seurat_object(seurat_obj, npcs = 30, save_path = NULL)  # or save_path = "path/to/preprocessed.rds"
```

### Step 1: Optimal Clustering and Annotation by Multi-Resolution Clustering

Identify major cell types through adaptive clustering at multiple resolutions with automated optimal resolution selection:

```r
# Setup ontology
cl <- get_ontology("http://purl.obolibrary.org/obo/cl.obo", extract_tags = "everything")
graph <- build_ontology_graph(cl)

# Step 1a: Run clustering and find markers (use same path for Step 1b marker_dir)
marker_dir <- "output/marker_genes"
cluster_results <- run_multi_resolution_clustering(
  seurat_obj = seurat_obj,
  resolutions = c(0.1, 0.3, 0.5, 0.7),
  result_dir = marker_dir,
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
  llm_config = list(provider = "openai", model = "gpt-5"),  # or "gpt-4", etc.
  marker_dir = marker_dir,
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

Depends on the research needs, subclustering annotated cell type (user specify or by default criteria: has decendants on CL and exceeds the minimum number of cells, 10,000 by default) and find marker genes:

```r
subcluster_results <- subcluster_and_find_markers(
  seurat_obj = seurat_obj,
  cl = cl,
  predicted_celltype_column = "celltype_parent",
  output_dir = "output/subclusters",
  resolutions = c(0.1, 0.2, 0.3),
  celltypes_to_subcluster = c("T cell", "B cell", "endothelial cell"), # if not specify it will select the cell types to be subclustered by default criteria
  dims = 1:30,
  reduction = "pca"
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
  llm_config = list(provider = "openai", model = "gpt-5"),  # or "gpt-4", etc.
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
  llm_config = list(provider = "openai", model = "gpt-5"),  # or "gpt-4", etc.
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

For other LLM providers or local models, pass the same `llm_config` into `run_subcluster_annotation_workflow()` and `assign_subcluster_annotations_to_full()`.

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

### Automated Cell Name Standardization

All cell type predictions are automatically cleaned and mapped to standardized Cell Ontology terms:
- "Helper T cells" → "T-helper cell"
- "Treg" → "regulatory T cell"
- "cardiomyocyte" → "cardiac muscle cell"

This ensures ontology-grounded evaluation and reproducible nomenclature.

### Evaluating annotations against manual labels

When you have a manual or reference annotation column, you can score agreement with GPTAnno predictions using the Cell Ontology. Build an ancestor map once, then call either the simplified or detailed scoring function:

- **`score_annotation_agreement_ontology()`** — Returns per-cell scores 1.0 (exact or manual is ancestor of predicted), 0.5 (manual is child of predicted), or 0.0 (unrelated). Useful for summary metrics.
- **`score_annotation_agreement_ontology_detailed()`** — Classifies each pair as exact, parent, child, sibling, or no_match and returns counts and pairwise tables.

Both require **`build_ancestor_type_map(cl)`** before use. For distance-based comparison use **`score_annotation_distance_ontology()`**; for a single manual vs predicted CL ID pair use **`check_cl_relationship()`**.

```r
ancestor_type_map <- build_ancestor_type_map(cl)
agreement <- score_annotation_agreement_ontology(
  seurat_obj,
  manual_col = "celltype_manual",
  predicted_col = "celltype_parent",
  ancestor_type_map = ancestor_type_map
)
print(agreement$summary)
```

### Two Subcluster Annotation Strategies

Choose between two complementary approaches for subcluster annotation:

| Strategy | When to Use | Advantages |
|----------|-------------|------------|
| **CL-restricted prompting** | General purpose, well-characterized cell types | Constrains predictions to Cell Ontology descendants, ensures biological validity through ontology structure |
| **Parent marker inheritance** | Cell types beyond CL coverage | Discovers cell types by combining parent and subcluster marker profiles |

**Ontology** strategy restricts predictions to CL descendants. If you also pass `user_restrict_to`, then `combine_restrictions` controls the allowed set: `TRUE` (default) = union of descendants and your list; `FALSE` = only your list. **Marker inheritance** combines parent and subcluster marker profiles and does not require CL descendants. Use the same `llm_config` in `run_subcluster_annotation_workflow()` for non-OpenAI or local models.

## Python Tools: PDF to Markers Extraction

GPTAnno includes a standalone Python pipeline for extracting cell type markers from scientific papers:

**Location**: `pdf2markers/`

**Features**:
- Extract cell types and marker genes from PDF papers using GPT
- Automated ontology filtering to focus on novel discoveries
- Batch processing with caching for cost efficiency
- Quality control and deduplication

**Quick Start**:
```bash
cd pdf2markers

# Install Python dependencies
pip install -r requirements.txt

# Set OpenAI API key
export OPENAI_API_KEY="your-api-key"

# Extract markers from a paper
python paper_extraction_cellNgenes.py --pdf "papers/example.pdf" --out "outputs"

# Filter out ontology-annotated cell types
python filterout_cell_ontology.py --input-dir ./outputs
```

See [`pdf2markers/README.md`](pdf2markers/README.md) for complete documentation.


## Example Dataset

The package has been tested on the **Aging mouse heart dataset** (non-cardiomyocytes):

```r
# Download dataset (example)
# aging <- readRDS("path/to/Aging_preprocessed.rds")

# Setup
cl <- get_ontology("http://purl.obolibrary.org/obo/cl.obo", extract_tags = "everything")
graph <- build_ontology_graph(cl)

# Run clustering (use same path for gptanno marker_dir)
marker_dir <- "output/marker_genes"
cluster_results <- run_multi_resolution_clustering(
  seurat_obj = aging,
  resolutions = c(0.1, 0.3, 0.5, 0.7),
  result_dir = marker_dir
)
aging <- cluster_results$seurat_obj

# Perform repeated GPT queries for reference-free annotation
results <- gptanno(
  seurat_obj = aging,
  resolutions = c(0.1, 0.3, 0.5, 0.7),
  cl = cl,
  graph = graph,
  tissue_name = "mouse adult heart non-cardiomyocytes",
  llm_config = list(provider = "openai", model = "gpt-5"),  # or "gpt-4", etc.
  marker_dir = marker_dir,
  n_runs = 10  # Multiple queries enable reproducibility assessment
)
```

## Using local LLMs (Ollama)

To use [Ollama](https://ollama.ai) locally instead of a cloud API, install Ollama and pull a model, then start the server:

```bash
# In a terminal:
ollama pull llama3    # or mistral, neural-chat, etc.
ollama serve          # default: http://localhost:11434
```

In R, check that Ollama is available and list models: `check_ollama_available()`, `list_ollama_models()`. Then pass the provider and model into any workflow that accepts `llm_config`:

```r
llm_config = list(
  provider = "ollama",
  model = "llama3",
  api_url = "http://localhost:11434"
)
# Use the same llm_config in gptanno(), run_subcluster_annotation_workflow(), etc.
```

## Citation

If you find GPTAnno useful in your research, please cite our bioRxiv preprint:

> Song, Y., Tang, M., Liu, Q., Wang, H., Qian, L., Zou, F., & Hou, W. (2025). GPTAnno: Ontology-tree-guided hierarchical cell type annotation based on GPT models for single-cell data. *bioRxiv*. https://doi.org/10.1101/2025.11.27.690951

- **bioRxiv**: [https://www.biorxiv.org/content/10.1101/2025.11.27.690951v1](https://www.biorxiv.org/content/10.1101/2025.11.27.690951v1)

## Related Work

This package builds upon GPT-4 cell type annotation methods described in:

> Hou, W., & Ji, Z. (2024). Assessing GPT-4 for cell type annotation in
> single-cell RNA-seq analysis. *Nature Methods*, 21(8), 1462-1465.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

## Issues and Support

- **Bug reports**: [GitHub Issues](https://github.com/yrsong001/GPTAnno/issues)
