# GPTAnno: PDF to Cell Type Markers Pipeline

A pipeline for extracting cell type and corresponding marker gene information from scientific papers via LLM (e.g., GPT), with quality control and optional ontology-based post-processing.

## Project Overview

This repo provides an automated workflow to extract cell type annotations and their associated marker genes from scientific literature. The pipeline combines PDF text extraction, AI-powered content analysis, and quality control to produce clean, structured datasets.

## Key Features

- **PDF Text Extraction**: Text extraction from scientific PDFs
- **AI-Powered Analysis**: Uses GPT models to identify cell types and marker genes
- **Batch Processing**: Concurrent processing for efficient handling of large documents
- **Caching System**: Reduces API costs with intelligent caching
- **Quality Control**: Multiple validation steps and filtering mechanisms

## Installation

### Prerequisites

1. **Python Environment**: Python 3.8+
2. **Required Packages**:
   ```bash
   pip install -r requirements.txt
   ```

   Or install manually:
   ```bash
   pip install pandas numpy pypdf tqdm openai requests
   ```

3. **API Credentials**: Set up OpenAI API access
   ```bash
   export OPENAI_API_KEY="your-api-key"
   export API_USER_ID="your-user-id"  # Optional
   export API_PASSWORD="your-api-key"  # Alternative to OPENAI_API_KEY
   ```

## Usage

### Step 1: Extract Cell Types and Markers from PDFs

The main extraction pipeline processes PDFs to identify cell types and their marker genes.

#### Interactive Mode
```bash
python paper_extraction_cellNgenes.py
```
This will prompt you to select a PDF from the `papers/` directory.

#### Direct Mode
```bash
python paper_extraction_cellNgenes.py --pdf "path/to/your/paper.pdf" --out "outputs"
```

#### Advanced Options
```bash
python paper_extraction_cellNgenes.py \
    --pdf "papers/example.pdf" \
    --out "outputs" \
    --model "gpt-4o-mini" \
    --limit 100
```

**Parameters:**
- `--pdf`: Path to input PDF file
- `--out`: Output directory (default: `./outputs`)
- `--model`: GPT model to use (default: `gpt-4o-mini`)
- `--limit`: Maximum number of text windows to analyze
- `--papers_dir`: Directory containing PDFs (default: `papers`)

## Pipeline Workflow

### Main Pipeline: PDF Processing (`paper_extraction_cellNgenes.py`)

1. **PDF Text Extraction**
   - Extracts text from all pages of the PDF
   - Handles various PDF formats and encoding issues

2. **Citation Metadata Extraction**
   - Uses GPT to extract paper metadata (author, journal, year)
   - Builds citation information for cell type naming

3. **Text Preprocessing**
   - Splits text into sentences
   - Identifies candidate sentences containing marker gene information
   - Creates overlapping text windows for analysis

4. **Abbreviation Mapping**
   - Builds mapping of cell type abbreviations to full names
   - Uses domain-specific heuristics for cardiac/vascular biology
   - Merges PDF-derived and heuristic mappings

5. **AI-Powered Extraction**
   - Uses GPT models to extract cell type and marker gene pairs
   - Processes text windows in batches for efficiency
   - Implements caching to reduce API costs

6. **Quality Control and Filtering**
   - Removes vague or incomplete entries
   - Filters out cell types with insufficient marker information
   - Applies domain-specific naming rules

7. **Deduplication**
   - Merges entries with identical cell type names
   - Combines marker gene lists for comprehensive coverage

## Output Files

The pipeline generates several output files for each processed paper:

### 1. Reference Data (`1-reference-*.csv`)
- Raw extracted data with source text
- Contains all initially identified cell type and marker pairs
- Includes provenance information

### 2. Removed Entries (`2-removed-*.csv`)
- Entries filtered out during quality control
- Includes reasons for removal
- Useful for manual review and pipeline improvement

### 3. Final Dataset (`3-final-*.csv`)
- Clean, deduplicated cell type and marker data
- Ready for downstream analysis
- Standardized format with consistent naming

### 4. Runtime Statistics (`runtime.json`)
- Processing time breakdown
- Performance metrics
- Configuration details

## Optional Utility

If needed, ontology-based post-processing can still be run separately:

```bash
python filterout_cell_ontology.py --input-dir ./outputs
```

This utility reads `3-final-*.csv` files and writes `4-final-filtered-*.csv` outputs, but it is not part of the main extraction pipeline.

## Configuration

### Environment Variables

```bash
# Required
export OPENAI_API_KEY="your-api-key"

# Optional
export API_USER_ID="your-user-id"
export API_PASSWORD="your-api-key"  # Alternative to OPENAI_API_KEY
export OPENAI_BASE_URL="https://api.openai.com/v1"  # Custom endpoint
export OPENAI_MODEL="gpt-4o-mini"  # Default model
```


## Advanced Features

### Caching System
- Automatic caching of API responses
- Reduces costs for repeated analyses
- Cache stored in `outputs/.cache/`

### Batch Processing
- Concurrent processing of multiple text windows
- Configurable batch sizes and worker counts
- Optimized for large documents

### Quality Control
- Multiple validation layers
- Domain-specific filtering rules
- Comprehensive error handling
