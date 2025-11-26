# Visium HD Spatial Transcriptomics Analysis Pipeline

## Overview

This repository contains R scripts for analyzing 10x Genomics Visium HD spatial transcriptomics data from mouse adrenal glands. The pipeline performs:

1. Data loading and preprocessing
2. Sample integration using Harmony
3. Cell type annotation
4. Cluster marker identification
5. Sex-differential expression analysis

## Study Description

Spatial transcriptomics analysis of mouse adrenal glands comparing male and female samples using 10x Genomics Visium HD technology. The analysis identifies zonally-specific cell populations and investigates sex-specific gene expression patterns across adrenal cortical and medullary regions.

## Directory Structure

```
visium_hd_analysis/
├── scripts/
│   ├── 01_data_loading_preprocessing.R    # Load Visium HD data, QC
│   ├── 02_merge_harmony_integration.R     # Sample merging, batch correction
│   ├── 03_cluster_annotation_visualization.R  # Cell type annotation, plots
│   ├── 04_cluster_marker_analysis.R       # Find cluster markers
│   ├── 05_sex_differential_expression.R   # Male vs Female DE analysis
│   ├── 06_cellchat_communication_analysis.R   # CellChat per sex
│   └── 07_cellchat_comparative_analysis.R     # CellChat comparison M vs F
├── data/                                  # Input data (not included)
├── results/                               # Output files
│   ├── *.rds                             # Seurat objects
│   ├── *.xlsx                            # DE results tables
│   ├── *.pdf                             # Visualizations
│   ├── volcano_plots_sex_DE/             # Volcano plots
│   └── cellchat/                         # CellChat results
│       ├── female/                       # Female CellChat
│       ├── male/                         # Male CellChat
│       ├── comparison/                   # Comparative analysis
│       └── diagnostics/                  # QC tables
└── README.md
```

## Requirements

### R Version
- R >= 4.3.0

### Required Packages

```r
# Core analysis
install.packages(c("Seurat", "SeuratObject", "Matrix"))

# Integration
install.packages("harmony")

# Visualization
install.packages(c("ggplot2", "patchwork", "pheatmap", 
                   "ComplexHeatmap", "viridis", "RColorBrewer"))

# Differential expression
BiocManager::install("edgeR")
install.packages("EnhancedVolcano")

# Cell-cell communication
devtools::install_github("jinworks/CellChat")  # CellChat v2.2.0+
install.packages(c("RANN", "future"))

# Data handling
install.packages(c("dplyr", "tibble", "openxlsx", "readr"))

# Other
install.packages(c("circlize", "cowplot", "ggrepel", "scales", "png"))
```

## Usage

### Quick Start

```r
# Set your data directory
base_dir <- "/path/to/visium_hd_data"

# Run scripts in order
source("scripts/01_data_loading_preprocessing.R")
source("scripts/02_merge_harmony_integration.R")
source("scripts/03_cluster_annotation_visualization.R")
source("scripts/04_cluster_marker_analysis.R")
source("scripts/05_sex_differential_expression.R")
```

### Input Data Structure

Each sample directory should contain:
```
sample_id/
├── segmented_outputs/
│   ├── filtered_feature_cell_matrix.h5
│   └── spatial/
│       ├── tissue_positions.csv
│       ├── scalefactors_json.json
│       ├── tissue_lowres_image.png
│       └── tissue_hires_image.png (optional)
└── spatial/
    └── (backup spatial files)
```

### Sample Configuration

Edit the sample configuration in `01_data_loading_preprocessing.R`:

```r
# Define sample identifiers and sex mapping
samples <- c("SAMPLE1", "SAMPLE2")
sample_to_sex <- c("SAMPLE1" = "Male", "SAMPLE2" = "Female")
```

## Analysis Pipeline

### Script 01: Data Loading and Preprocessing
- Loads 10x Genomics Visium HD data (H5 format)
- Creates Seurat objects with spatial images
- Calculates QC metrics (mitochondrial %, ribosomal %)
- Handles barcode format differences between count matrix and spatial files

### Script 02: Sample Merging and Harmony Integration
- Merges multiple samples into single Seurat object
- Handles Seurat v5 layer structure
- Performs Harmony batch correction on sample variable
- Runs PCA, UMAP, and Louvain clustering

### Script 03: Cluster Annotation and Visualization
- Filters low-quality clusters
- Annotates clusters with cell type names based on marker expression
- Generates UMAP and spatial visualizations
- Creates cell composition plots

### Script 04: Cluster Marker Analysis
- Identifies marker genes for each cluster using Wilcoxon rank-sum test
- Exports marker gene lists to Excel
- Creates expression heatmaps

### Script 05: Sex-Differential Expression Analysis
- Performs pseudobulk differential expression (Male vs Female)
- Analyzes both whole-adrenal and per-cluster differences
- Uses edgeR with quasi-likelihood F-test
- Generates volcano plots for each comparison

### Script 06: CellChat Communication Analysis
- Builds CellChat objects for male and female samples separately
- Implements intelligent downsampling with biological exceptions
- Calibrates spatial parameters for Visium HD cell segmentation
- Computes distance-weighted communication probabilities
- Uses mouse secreted signaling database

### Script 07: CellChat Comparative Analysis
- Compares cell-cell communication between sexes
- Identifies conserved vs sex-specific signaling pathways
- Generates differential interaction heatmaps
- Visualizes pathway-level differences
- Exports communication statistics tables

## Cell Type Annotations

The following adrenal gland cell populations are identified:

| Cluster | Cell Type | Description |
|---------|-----------|-------------|
| ZG | Zona Glomerulosa | Mineralocorticoid-producing cells |
| Subcapsular ZG | Subcapsular Zona Glomerulosa | Outer ZG cells |
| outer ZF | Outer Zona Fasciculata | Glucocorticoid-producing cells |
| inner ZF | Inner Zona Fasciculata | Glucocorticoid-producing cells |
| X/JMZ | X-zone/Juxtamedullary Zone | Specialized cortical zone |
| medulla | Adrenal Medulla | Chromaffin cells (catecholamines) |
| CT capsule | Connective Tissue Capsule | Structural support |
| WAT | White Adipose Tissue | Energy storage |
| BAT | Brown Adipose Tissue | Thermogenesis |
| fibroblasts | Fibroblasts | Stromal cells |

## Key Parameters

### Quality Control
- `min.cells = 3`: Genes expressed in at least 3 cells
- `min.features = 200`: Cells with at least 200 genes

### Integration
- `nfeatures = 3000`: Variable features for integration
- `resolution = 0.35`: Louvain clustering resolution

### Marker Detection
- `min.pct = 0.25`: Minimum percentage of cells expressing marker
- `logfc.threshold = 0.25`: Minimum log2 fold-change
- `p_val_adj < 0.05`: Significance threshold (Bonferroni-adjusted)

### Differential Expression
- `min_pct_within_cluster = 0.20`: Per-cluster expression filter
- `specificity_cutoff = 0.20`: Cluster specificity filter
- `disp_fallback = 0.1`: Fixed dispersion for 2-sample comparisons

### CellChat Analysis
- `contact_multiplier = 2.5`: Contact range (~2.5x cell diameter)
- `interaction_multiplier = 15.0`: Paracrine range (~15x cell diameter)
- `min_cells_per_cluster = 50`: Minimum cells for downsampling
- Database: CellChatDB.mouse (Secreted Signaling)

## Output Files

### Seurat Objects
- `combined_integrated_harmony.rds`: After integration
- `combined_annotated.rds`: With cell type annotations

### Excel Files
- `cluster_markers_all_cells.xlsx`: Cluster marker genes
- `DE_sex_per_cluster.xlsx`: Per-cluster sex DE results
- `DE_sex_whole_adrenal.xlsx`: Whole-adrenal sex DE results

### Visualizations
- `umap_clusters.pdf`: UMAP with cluster annotations
- `umap_by_sex.pdf`: UMAP colored by sex
- `spatial_clusters.pdf`: Spatial distribution of cell types
- `heatmap_top_markers.pdf`: Marker gene expression heatmap
- `volcano_plots_sex_DE/`: Volcano plots for DE results

### CellChat Outputs
- `cellchat/cellchat_female_final.rds`: Female CellChat object
- `cellchat/cellchat_male_final.rds`: Male CellChat object
- `cellchat/comparison/`: Comparative analysis results
- `cellchat/comparison/pathway_statistics.csv`: Pathway-level statistics
- `cellchat/comparison/communication_*.csv`: L-R pair communication tables


1. Korsunsky I, et al. (2019). Fast, sensitive and accurate integration of single-cell data with Harmony. Nature Methods, 16(12), 1289-1296.

2. Hao Y, et al. (2021). Integrated analysis of multimodal single-cell data. Cell, 184(13), 3573-3587.

3. Robinson MD, et al. (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics, 26(1), 139-140.

4. Jin S, et al. (2021). Inference and analysis of cell-cell communication using CellChat. Nature Communications, 12(1), 1088.

5. Jin S, et al. (2024). CellChat for systematic analysis of cell-cell communication from single-cell and spatially resolved transcriptomics. Nature Protocols.
