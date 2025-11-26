# =============================================================================
# Visium HD Spatial Transcriptomics Analysis Pipeline
# Script 04: Cluster Marker Gene Identification
# =============================================================================
#
# Description:
#   This script identifies marker genes specific to each cell type cluster
#   using the Wilcoxon rank-sum test. It generates marker gene tables,
#   expression heatmaps, and exports results for publication.
#
# Input:
#   - Annotated Seurat object from Script 03
#
# Output:
#   - Excel file with cluster markers (all, top N per cluster)
#   - Heatmap of top marker genes
#   - Marker gene statistics
#
# Methods:
#   - FindAllMarkers with Wilcoxon rank-sum test
#   - One-vs-all comparison for each cluster
#   - Multiple testing correction using Bonferroni method
#
# Dependencies:
#   - Seurat (>= 5.0)
#   - dplyr
#   - openxlsx
#   - pheatmap
#   - ComplexHeatmap
#   - circlize
#
# =============================================================================

# -----------------------------------------------------------------------------
# Load Required Packages
# -----------------------------------------------------------------------------

library(Seurat)
library(dplyr)
library(openxlsx)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(Matrix)

set.seed(1234)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Output directory
output_dir <- "results"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Analysis parameters
params <- list(
  # Filtering criteria for marker detection
  min_pct = 0.25,           # Minimum % of cells expressing gene in either group
  logfc_threshold = 0.25,   # Minimum log2 fold-change
  padj_cutoff = 0.05,       # Adjusted p-value threshold
  
  # Output parameters
  top_n_markers = 10,       # Number of top markers per cluster for summary
  top_n_heatmap = 5,        # Number of top markers per cluster for heatmap
  
  # Heatmap parameters
  cells_per_cluster = 100   # Maximum cells per cluster for heatmap
)

# Color palette for cell types (must match Script 03)
cell_type_colors <- c(
  "outer ZF"        = "#E41A1C",
  "ZG"              = "#377EB8",
  "inner ZF"        = "#4DAF4A",
  "X/JMZ"           = "#984EA3",
  "BAT"             = "#FF7F00",
  "CT capsule"      = "#A65628",
  "WAT"             = "#F781BF",
  "medulla"         = "#00CED1",
  "Subcapsular ZG"  = "#FFD92F",
  "fibroblasts"     = "grey"
)

# Sex colors
sex_colors <- c("Male" = "#377EB8", "Female" = "#E41A1C")

# Cluster order for visualization
cluster_order <- c(
  "fibroblasts", "CT capsule", "Subcapsular ZG", "ZG",
  "outer ZF", "inner ZF", "X/JMZ",
  "medulla", "WAT", "BAT"
)

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

#' Create Flat DE Assay from Layered Spatial Assay
#'
#' In Seurat v5, the Spatial assay may have multiple layers. This function
#' creates a single-layer assay suitable for FindMarkers functions.
#'
#' @param obj Seurat object
#' @param spatial_assay Name of the spatial assay
#' @return Assay object with combined counts and normalized data
make_DE_assay <- function(obj, spatial_assay = "Spatial") {
  stopifnot(spatial_assay %in% Assays(obj))
  
  layers <- Layers(obj[[spatial_assay]])
  
  # Get counts layers
  counts_layers <- grep("^counts", layers, value = TRUE)
  if (length(counts_layers) == 0) {
    stop("No 'counts.*' layers found in Spatial assay.")
  }
  
  # Combine counts matrices if multiple layers exist
  counts_list <- lapply(counts_layers, function(layer) {
    LayerData(obj[[spatial_assay]], layer = layer)
  })
  
  if (length(counts_list) == 1) {
    counts_flat <- counts_list[[1]]
  } else {
    counts_flat <- Reduce(`+`, counts_list)
  }
  
  # Check for data layer
  data_layers <- grep("^data$", layers, value = TRUE)
  
  if (length(data_layers) == 1) {
    # Use existing normalized data
    data_mat <- LayerData(obj[[spatial_assay]], layer = "data")
    DE_assay <- CreateAssayObject(counts = counts_flat)
    DE_assay <- SetAssayData(DE_assay, slot = "data", new.data = data_mat)
  } else {
    # Normalize de novo
    DE_assay <- CreateAssayObject(counts = counts_flat)
    obj[["DE"]] <- DE_assay
    obj <- NormalizeData(obj, assay = "DE", normalization.method = "LogNormalize",
                         scale.factor = 10000, verbose = FALSE)
    DE_assay <- obj[["DE"]]
    obj[["DE"]] <- NULL
  }
  
  return(DE_assay)
}

#' Downsample Cells Per Cluster
#'
#' Subsamples cells to a maximum number per cluster for visualization.
#' Maintains cluster order and balances sexes when possible.
#'
#' @param obj Seurat object
#' @param per_cluster Maximum cells per cluster
#' @param order_levels Optional vector specifying cluster order
#' @return Character vector of selected cell names
downsample_cells <- function(obj, per_cluster = 500L, order_levels = NULL) {
  cells_by_ident <- CellsByIdentities(obj)
  cells_by_ident <- cells_by_ident[vapply(cells_by_ident, length, 1L) > 0]
  
  if (!is.null(order_levels)) {
    cells_by_ident <- cells_by_ident[intersect(order_levels, names(cells_by_ident))]
  }
  
  if (!length(cells_by_ident)) return(character(0))
  
  selected <- lapply(cells_by_ident, function(v) {
    if (length(v) > per_cluster) {
      sample(v, per_cluster)
    } else {
      v
    }
  })
  
  unlist(selected, use.names = FALSE)
}

#' Downsample and Order Cells by Sex
#'
#' Subsamples cells while maintaining sex balance within each cluster.
#' Orders cells: Male first, then Female within each cluster.
#'
#' @param obj Seurat object
#' @param n_per_cluster Maximum cells per cluster
#' @param cluster_order Optional vector specifying cluster order
#' @return Character vector of selected cell names (ordered)
downsample_and_order_cells <- function(obj, n_per_cluster = 500,
                                       cluster_order = NULL) {
  current_idents <- Idents(obj)
  sex_vec <- obj$sex
  
  cells_by_cluster <- split(names(current_idents), current_idents)
  
  if (!is.null(cluster_order)) {
    cells_by_cluster <- cells_by_cluster[intersect(cluster_order, 
                                                   names(cells_by_cluster))]
  }
  
  ordered_cells <- lapply(names(cells_by_cluster), function(cl) {
    cells <- cells_by_cluster[[cl]]
    
    # Split by sex
    male_cells <- cells[sex_vec[cells] == "Male"]
    female_cells <- cells[sex_vec[cells] == "Female"]
    
    n_male <- length(male_cells)
    n_female <- length(female_cells)
    total <- n_male + n_female
    
    # Proportional downsampling if needed
    if (total > n_per_cluster) {
      n_male_keep <- round(n_per_cluster * n_male / total)
      n_female_keep <- n_per_cluster - n_male_keep
      
      if (n_male > n_male_keep) {
        male_cells <- sample(male_cells, n_male_keep)
      }
      if (n_female > n_female_keep) {
        female_cells <- sample(female_cells, n_female_keep)
      }
    }
    
    # Order: Male first, then Female
    c(male_cells, female_cells)
  })
  
  unlist(ordered_cells, use.names = FALSE)
}

# -----------------------------------------------------------------------------
# Load Data
# -----------------------------------------------------------------------------

message("\n", rep("=", 80))
message("LOADING ANNOTATED DATA")
message(rep("=", 80))

combined <- readRDS(file.path(output_dir, "combined_annotated.rds"))
message("Loaded object with ", ncol(combined), " cells")

# =============================================================================
# PREPARE FOR MARKER DETECTION
# =============================================================================

message("\n", rep("=", 80))
message("PREPARING FOR MARKER DETECTION")
message(rep("=", 80))

# Create flat DE assay for FindAllMarkers compatibility
combined[["DE"]] <- make_DE_assay(combined, spatial_assay = "Spatial")
DefaultAssay(combined) <- "DE"

# Get clusters present in data
present_clusters <- intersect(cluster_order, levels(Idents(combined)))
if (length(present_clusters) == 0) {
  stop("No clusters from cluster_order found in object!")
}

message("Clusters present: ", paste(present_clusters, collapse = ", "))
message("Total cells: ", ncol(combined))

# Report cell counts per cluster
cells_per_cluster <- table(Idents(combined))
message("\nCells per cluster:")
print(cells_per_cluster)

# =============================================================================
# FIND CLUSTER MARKERS
# =============================================================================

message("\n", rep("=", 80))
message("FINDING CLUSTER-SPECIFIC MARKERS")
message(rep("=", 80))

# FindAllMarkers performs one-vs-all comparisons for each cluster
# Uses Wilcoxon rank-sum test (non-parametric)
message("Running FindAllMarkers (this may take several minutes)...")

all_cluster_markers <- FindAllMarkers(
  combined,
  assay = "DE",
  only.pos = TRUE,              # Only find upregulated markers
  min.pct = params$min_pct,     # Minimum % expressing
  logfc.threshold = params$logfc_threshold,
  test.use = "wilcox",          # Wilcoxon rank-sum test
  verbose = TRUE
)

message("Raw markers found: ", nrow(all_cluster_markers))

# Filter by adjusted p-value
cluster_markers_filtered <- all_cluster_markers %>%
  filter(p_val_adj < params$padj_cutoff) %>%
  group_by(cluster) %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  ungroup()

message("After filtering (padj < ", params$padj_cutoff, "): ", 
        nrow(cluster_markers_filtered))

# Get top N markers per cluster
top_cluster_markers <- cluster_markers_filtered %>%
  group_by(cluster) %>%
  mutate(rank = row_number()) %>%
  filter(rank <= params$top_n_markers) %>%
  ungroup()

message("\nTop ", params$top_n_markers, " markers per cluster:")
print(table(top_cluster_markers$cluster))

# =============================================================================
# EXPORT MARKER RESULTS
# =============================================================================

message("\n", rep("=", 80))
message("EXPORTING MARKER RESULTS")
message(rep("=", 80))

wb_markers <- createWorkbook()

# Sheet 1: Summary statistics
marker_summary <- cluster_markers_filtered %>%
  group_by(cluster) %>%
  summarise(
    n_markers = n(),
    n_cells = as.numeric(cells_per_cluster[cluster[1]]),
    mean_log2FC = round(mean(avg_log2FC), 3),
    median_log2FC = round(median(avg_log2FC), 3),
    min_pval_adj = format(min(p_val_adj), scientific = TRUE, digits = 3),
    .groups = "drop"
  )

addWorksheet(wb_markers, "Summary")
writeData(wb_markers, "Summary", marker_summary)

# Sheet 2: All filtered markers
addWorksheet(wb_markers, "All_markers")
writeData(wb_markers, "All_markers", cluster_markers_filtered)

# Sheet 3: Top N markers per cluster
addWorksheet(wb_markers, paste0("Top_", params$top_n_markers))
writeData(wb_markers, paste0("Top_", params$top_n_markers), top_cluster_markers)

# Individual sheets for each cluster
for (clust in present_clusters) {
  clust_markers <- cluster_markers_filtered %>% filter(cluster == clust)
  if (nrow(clust_markers) > 0) {
    # Sanitize sheet name (Excel has 31 character limit)
    sheet_name <- gsub("[^A-Za-z0-9_]", "_", clust)
    sheet_name <- substr(sheet_name, 1, 31)
    addWorksheet(wb_markers, sheet_name)
    writeData(wb_markers, sheet_name, clust_markers)
  }
}

markers_file <- file.path(output_dir, "cluster_markers_all_cells.xlsx")
saveWorkbook(wb_markers, markers_file, overwrite = TRUE)
message("Exported cluster markers: ", markers_file)

# =============================================================================
# CREATE MARKER HEATMAP
# =============================================================================

message("\n", rep("=", 80))
message("CREATING MARKER HEATMAP")
message(rep("=", 80))

# Get extended marker list (more than needed to handle duplicates)
top_markers_extended <- cluster_markers_filtered %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 15) %>%  # Get top 15 for backup genes
  ungroup()

# Order by cluster
top_markers_extended$cluster <- factor(
  top_markers_extended$cluster, 
  levels = cluster_order
)
top_markers_extended <- top_markers_extended %>% arrange(cluster)

# Select unique genes (avoid duplicates across clusters)
genes_for_heatmap <- c()
cluster_for_gene <- c()
used_genes <- c()

for (clust in present_clusters) {
  cluster_genes <- top_markers_extended %>%
    filter(cluster == clust) %>%
    pull(gene)
  
  # Select up to top_n_heatmap unique genes
  selected <- 0
  for (gene in cluster_genes) {
    if (selected >= params$top_n_heatmap) break
    if (!gene %in% used_genes) {
      genes_for_heatmap <- c(genes_for_heatmap, gene)
      cluster_for_gene <- c(cluster_for_gene, as.character(clust))
      used_genes <- c(used_genes, gene)
      selected <- selected + 1
    }
  }
  
  message("  ", clust, ": selected ", selected, " unique genes")
}

message("Total unique genes for heatmap: ", length(genes_for_heatmap))

# Downsample cells for visualization
cells_for_heatmap <- downsample_and_order_cells(
  combined,
  n_per_cluster = params$cells_per_cluster,
  cluster_order = present_clusters
)
message("Cells for heatmap: ", length(cells_for_heatmap))

# Extract expression matrix
expr_mat <- GetAssayData(
  combined, 
  assay = "Spatial", 
  slot = "data"
)[genes_for_heatmap, cells_for_heatmap, drop = FALSE]

message("Expression matrix: ", nrow(expr_mat), " genes x ", ncol(expr_mat), " cells")

# Calculate z-scores (row-wise standardization)
row_means <- Matrix::rowMeans(expr_mat)
row_sds <- sqrt(Matrix::rowMeans((expr_mat - row_means)^2))
row_sds[row_sds == 0 | !is.finite(row_sds)] <- 1
z_mat <- as.matrix((expr_mat - row_means) / row_sds)

# Create column annotations
cluster_vec <- as.character(Idents(combined)[cells_for_heatmap])
sex_vec <- as.character(combined$sex[cells_for_heatmap])

annotation_col <- data.frame(
  Cluster = factor(cluster_vec, levels = present_clusters),
  Sex = factor(sex_vec, levels = c("Male", "Female")),
  row.names = colnames(z_mat)
)

# Create row annotations
annotation_row <- data.frame(
  Marker_for = factor(cluster_for_gene, levels = present_clusters),
  row.names = genes_for_heatmap
)

# Define colors
annotation_colors <- list(
  Cluster = cell_type_colors[present_clusters],
  Sex = sex_colors,
  Marker_for = cell_type_colors[present_clusters]
)

# Color scale (capped at 98th percentile)
z_limit <- quantile(abs(z_mat[is.finite(z_mat)]), 0.98, na.rm = TRUE)
z_limit <- max(z_limit, 1e-3)
heatmap_colors <- colorRampPalette(c("blue", "black", "yellow"))(100)
heatmap_breaks <- seq(-z_limit, z_limit, length.out = 101)

# Calculate column gaps (between clusters)
cluster_table <- table(annotation_col$Cluster)[present_clusters]
gaps_col <- cumsum(cluster_table)

# Create heatmap
message("Generating heatmap...")

pdf(file.path(output_dir, "heatmap_top_markers.pdf"), width = 10, height = 8)

pheatmap(
  z_mat,
  color = heatmap_colors,
  breaks = heatmap_breaks,
  cluster_rows = FALSE,          # Preserve gene order by cluster
  cluster_cols = FALSE,          # Preserve cell order by cluster
  show_colnames = FALSE,         # Too many cells to show names
  show_rownames = TRUE,
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = annotation_colors,
  fontsize_row = 7,
  fontsize_col = 6,
  gaps_col = gaps_col,           # Visual separation between clusters
  main = "Top Cluster Markers (Z-score)",
  border_color = NA,
  legend = TRUE
)

dev.off()

message("Saved: heatmap_top_markers.pdf")

# =============================================================================
# SUMMARY
# =============================================================================

message("\n", rep("=", 80))
message("MARKER ANALYSIS COMPLETE")
message(rep("=", 80))

message("\nGenerated files:")
message("  - cluster_markers_all_cells.xlsx (marker tables)")
message("  - heatmap_top_markers.pdf")

# Print final gene selection summary
message("\n=== SELECTED MARKER GENES ===")
for (clust in present_clusters) {
  genes_in_clust <- genes_for_heatmap[cluster_for_gene == clust]
  message(clust, ": ", paste(genes_in_clust, collapse = ", "))
}