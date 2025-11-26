# =============================================================================
# Visium HD Spatial Transcriptomics Analysis Pipeline
# Script 02: Sample Merging and Harmony Integration
# =============================================================================
#
# Description:
#   This script merges individual Seurat objects from multiple samples,
#   performs layer joining for Seurat v5 compatibility, normalizes data,
#   and runs Harmony batch correction to integrate samples while preserving
#   biological variation.
#
# Input:
#   - List of Seurat objects (objects_list) from Script 01
#   - Or: Individual RDS files for each sample
#
# Output:
#   - Merged and integrated Seurat object with Harmony embeddings
#   - UMAP coordinates for visualization
#
# Dependencies:
#   - Seurat (>= 5.0)
#   - harmony
#   - Matrix
#   - dplyr
#
# References:
#   - Korsunsky et al. (2019) Fast, sensitive, and accurate integration of 
#     single cell data with Harmony. Nature Methods.
#

# =============================================================================

# -----------------------------------------------------------------------------
# Load Required Packages
# -----------------------------------------------------------------------------

library(Seurat)
library(SeuratObject)
library(harmony)
library(Matrix)
library(dplyr)

set.seed(1234)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Define sample order and sex mapping
samples <- c("CCAA2905", "CCAA2105")
sample_to_sex <- c("CCAA2905" = "Male", "CCAA2105" = "Female")

# Output directory
output_dir <- "results"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Harmony parameters
harmony_dims <- 1:50       # PCA dimensions to use for Harmony
harmony_theta <- 2         # Diversity clustering penalty parameter
harmony_lambda <- 1        # Ridge regression penalty parameter

# Clustering parameters
cluster_resolution <- 0.35  # Louvain clustering resolution

# -----------------------------------------------------------------------------
# Load Data (if not already in memory)
# -----------------------------------------------------------------------------

# Uncomment if loading from saved RDS files:
# objects_list <- lapply(samples, function(s) {
#   readRDS(file.path(output_dir, paste0(s, "_seurat.rds")))
# })
# names(objects_list) <- samples

# Verify objects_list exists
stopifnot(exists("objects_list") && length(objects_list) > 0)

# =============================================================================
# STEP 2: MERGE SAMPLES
# =============================================================================

message("\n", rep("=", 80))
message("STEP 2: MERGING SAMPLES")
message(rep("=", 80))

# Merge all Seurat objects
# add.cell.ids: Prefixes cell barcodes with sample ID to ensure uniqueness
combined <- merge(
  x = objects_list[[1]],
  y = objects_list[-1],
  add.cell.ids = samples,
  project = "VisiumHD_AdrenalGland_Combined"
)

# Clean up memory
rm(objects_list)
gc()

# Print merge summary
message("\nMerged object summary:")
message("  - Total cells: ", ncol(combined))
message("  - Total features: ", nrow(combined))
message("  - Assays: ", paste(Assays(combined), collapse = ", "))
message("  - Images: ", paste(Images(combined), collapse = ", "))

# Verify sex distribution
message("\nSex distribution across samples:")
print(table(combined$sex, combined$sample))

# =============================================================================
# STEP 3: JOIN LAYERS (SEURAT V5 SPECIFIC)
# =============================================================================

message("\n", rep("=", 80))
message("STEP 3: JOINING LAYERS")
message(rep("=", 80))

# In Seurat v5, merged objects maintain separate layers per sample.
# We need to join these layers for downstream analysis.

DefaultAssay(combined) <- "Spatial"

# Get current layer structure
layers <- Layers(combined[["Spatial"]])
message("Current layers: ", paste(layers, collapse = ", "))

# Identify counts layers (one per sample after merge)
counts_layers <- grep("^counts", layers, value = TRUE)
message("Found ", length(counts_layers), " counts layers to join")

# Extract all count matrices
counts_list <- lapply(counts_layers, function(layer) {
  LayerData(combined[["Spatial"]], layer = layer)
})

# Report layer dimensions
message("\nLayer dimensions:")
for (i in seq_along(counts_list)) {
  message("  ", counts_layers[i], ": ", 
          nrow(counts_list[[i]]), " genes x ", 
          ncol(counts_list[[i]]), " cells")
}

# Get union of all genes across layers
# Some genes may be present in one sample but not another
all_genes <- Reduce(union, lapply(counts_list, rownames))
message("\nTotal unique genes across all layers: ", length(all_genes))

# Align all matrices to have the same gene set
# Missing genes are filled with zeros (sparse matrix for efficiency)
message("Aligning gene sets across layers...")

aligned_counts <- lapply(seq_along(counts_list), function(i) {
  mat <- counts_list[[i]]
  
  # Create empty sparse matrix with all genes
  aligned_mat <- Matrix::Matrix(
    data = 0,
    nrow = length(all_genes),
    ncol = ncol(mat),
    sparse = TRUE
  )
  rownames(aligned_mat) <- all_genes
  colnames(aligned_mat) <- colnames(mat)
  
  # Copy data for genes present in this layer
  common_genes <- intersect(all_genes, rownames(mat))
  aligned_mat[common_genes, ] <- mat[common_genes, ]
  
  message("  Aligned layer ", i, ": ", nrow(aligned_mat), " x ", ncol(aligned_mat))
  
  return(aligned_mat)
})

# Combine aligned matrices column-wise
message("\nCombining aligned matrices...")
all_counts <- do.call(cbind, aligned_counts)
message("Combined counts matrix: ", nrow(all_counts), " genes x ", ncol(all_counts), " cells")

# Verify cell count matches
if (ncol(all_counts) != ncol(combined)) {
  warning("Cell count mismatch between combined matrix and object")
  
  # Get common cells and reorder
  obj_cells <- colnames(combined)
  mat_cells <- colnames(all_counts)
  common_cells <- intersect(obj_cells, mat_cells)
  
  message("  Using ", length(common_cells), " common cells")
  all_counts <- all_counts[, common_cells]
  combined <- subset(combined, cells = common_cells)
}

# Ensure matrix column order matches object
all_counts <- all_counts[, colnames(combined)]

# Create new unified Spatial assay
message("\nCreating unified Spatial assay...")
combined[["Spatial"]] <- CreateAssay5Object(counts = all_counts)

# Verify metadata preserved
if (!"sex" %in% colnames(combined@meta.data)) {
  stop("ERROR: Metadata lost during assay recreation!")
}

# Normalize data using log-normalization
# scale.factor = 10000: Counts per 10,000 (CP10K) normalization
message("Normalizing data...")
combined <- NormalizeData(
  combined,
  assay = "Spatial",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

# Verify final layer structure
final_layers <- Layers(combined[["Spatial"]])
message("\nFinal layers: ", paste(final_layers, collapse = ", "))

# Final verification
message("\nFinal verification:")
message("  - Total cells: ", ncol(combined))
message("  - Cells in counts layer: ", ncol(LayerData(combined[["Spatial"]], layer = "counts")))
message("  - Cells in data layer: ", ncol(LayerData(combined[["Spatial"]], layer = "data")))
message("  - Total genes: ", nrow(combined))

# Verify no duplicate cell names
if (any(duplicated(colnames(combined)))) {
  stop("ERROR: Duplicate cell names detected!")
}

message("\nLayers successfully unified")

# Clean up
rm(counts_list, aligned_counts, all_counts)
gc()

# =============================================================================
# STEP 4: PREPROCESSING FOR INTEGRATION
# =============================================================================

message("\n", rep("=", 80))
message("STEP 4: PREPROCESSING")
message(rep("=", 80))

# Find highly variable features
# VST: Variance stabilizing transformation
# nfeatures = 3000: Select top 3000 variable genes for integration
message("Finding variable features...")
combined <- FindVariableFeatures(
  combined,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)
message("  - Variable features selected: ", length(VariableFeatures(combined)))

# Scale data
# Centers and scales expression values for PCA
message("Scaling data...")
combined <- ScaleData(
  combined,
  features = VariableFeatures(combined),
  verbose = FALSE
)

# Run PCA
# npcs = 50: Compute 50 principal components
message("Running PCA...")
combined <- RunPCA(
  combined,
  features = VariableFeatures(combined),
  npcs = 50,
  verbose = FALSE
)

# =============================================================================
# STEP 5: HARMONY INTEGRATION
# =============================================================================

message("\n", rep("=", 80))
message("STEP 5: HARMONY BATCH CORRECTION")
message(rep("=", 80))

# Determine optimal number of PCA dimensions using elbow method
# Kneedle algorithm finds the "elbow" point in variance explained curve
kneedle <- function(y) {
  x <- seq_along(y)
  p1 <- c(1, y[1])
  p2 <- c(length(y), y[length(y)])
  
  # Calculate perpendicular distance from each point to line p1-p2
  dist <- abs((p2[2] - p1[2]) * x - (p2[1] - p1[1]) * y + 
                p2[1] * p1[2] - p2[2] * p1[1]) /
    sqrt((p2[2] - p1[2])^2 + (p2[1] - p1[1])^2)
  
  return(which.max(dist))
}

# Get standard deviations of principal components
pc_sdev <- combined[["pca"]]@stdev
elbow_pc <- kneedle(pc_sdev)

# Use dimensions up to elbow + buffer
dims_use <- seq_len(min(50, elbow_pc + 5))
message("Elbow detected at PC", elbow_pc)
message("Using PCA dimensions: 1-", max(dims_use))

# Run Harmony integration
# Corrects for batch effects between samples while preserving biological variation
message("\nRunning Harmony integration...")
combined <- RunHarmony(
  object = combined,
  group.by.vars = "sample",  # Variable to integrate on
  dims.use = dims_use,
  assay.use = "Spatial",
  theta = harmony_theta,     # Diversity penalty
  lambda = harmony_lambda,   # Ridge regression penalty
  verbose = TRUE
)

# Filter Harmony dimensions by variance
# Remove low-variance dimensions that may be noise
harmony_embeddings <- Embeddings(combined, "harmony")
harmony_sdev <- apply(harmony_embeddings, 2, sd)
variance_threshold <- quantile(harmony_sdev, 0.20)
dims_final <- which(harmony_sdev > variance_threshold)
message("\nUsing ", length(dims_final), " Harmony dimensions after variance filtering")

# =============================================================================
# STEP 6: CLUSTERING
# =============================================================================

message("\n", rep("=", 80))
message("STEP 6: CLUSTERING")
message(rep("=", 80))

# Build shared nearest neighbor graph
message("Building SNN graph...")
combined <- FindNeighbors(
  combined,
  reduction = "harmony",
  dims = dims_final,
  verbose = FALSE
)

# Cluster cells using Louvain algorithm
message("Finding clusters (resolution = ", cluster_resolution, ")...")
combined <- FindClusters(
  combined,
  resolution = cluster_resolution,
  verbose = FALSE
)

# Run UMAP for visualization
message("Running UMAP...")
combined <- RunUMAP(
  combined,
  reduction = "harmony",
  dims = dims_final,
  verbose = FALSE
)

# Report clustering results
n_clusters <- length(unique(Idents(combined)))
message("\nClusters identified: ", n_clusters)
message("\nCluster sizes:")
print(sort(table(Idents(combined)), decreasing = TRUE))

# =============================================================================
# VERIFICATION
# =============================================================================

message("\n", rep("=", 80))
message("FINAL VERIFICATION")
message(rep("=", 80))

message("\nObject summary:")
message("  - Total cells: ", ncol(combined))
message("  - Total features: ", nrow(combined))
message("  - Assays: ", paste(Assays(combined), collapse = ", "))
message("  - Reductions: ", paste(names(combined@reductions), collapse = ", "))
message("  - Clusters: ", n_clusters)

message("\nCells per sex:")
print(table(combined$sex))

message("\nCells per cluster per sex:")
print(table(Idents(combined), combined$sex))

# =============================================================================
# SAVE INTEGRATED OBJECT
# =============================================================================

message("\n", rep("=", 80))
message("SAVING INTEGRATED OBJECT")
message(rep("=", 80))

saveRDS(combined, file.path(output_dir, "combined_integrated_harmony.rds"))
message("Saved: combined_integrated_harmony.rds")

message("\n", rep("=", 80))
message("INTEGRATION COMPLETE")
message(rep("=", 80))