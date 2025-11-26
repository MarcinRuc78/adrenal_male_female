# =============================================================================
# Visium HD Spatial Transcriptomics Analysis Pipeline
# Script 03: Cluster Annotation and Visualization
# =============================================================================
#
# Description:
#   This script performs cluster filtering, annotation with biologically
#   meaningful cell type names, and generates publication-quality
#   visualizations including UMAP plots and spatial maps.
#
# Input:
#   - Integrated Seurat object from Script 02
#
# Output:
#   - Annotated Seurat object with cell type identities
#   - UMAP visualizations (clusters, sex)
#   - Spatial dimension plots
#   - PDF and PNG figures for publication
#
# Dependencies:
#   - Seurat (>= 5.0)
#   - ggplot2
#   - patchwork
#   - viridis
#   - RColorBrewer
#
# =============================================================================

# -----------------------------------------------------------------------------
# Load Required Packages
# -----------------------------------------------------------------------------

library(Seurat)
library(ggplot2)
library(patchwork)
library(viridis)
library(RColorBrewer)
library(scales)

set.seed(1234)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Output directory
output_dir <- "results"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Color palette for adrenal gland cell populations
# Colors are carefully chosen to distinguish adjacent zones
cell_type_colors <- c(
  "outer ZF"        = "#E41A1C",   # Red - outer zona fasciculata
  
  "ZG"              = "#377EB8",   # Blue - zona glomerulosa
  "inner ZF"        = "#4DAF4A",   # Green - inner zona fasciculata
  "X/JMZ"           = "#984EA3",   # Purple - X-zone/juxtamedullary zone
  "BAT"             = "#FF7F00",   # Orange - brown adipose tissue
  "CT capsule"      = "#A65628",   # Brown - connective tissue capsule
  "WAT"             = "#F781BF",   # Pink - white adipose tissue
  "medulla"         = "#00CED1",   # Cyan - adrenal medulla
  "Subcapsular ZG"  = "#FFD92F",   # Yellow - subcapsular zona glomerulosa
  "fibroblasts"     = "grey"       # Grey - fibroblasts
)

# Sex colors
sex_colors <- c("Male" = "#377EB8", "Female" = "#E41A1C")

# Clusters to exclude from analysis (e.g., low quality, artifacts)
clusters_to_remove <- c("10", "11")

# Publication-quality theme for plots
theme_publication <- theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right"
  )

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

#' Save Plot as PDF
#'
#' Saves a ggplot object as a high-quality PDF using cairo_pdf device
#' for proper font embedding and transparency support.
#'
#' @param plot ggplot object
#' @param filename Output filename
#' @param width Plot width in inches
#' @param height Plot height in inches
save_pdf <- function(plot, filename, width = 10, height = 5) {
  ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    device = cairo_pdf
  )
}

# -----------------------------------------------------------------------------
# Load Data
# -----------------------------------------------------------------------------

message("\n", rep("=", 80))
message("LOADING INTEGRATED DATA")
message(rep("=", 80))

# Load integrated object (from Script 02)
combined <- readRDS(file.path(output_dir, "combined_integrated_harmony.rds"))

message("Loaded object with ", ncol(combined), " cells")

# =============================================================================
# STEP 7: CLUSTER FILTERING
# =============================================================================

message("\n", rep("=", 80))
message("STEP 7: CLUSTER FILTERING")
message(rep("=", 80))

# Check cluster sizes before filtering
cluster_sizes <- table(Idents(combined))
message("Cluster sizes before filtering:")
print(sort(cluster_sizes, decreasing = TRUE))

# Identify clusters present in data that we want to remove
clusters_present <- intersect(clusters_to_remove, names(cluster_sizes))

if (length(clusters_present) > 0) {
  message("\nRemoving clusters: ", paste(clusters_present, collapse = ", "))
  
  # Get cells to keep (all clusters except those being removed)
  clusters_keep <- setdiff(levels(Idents(combined)), clusters_present)
  cells_keep <- WhichCells(combined, idents = clusters_keep)
  
  # Subset object
  combined <- subset(combined, cells = cells_keep)
  
  # Reset factor levels to remove empty levels
  Idents(combined) <- factor(as.character(Idents(combined)))
  
  message("Cells remaining after filtering: ", ncol(combined))
} else {
  message("No clusters to remove")
}

# =============================================================================
# STEP 8: CLUSTER ANNOTATION
# =============================================================================

message("\n", rep("=", 80))
message("STEP 8: CLUSTER ANNOTATION")
message(rep("=", 80))

# Map cluster numbers to biologically meaningful cell type names
# These annotations are based on marker gene expression analysis
# and known adrenal gland anatomy

cluster_annotations <- c(
  "0" = "inner ZF",        # Inner zona fasciculata
  "1" = "outer ZF",        # Outer zona fasciculata
  "2" = "X/JMZ",           # X-zone / juxtamedullary zone
  "3" = "ZG",              # Zona glomerulosa
  "4" = "medulla",         # Adrenal medulla (chromaffin cells)
  "5" = "CT capsule",      # Connective tissue capsule
  "6" = "WAT",             # White adipose tissue
  "7" = "BAT",             # Brown adipose tissue
  "8" = "Subcapsular ZG",  # Subcapsular zona glomerulosa
  "9" = "fibroblasts"      # Fibroblasts
)

# Apply annotations
combined <- RenameIdents(combined, cluster_annotations)

# Store cell type in metadata for easy access
combined$cell_type <- Idents(combined)

# Define cluster order (from capsule inward, then medulla and adipose)
cluster_order <- c(
  "fibroblasts", "CT capsule", "Subcapsular ZG", "ZG",
  "outer ZF", "inner ZF", "X/JMZ",
  "medulla", "WAT", "BAT"
)

# Reorder factor levels
combined$cell_type <- factor(
  combined$cell_type,
  levels = intersect(cluster_order, unique(combined$cell_type))
)
Idents(combined) <- combined$cell_type

# Print final cluster distribution
message("\nFinal cluster distribution:")
print(table(Idents(combined)))

message("\nCluster distribution by sex:")
print(table(Idents(combined), combined$sex))

# =============================================================================
# STEP 9: VISUALIZATION
# =============================================================================

message("\n", rep("=", 80))
message("STEP 9: CREATING VISUALIZATIONS")
message(rep("=", 80))

# -------------------------------------------------------------------------
# UMAP: Cell populations
# -------------------------------------------------------------------------
message("Creating UMAP plot - cell populations...")

p_umap_clusters <- DimPlot(
  combined,
  reduction = "umap",
  label = TRUE,
  label.size = 4,
  cols = cell_type_colors,
  pt.size = 0.5
) +
  ggtitle("UMAP - Adrenal Gland Cell Populations") +
  theme_publication +
  NoLegend()

save_pdf(p_umap_clusters, file.path(output_dir, "umap_clusters.pdf"), 
         width = 6, height = 5)
message("  Saved: umap_clusters.pdf")

# -------------------------------------------------------------------------
# UMAP: By sex
# -------------------------------------------------------------------------
message("Creating UMAP plot - by sex...")

p_umap_sex <- DimPlot(
  combined,
  reduction = "umap",
  group.by = "sex",
  cols = sex_colors,
  pt.size = 0.5
) +
  ggtitle("UMAP - By Sex") +
  theme_publication

save_pdf(p_umap_sex, file.path(output_dir, "umap_by_sex.pdf"), 
         width = 6, height = 5)
message("  Saved: umap_by_sex.pdf")

# -------------------------------------------------------------------------
# UMAP: Split by sex
# -------------------------------------------------------------------------
message("Creating UMAP plot - split by sex...")

p_umap_split <- DimPlot(
  combined,
  reduction = "umap",
  split.by = "sex",
  cols = cell_type_colors,
  label = TRUE,
  label.size = 3,
  pt.size = 0.3
) +
  ggtitle("UMAP - Cell Populations Split by Sex") +
  theme_publication

save_pdf(p_umap_split, file.path(output_dir, "umap_split_by_sex.pdf"), 
         width = 12, height = 5)
message("  Saved: umap_split_by_sex.pdf")

# -------------------------------------------------------------------------
# Combined UMAP panel
# -------------------------------------------------------------------------
message("Creating combined UMAP panel...")

p_umap_combined <- (p_umap_clusters | p_umap_sex) +
  plot_annotation(
    title = "Visium HD Mouse Adrenal Gland",
    theme = theme(plot.title = element_text(face = "bold", size = 16))
  )

save_pdf(p_umap_combined, file.path(output_dir, "umap_combined_panel.pdf"), 
         width = 14, height = 6)
message("  Saved: umap_combined_panel.pdf")

# -------------------------------------------------------------------------
# Spatial dimension plot
# -------------------------------------------------------------------------
message("Creating spatial dimension plot...")

p_spatial <- SpatialDimPlot(
  combined,
  images = Images(combined),
  cols = cell_type_colors,
  label = TRUE,
  label.size = 6,
  pt.size.factor = 3,
  image.alpha = 1,
  crop = TRUE,
  combine = TRUE
) +
  ggtitle("Spatial Distribution - Adrenal Cell Populations") +
  theme(plot.title = element_text(face = "bold", size = 14))

save_pdf(p_spatial, file.path(output_dir, "spatial_clusters.pdf"), 
         width = 18, height = 16)
message("  Saved: spatial_clusters.pdf")

# -------------------------------------------------------------------------
# Cell type composition bar plot
# -------------------------------------------------------------------------
message("Creating cell type composition plot...")

# Calculate proportions
composition_data <- combined@meta.data %>%
  group_by(sex, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(sex) %>%
  mutate(
    total = sum(n),
    proportion = n / total * 100
  )

p_composition <- ggplot(
  composition_data,
  aes(x = sex, y = proportion, fill = cell_type)
) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = cell_type_colors, name = "Cell Type") +
  labs(
    x = "Sex",
    y = "Proportion (%)",
    title = "Cell Type Composition by Sex"
  ) +
  theme_publication +
  theme(legend.position = "right")

save_pdf(p_composition, file.path(output_dir, "cell_type_composition.pdf"), 
         width = 8, height = 6)
message("  Saved: cell_type_composition.pdf")

# =============================================================================
# SAVE ANNOTATED OBJECT
# =============================================================================

message("\n", rep("=", 80))
message("SAVING ANNOTATED OBJECT")
message(rep("=", 80))

saveRDS(combined, file.path(output_dir, "combined_annotated.rds"))
message("Saved: combined_annotated.rds")

# =============================================================================
# SUMMARY
# =============================================================================

message("\n", rep("=", 80))
message("ANNOTATION AND VISUALIZATION COMPLETE")
message(rep("=", 80))

message("\nGenerated files:")
message("  - combined_annotated.rds (annotated Seurat object)")
message("  - umap_clusters.pdf")
message("  - umap_by_sex.pdf")
message("  - umap_split_by_sex.pdf")
message("  - umap_combined_panel.pdf")
message("  - spatial_clusters.pdf")
message("  - cell_type_composition.pdf")

message("\nNext steps:")
message("  1. Find cluster marker genes (Script 04)")
message("  2. Perform sex-differential expression analysis (Script 05)")
message("  3. Conduct pathway enrichment analysis")