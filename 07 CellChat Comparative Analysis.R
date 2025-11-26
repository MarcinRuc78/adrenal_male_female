# =============================================================================
# Visium HD Spatial Transcriptomics Analysis Pipeline
# Script 07: CellChat Comparative Analysis (Male vs Female)
# =============================================================================
#
# Description:
#   This script performs comparative analysis of cell-cell communication
#   between male and female mouse adrenal glands using CellChat.
#   It identifies sex-specific signaling patterns and generates
#   publication-quality visualizations.
#
# Input:
#   - CellChat objects from Script 06 (male and female)
#
# Output:
#   - Comparative communication analysis
#   - Differential signaling heatmaps
#   - Network visualizations
#   - Information flow comparisons
#   - Pathway-level differential analysis
#
# Dependencies:
#   - CellChat (>= 2.1.0)
#   - ggplot2
#   - patchwork
#   - ComplexHeatmap
#
# References:
#   - Jin S, et al. (2021). CellChat: Inference and analysis of cell-cell
#     communication. Nature Communications.
#
# =============================================================================

# -----------------------------------------------------------------------------
# Load Required Packages
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
})

set.seed(42)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Directories
output_dir <- "results"
cellchat_dir <- file.path(output_dir, "cellchat")
comparison_dir <- file.path(cellchat_dir, "comparison")
dir.create(comparison_dir, showWarnings = FALSE, recursive = TRUE)

# Visualization parameters
plot_params <- list(
  width = 10,
  height = 8,
  dpi = 300
)

# Color schemes
sex_colors <- c("Female" = "#E41A1C", "Male" = "#377EB8")

# =============================================================================
# LOAD CELLCHAT OBJECTS
# =============================================================================

message("\n", rep("=", 80))
message("LOADING CELLCHAT OBJECTS")
message(rep("=", 80))

# Load CellChat objects
cellchat_female <- readRDS(file.path(cellchat_dir, "cellchat_female_final.rds"))
cellchat_male <- readRDS(file.path(cellchat_dir, "cellchat_male_final.rds"))

message("Female CellChat:")
message("  Cells: ", nrow(cellchat_female@meta))
message("  Cell types: ", length(unique(cellchat_female@meta$labels)))
message("  Pathways: ", length(cellchat_female@netP))

message("\nMale CellChat:")
message("  Cells: ", nrow(cellchat_male@meta))
message("  Cell types: ", length(unique(cellchat_male@meta$labels)))
message("  Pathways: ", length(cellchat_male@netP))

# =============================================================================
# LIFT CELLCHAT OBJECTS FOR COMPARISON
# =============================================================================

message("\n", rep("=", 80))
message("PREPARING OBJECTS FOR COMPARISON")
message(rep("=", 80))

# Lift CellChat objects to have the same cell types
# This is necessary when comparing datasets with potentially different cell types
group_new <- union(
  levels(cellchat_female@idents),
  levels(cellchat_male@idents)
)

message("Unified cell types: ", paste(group_new, collapse = ", "))

cellchat_female <- liftCellChat(cellchat_female, group.new = group_new)
cellchat_male <- liftCellChat(cellchat_male, group.new = group_new)

message("CellChat objects lifted to common cell type space")

# =============================================================================
# MERGE CELLCHAT OBJECTS
# =============================================================================

message("\n", rep("=", 80))
message("MERGING CELLCHAT OBJECTS")
message(rep("=", 80))

# Create merged CellChat object for comparison
cellchat_merged <- mergeCellChat(
  list(Female = cellchat_female, Male = cellchat_male),
  add.names = c("Female", "Male"),
  cell.prefix = TRUE
)

message("Merged CellChat object created")

# Save merged object
saveRDS(cellchat_merged, file.path(comparison_dir, "cellchat_merged.rds"))

# =============================================================================
# COMPARE TOTAL INTERACTION NUMBER AND STRENGTH
# =============================================================================

message("\n", rep("=", 80))
message("COMPARING TOTAL INTERACTIONS")
message(rep("=", 80))

# Compare number of interactions and interaction strength
pdf(file.path(comparison_dir, "interaction_comparison_barplot.pdf"),
    width = 8, height = 5)

gg_interaction <- compareInteractions(
  cellchat_merged,
  show.legend = TRUE,
  group = c(1, 2),
  measure = c("count", "weight")
)
print(gg_interaction)

dev.off()
message("Saved: interaction_comparison_barplot.pdf")

# Print summary statistics
message("\nInteraction summary:")
message("  Female - Count: ", sum(cellchat_female@net$count > 0),
        ", Weight: ", round(sum(cellchat_female@net$weight), 4))
message("  Male - Count: ", sum(cellchat_male@net$count > 0),
        ", Weight: ", round(sum(cellchat_male@net$weight), 4))

# =============================================================================
# COMPARE NUMBER OF INTERACTIONS BY CELL TYPE PAIR
# =============================================================================

message("\n", rep("=", 80))
message("DIFFERENTIAL NUMBER OF INTERACTIONS")
message(rep("=", 80))

# Heatmap showing differential number of interactions
pdf(file.path(comparison_dir, "diff_interactions_heatmap.pdf"),
    width = 12, height = 10)

# Number of interactions
netVisual_heatmap(
  cellchat_merged,
  measure = "count",
  comparison = c(1, 2),
  title.name = "Differential number of interactions"
)

dev.off()

# Interaction strength
pdf(file.path(comparison_dir, "diff_strength_heatmap.pdf"),
    width = 12, height = 10)

netVisual_heatmap(
  cellchat_merged,
  measure = "weight",
  comparison = c(1, 2),
  title.name = "Differential interaction strength"
)

dev.off()

message("Saved: diff_interactions_heatmap.pdf")
message("Saved: diff_strength_heatmap.pdf")

# =============================================================================
# COMPARE OUTGOING AND INCOMING SIGNALING
# =============================================================================

message("\n", rep("=", 80))
message("COMPARING SIGNALING PATTERNS")
message(rep("=", 80))

# Compare outgoing signaling patterns
pdf(file.path(comparison_dir, "outgoing_signaling_comparison.pdf"),
    width = 14, height = 8)

# Identify signaling groups based on functional similarity
cellchat_merged <- computeNetSimilarityPairwise(
  cellchat_merged,
  type = "functional"
)

cellchat_merged <- netEmbedding(
  cellchat_merged,
  type = "functional",
  umap.method = "umap"
)

# Clustering of signaling pathways
cellchat_merged <- netClustering(
  cellchat_merged,
  type = "functional",
  do.parallel = FALSE
)

# Visualization
netVisual_embeddingPairwise(
  cellchat_merged,
  type = "functional",
  label.size = 3.5
)

dev.off()
message("Saved: outgoing_signaling_comparison.pdf")

# =============================================================================
# IDENTIFY DIFFERENTIAL SIGNALING PATHWAYS
# =============================================================================

message("\n", rep("=", 80))
message("IDENTIFYING DIFFERENTIAL PATHWAYS")
message(rep("=", 80))

# Rank signaling pathways by differences
pdf(file.path(comparison_dir, "pathway_ranking.pdf"),
    width = 10, height = 12)

# Overall information flow for each signaling pathway
gg_info_flow <- rankNet(
  cellchat_merged,
  mode = "comparison",
  stacked = TRUE,
  do.stat = TRUE
)
print(gg_info_flow)

dev.off()
message("Saved: pathway_ranking.pdf")

# Detailed pathway comparison
pdf(file.path(comparison_dir, "pathway_ranking_unstacked.pdf"),
    width = 10, height = 12)

gg_info_flow_unstacked <- rankNet(
  cellchat_merged,
  mode = "comparison",
  stacked = FALSE,
  do.stat = TRUE
)
print(gg_info_flow_unstacked)

dev.off()
message("Saved: pathway_ranking_unstacked.pdf")

# =============================================================================
# COMPARE SIGNALING ASSOCIATED WITH CELL POPULATIONS
# =============================================================================

message("\n", rep("=", 80))
message("CELL TYPE-SPECIFIC SIGNALING")
message(rep("=", 80))

# Outgoing signaling patterns (each cell type as signal sender)
pdf(file.path(comparison_dir, "outgoing_signaling_by_celltype.pdf"),
    width = 16, height = 10)

# Female
gg_out_f <- netAnalysis_signalingRole_heatmap(
  cellchat_female,
  pattern = "outgoing",
  signaling = NULL,
  title = "Outgoing signaling - Female",
  width = 10,
  height = 14
)

# Male
gg_out_m <- netAnalysis_signalingRole_heatmap(
  cellchat_male,
  pattern = "outgoing",
  signaling = NULL,
  title = "Outgoing signaling - Male",
  width = 10,
  height = 14
)

# Combine using draw from ComplexHeatmap
ht_list <- gg_out_f + gg_out_m
draw(ht_list, ht_gap = unit(0.5, "cm"))

dev.off()
message("Saved: outgoing_signaling_by_celltype.pdf")

# Incoming signaling patterns (each cell type as signal receiver)
pdf(file.path(comparison_dir, "incoming_signaling_by_celltype.pdf"),
    width = 16, height = 10)

gg_in_f <- netAnalysis_signalingRole_heatmap(
  cellchat_female,
  pattern = "incoming",
  signaling = NULL,
  title = "Incoming signaling - Female",
  width = 10,
  height = 14
)

gg_in_m <- netAnalysis_signalingRole_heatmap(
  cellchat_male,
  pattern = "incoming",
  signaling = NULL,
  title = "Incoming signaling - Male",
  width = 10,
  height = 14
)

ht_list_in <- gg_in_f + gg_in_m
draw(ht_list_in, ht_gap = unit(0.5, "cm"))

dev.off()
message("Saved: incoming_signaling_by_celltype.pdf")

# =============================================================================
# IDENTIFY CONSERVED AND CONTEXT-SPECIFIC SIGNALING
# =============================================================================

message("\n", rep("=", 80))
message("CONSERVED VS SEX-SPECIFIC SIGNALING")
message(rep("=", 80))

# Identify pathways present in both conditions
pathways_female <- names(cellchat_female@netP)
pathways_male <- names(cellchat_male@netP)

pathways_common <- intersect(pathways_female, pathways_male)
pathways_female_only <- setdiff(pathways_female, pathways_male)
pathways_male_only <- setdiff(pathways_male, pathways_female)

message("\nPathway comparison:")
message("  Common pathways: ", length(pathways_common))
message("  Female-specific: ", length(pathways_female_only))
message("  Male-specific: ", length(pathways_male_only))

if (length(pathways_female_only) > 0) {
  message("  Female-only: ", paste(pathways_female_only, collapse = ", "))
}
if (length(pathways_male_only) > 0) {
  message("  Male-only: ", paste(pathways_male_only, collapse = ", "))
}

# Save pathway lists
pathway_comparison <- list(
  common = pathways_common,
  female_specific = pathways_female_only,
  male_specific = pathways_male_only
)
saveRDS(pathway_comparison, file.path(comparison_dir, "pathway_comparison.rds"))

# =============================================================================
# VISUALIZE SPECIFIC PATHWAYS OF INTEREST
# =============================================================================

message("\n", rep("=", 80))
message("PATHWAY-SPECIFIC VISUALIZATIONS")
message(rep("=", 80))

# Define pathways of interest for adrenal gland
pathways_of_interest <- c("WNT", "HH", "BMP", "NOTCH", "IGF", "VEGF", "PDGF")
pathways_to_plot <- intersect(pathways_of_interest, pathways_common)

if (length(pathways_to_plot) > 0) {
  
  for (pathway in pathways_to_plot) {
    message("  Visualizing pathway: ", pathway)
    
    # Circle plot comparison
    pdf(file.path(comparison_dir, paste0("pathway_", pathway, "_circle.pdf")),
        width = 14, height = 6)
    
    par(mfrow = c(1, 2))
    
    # Female
    netVisual_aggregate(
      cellchat_female,
      signaling = pathway,
      layout = "circle",
      title.name = paste0(pathway, " - Female")
    )
    
    # Male
    netVisual_aggregate(
      cellchat_male,
      signaling = pathway,
      layout = "circle",
      title.name = paste0(pathway, " - Male")
    )
    
    dev.off()
    
    # Chord diagram comparison
    pdf(file.path(comparison_dir, paste0("pathway_", pathway, "_chord.pdf")),
        width = 14, height = 6)
    
    par(mfrow = c(1, 2))
    
    netVisual_aggregate(
      cellchat_female,
      signaling = pathway,
      layout = "chord",
      title.name = paste0(pathway, " - Female")
    )
    
    netVisual_aggregate(
      cellchat_male,
      signaling = pathway,
      layout = "chord",
      title.name = paste0(pathway, " - Male")
    )
    
    dev.off()
    
    # Heatmap
    pdf(file.path(comparison_dir, paste0("pathway_", pathway, "_heatmap.pdf")),
        width = 12, height = 5)
    
    par(mfrow = c(1, 2))
    
    netVisual_heatmap(
      cellchat_female,
      signaling = pathway,
      title.name = paste0(pathway, " - Female"),
      color.heatmap = "Reds"
    )
    
    netVisual_heatmap(
      cellchat_male,
      signaling = pathway,
      title.name = paste0(pathway, " - Male"),
      color.heatmap = "Blues"
    )
    
    dev.off()
  }
  
  message("Saved pathway-specific visualizations")
}

# =============================================================================
# LIGAND-RECEPTOR PAIR ANALYSIS
# =============================================================================

message("\n", rep("=", 80))
message("LIGAND-RECEPTOR PAIR ANALYSIS")
message(rep("=", 80))

# Compare specific L-R pairs between conditions
pdf(file.path(comparison_dir, "LR_pair_comparison_bubble.pdf"),
    width = 14, height = 16)

# Bubble plot for all pathways
gg_bubble <- netVisual_bubble(
  cellchat_merged,
  sources.use = NULL,
  targets.use = NULL,
  comparison = c(1, 2),
  angle.x = 45,
  remove.isolate = TRUE,
  max.dataset = 2,
  title.name = "L-R pairs: Female vs Male"
)
print(gg_bubble)

dev.off()
message("Saved: LR_pair_comparison_bubble.pdf")

# =============================================================================
# EXPORT COMMUNICATION PROBABILITY TABLES
# =============================================================================

message("\n", rep("=", 80))
message("EXPORTING COMMUNICATION TABLES")
message(rep("=", 80))

# Extract communication data frames
df_female <- subsetCommunication(cellchat_female)
df_male <- subsetCommunication(cellchat_male)

# Add sex label
df_female$sex <- "Female"
df_male$sex <- "Male"

# Combine
df_combined <- rbind(df_female, df_male)

# Save
write.csv(df_female, 
          file.path(comparison_dir, "communication_female.csv"),
          row.names = FALSE)
write.csv(df_male, 
          file.path(comparison_dir, "communication_male.csv"),
          row.names = FALSE)
write.csv(df_combined, 
          file.path(comparison_dir, "communication_combined.csv"),
          row.names = FALSE)

message("Saved communication tables")

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

message("\n", rep("=", 80))
message("GENERATING SUMMARY STATISTICS")
message(rep("=", 80))

# Calculate pathway-level statistics
pathway_stats <- data.frame(
  pathway = pathways_common,
  female_interactions = sapply(pathways_common, function(p) {
    sum(cellchat_female@netP[[p]]$prob > 0, na.rm = TRUE)
  }),
  male_interactions = sapply(pathways_common, function(p) {
    sum(cellchat_male@netP[[p]]$prob > 0, na.rm = TRUE)
  }),
  female_weight = sapply(pathways_common, function(p) {
    sum(cellchat_female@netP[[p]]$prob, na.rm = TRUE)
  }),
  male_weight = sapply(pathways_common, function(p) {
    sum(cellchat_male@netP[[p]]$prob, na.rm = TRUE)
  })
)

pathway_stats$diff_interactions <- pathway_stats$male_interactions - 
  pathway_stats$female_interactions
pathway_stats$diff_weight <- pathway_stats$male_weight - pathway_stats$female_weight
pathway_stats$log2FC_weight <- log2((pathway_stats$male_weight + 0.001) / 
                                      (pathway_stats$female_weight + 0.001))

# Sort by absolute difference
pathway_stats <- pathway_stats[order(-abs(pathway_stats$log2FC_weight)), ]

# Save
write.csv(pathway_stats,
          file.path(comparison_dir, "pathway_statistics.csv"),
          row.names = FALSE)

message("\nTop 10 differentially active pathways (by weight):")
print(head(pathway_stats[, c("pathway", "female_weight", "male_weight", "log2FC_weight")], 10))

# =============================================================================
# ANALYSIS COMPLETE
# =============================================================================

message("\n", rep("=", 80))
message("COMPARATIVE ANALYSIS COMPLETED")
message(rep("=", 80))

message("\nOutput files saved to: ", comparison_dir)
message("\nGenerated files:")
message("  - cellchat_merged.rds (merged CellChat object)")
message("  - interaction_comparison_barplot.pdf")
message("  - diff_interactions_heatmap.pdf")
message("  - diff_strength_heatmap.pdf")
message("  - pathway_ranking.pdf")
message("  - outgoing/incoming_signaling_by_celltype.pdf")
message("  - pathway_comparison.rds")
message("  - Pathway-specific visualizations")
message("  - communication_*.csv (communication tables)")
message("  - pathway_statistics.csv")

message("\nKey findings:")
message("  Common pathways: ", length(pathways_common))
message("  Female-specific pathways: ", length(pathways_female_only))
message("  Male-specific pathways: ", length(pathways_male_only))