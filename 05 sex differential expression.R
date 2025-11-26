# =============================================================================
# Visium HD Spatial Transcriptomics Analysis Pipeline
# Script 05: Sex-Differential Expression Analysis
# =============================================================================
#
# Description:
#   This script performs differential expression analysis between male and
#   female samples, both whole-adrenal and per-cluster. Uses pseudobulk
#   approach with edgeR for robust statistical testing.
#
# Input:
#   - Annotated Seurat object from Script 03
#
# Output:
#   - Excel files with DE results (whole adrenal and per-cluster)
#   - Volcano plots for each comparison
#
# Methods:
#   - Pseudobulk aggregation per sample/replicate
#   - edgeR quasi-likelihood F-test (QLF) for differential expression
#   - Multiple testing correction using Benjamini-Hochberg FDR
#   - Gene filtering based on minimum expression percentage and specificity
#
# Statistical approach:
#   - When replicates are available (≥2 per sex): edgeR QLF with estimated dispersion
#   - When replicates are limited: fallback to fixed dispersion (0.1)
#
# Dependencies:
#   - Seurat (>= 5.0)
#   - edgeR
#   - dplyr
#   - tibble
#   - openxlsx
#   - EnhancedVolcano
#   - ggplot2
#
# =============================================================================

# -----------------------------------------------------------------------------
# Load Required Packages
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(dplyr)
  library(tibble)
  library(openxlsx)
  library(EnhancedVolcano)
  library(ggplot2)
})

set.seed(1234)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Output directory
output_dir <- "results"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Create subdirectory for plots
plots_dir <- file.path(output_dir, "volcano_plots_sex_DE")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# Analysis parameters
params <- list(
  # Gene filtering
  min_pct_within_cluster = 0.20,  # Min % cells expressing gene in cluster
  specificity_cutoff = 0.20,       # Min specificity score for cluster
  
  # Marker detection (for gene filtering)
  markers_min_pct = 0.20,
  markers_logfc_threshold = 0.20,
  markers_padj_threshold = 0.05,
  markers_min_cells_per_ident = 20,
  markers_max_cells_per_ident = 4000,
  markers_hvg_n = 16000,
  
  # edgeR parameters
  disp_fallback = 0.1,   # Fixed dispersion when replicates insufficient
  
  # Plotting thresholds
  p_adj_threshold = 0.05,
  lfc_threshold = 0.25,
  
  # Plot dimensions
  plot_width = 6,
  plot_height = 5
)

# Color palette (colorblind-friendly Okabe-Ito)
okabe_ito <- c(
  "#0072B2", "#D55E00", "#009E73", "#CC79A7",
  "#F0E442", "#56B4E9", "#E69F00", "#000000"
)

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

#' Run edgeR Differential Expression (Two-Group Comparison)
#'
#' Performs differential expression analysis using edgeR with either
#' estimated or fixed dispersion.
#'
#' @param counts_mat Matrix of counts (genes x samples)
#' @param group Factor specifying group membership (Female/Male)
#' @param use_disp_fallback Logical, use fixed dispersion instead of estimation
#' @param disp_fallback Fixed dispersion value when use_disp_fallback = TRUE
#' @return Data frame with DE results (gene, log2FC, p-value, FDR)
edger_run <- function(counts_mat, group, 
                      use_disp_fallback = FALSE, 
                      disp_fallback = 0.1) {
  
  # Create DGEList object
  y <- DGEList(counts = counts_mat, group = group)
  
  # Filter lowly expressed genes
  keep <- filterByExpr(y, group = group)
  if (!any(keep)) return(NULL)
  y <- y[keep, , keep.lib.sizes = FALSE]
  
  # Calculate normalization factors
  y <- calcNormFactors(y)
  
  # Design matrix: ~0 + group for coefficient-based contrast
  design <- model.matrix(~ 0 + group)
  colnames(design) <- sub("^group", "", colnames(design))
  
  if (use_disp_fallback) {
    # Use fixed dispersion (for limited replicates)
    fit <- glmFit(y, design = design, dispersion = disp_fallback)
    tt <- glmLRT(fit, contrast = c(-1, 1))  # Male - Female
  } else {
    # Estimate dispersion from data
    y <- estimateDisp(y, design = design)
    fit <- glmQLFit(y, design = design)
    tt <- glmQLFTest(fit, contrast = c(-1, 1))  # Male - Female
  }
  
  # Extract results
  results <- topTags(tt, n = Inf)$table %>%
    rownames_to_column("gene") %>%
    rename(
      avg_log2FC = logFC,
      p_val = PValue,
      p_val_adj = FDR
    ) %>%
    arrange(p_val_adj, desc(abs(avg_log2FC)))
  
  return(results)
}

#' Create and Save Volcano Plot
#'
#' Generates publication-quality volcano plot using EnhancedVolcano.
#'
#' @param df Data frame with gene, avg_log2FC, p_val_adj columns
#' @param title Plot title
#' @param outbase Output file path (without extension)
#' @param p_adj_thr Adjusted p-value threshold for significance
#' @param lfc_thr Log2 fold-change threshold for significance
#' @param caption Optional caption text
make_volcano <- function(df, title, outbase,
                         p_adj_thr = 0.05,
                         lfc_thr = 0.25,
                         caption = NULL) {
  
  if (is.null(df) || !nrow(df)) return(invisible(NULL))
  
  # Handle extreme p-values
  df$p_val_adj[!is.finite(df$p_val_adj)] <- 1
  df$p_val_adj[df$p_val_adj < .Machine$double.xmin] <- .Machine$double.xmin
  
  # Create labels
  labels <- df$gene
  names(labels) <- df$gene
  
  # Create volcano plot
  p <- EnhancedVolcano(
    df,
    lab = labels,
    x = "avg_log2FC",
    y = "p_val_adj",
    xlab = "log2 Fold Change (Male vs Female)",
    ylab = "-log10(Adjusted P-value)",
    pCutoff = p_adj_thr,
    FCcutoff = lfc_thr,
    col = c("grey80", okabe_ito[1], okabe_ito[2], okabe_ito[3]),
    drawConnectors = TRUE,
    widthConnectors = 0.4,
    labSize = 3.5,
    title = title,
    subtitle = NULL,
    legendPosition = "none"
  ) +
    theme(
      text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 12),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(face = "bold"),
      plot.caption = element_text(face = "bold", hjust = 0, size = 8)
    )
  
  # Add caption if provided
  if (!is.null(caption)) {
    p <- p + labs(caption = caption)
  }
  
  # Save plots
  ggsave(paste0(outbase, ".pdf"), p, 
         width = params$plot_width, height = params$plot_height, 
         device = cairo_pdf)
  ggsave(paste0(outbase, ".png"), p, 
         width = params$plot_width, height = params$plot_height, 
         dpi = 300)
  
  invisible(p)
}

#' Normalize Sex Labels
#'
#' Standardizes various sex label formats to "Male" or "Female".
#'
#' @param x Character vector of sex labels
#' @return Factor with levels c("Female", "Male")
normalize_sex <- function(x) {
  x0 <- tolower(trimws(as.character(x)))
  result <- ifelse(
    x0 %in% c("male", "m"),
    "Male",
    ifelse(x0 %in% c("female", "f"), "Female", NA)
  )
  factor(result, levels = c("Female", "Male"))
}

# -----------------------------------------------------------------------------
# Load Data
# -----------------------------------------------------------------------------

message("\n", rep("=", 80))
message("LOADING DATA")
message(rep("=", 80))

combined <- readRDS(file.path(output_dir, "combined_annotated.rds"))

# Verify required metadata
stopifnot("Spatial" %in% Assays(combined))
stopifnot("sex" %in% colnames(combined@meta.data))

message("Loaded object with ", ncol(combined), " cells")

DefaultAssay(combined) <- "Spatial"

# -----------------------------------------------------------------------------
# Prepare Sex Vector
# -----------------------------------------------------------------------------

# Normalize sex labels
sex_vec <- normalize_sex(combined@meta.data$sex)
names(sex_vec) <- colnames(combined)

message("\nSex distribution:")
print(table(sex_vec, useNA = "ifany"))

# Check for NA values
n_na <- sum(is.na(sex_vec))
if (n_na > 0) {
  warning("Found ", n_na, " cells with NA sex annotation")
}

# =============================================================================
# CALCULATE GENE SPECIFICITY PER CLUSTER
# =============================================================================

message("\n", rep("=", 80))
message("CALCULATING GENE SPECIFICITY")
message(rep("=", 80))

# Get layer information
layers <- Layers(combined[["Spatial"]])
counts_layers <- grep("^counts", layers, value = TRUE)
data_layers <- grep("^data", layers, value = TRUE)

message("Counts layers: ", paste(counts_layers, collapse = ", "))
message("Data layers: ", paste(data_layers, collapse = ", "))

# Get cluster assignments
clusters_all <- as.character(Idents(combined))
names(clusters_all) <- colnames(combined)
cluster_levels <- unique(clusters_all)

# Calculate mean expression per cluster per gene using data layers
if (length(data_layers) > 0) {
  message("Calculating specificity from normalized data...")
  
  # Initialize
  X_ref <- LayerData(combined[["Spatial"]], layer = data_layers[1])
  genes_ref <- rownames(X_ref)
  rm(X_ref)
  
  MEAN <- matrix(
    0, 
    nrow = length(genes_ref), 
    ncol = length(cluster_levels),
    dimnames = list(genes_ref, cluster_levels)
  )
  N_cells <- setNames(integer(length(cluster_levels)), cluster_levels)
  
  # Accumulate means across layers
  for (layer in data_layers) {
    X_data <- LayerData(combined[["Spatial"]], layer = layer)
    
    # Align genes if needed
    if (!identical(rownames(X_data), genes_ref)) {
      X_data <- X_data[genes_ref, , drop = FALSE]
    }
    
    cluster_vec <- clusters_all[colnames(X_data)]
    keep <- !is.na(cluster_vec)
    
    if (!any(keep)) {
      rm(X_data)
      next
    }
    
    X_data <- X_data[, keep, drop = FALSE]
    cluster_vec <- cluster_vec[keep]
    
    for (cl in cluster_levels) {
      idx <- which(cluster_vec == cl)
      if (!length(idx)) next
      
      MEAN[, cl] <- MEAN[, cl] + 
        as.numeric(Matrix::rowMeans(X_data[, idx, drop = FALSE])) * length(idx)
      N_cells[cl] <- N_cells[cl] + length(idx)
    }
    
    rm(X_data)
    gc()
  }
  
  # Calculate final means and specificity
  for (cl in cluster_levels) {
    if (N_cells[cl] > 0) {
      MEAN[, cl] <- MEAN[, cl] / N_cells[cl]
    }
  }
  
  # Specificity = cluster_mean / row_sum (proportion of expression in cluster)
  row_sums <- rowSums(MEAN)
  row_sums[row_sums == 0] <- 1
  SPECIFICITY <- MEAN / row_sums
  
} else {
  warning("No data layers found - using uniform specificity")
  X_tmp <- LayerData(combined[["Spatial"]], layer = counts_layers[1])
  SPECIFICITY <- matrix(
    1, 
    nrow = nrow(X_tmp), 
    ncol = length(cluster_levels),
    dimnames = list(rownames(X_tmp), cluster_levels)
  )
  rm(X_tmp)
}

message("Specificity matrix: ", nrow(SPECIFICITY), " genes x ", 
        ncol(SPECIFICITY), " clusters")

# =============================================================================
# PER-CLUSTER DIFFERENTIAL EXPRESSION
# =============================================================================

message("\n", rep("=", 80))
message("DIFFERENTIAL EXPRESSION: PER-CLUSTER (Male vs Female)")
message(rep("=", 80))

# Get reference gene list
X_ref <- LayerData(combined[["Spatial"]], layer = counts_layers[1])
genes_counts_ref <- rownames(X_ref)
rm(X_ref)

# Initialize Excel workbook
wb_clusters <- createWorkbook()

# Process each cluster
for (cluster in cluster_levels) {
  message("\n  Processing cluster: ", cluster)
  
  # Initialize accumulators
  female_sum <- NULL
  male_sum <- NULL
  nz_female_sum <- NULL
  nz_male_sum <- NULL
  n_female_total <- 0L
  n_male_total <- 0L
  
  # Accumulate counts across layers
  for (layer in counts_layers) {
    X <- LayerData(combined[["Spatial"]], layer = layer)
    
    # Align genes to reference
    if (!identical(rownames(X), genes_counts_ref)) {
      common <- intersect(genes_counts_ref, rownames(X))
      X2 <- Matrix::Matrix(
        0, 
        nrow = length(genes_counts_ref), 
        ncol = ncol(X), 
        sparse = TRUE
      )
      rownames(X2) <- genes_counts_ref
      colnames(X2) <- colnames(X)
      X2[match(common, genes_counts_ref), ] <- X[common, , drop = FALSE]
      X <- X2
      rm(X2)
    }
    
    # Get sex and cluster for cells in this layer
    sx <- sex_vec[colnames(X)]
    cl_vec <- clusters_all[colnames(X)]
    
    # Filter to current cluster
    keep <- !is.na(sx) & !is.na(cl_vec) & cl_vec == cluster
    if (!any(keep)) {
      rm(X)
      next
    }
    
    sx <- sx[keep]
    X <- X[, keep, drop = FALSE]
    
    idx_male <- which(sx == "Male")
    idx_female <- which(sx == "Female")
    
    # Accumulate female counts
    if (length(idx_female) > 0) {
      s_female <- if (length(idx_female) == 1) {
        as.numeric(X[, idx_female])
      } else {
        as.numeric(Matrix::rowSums(X[, idx_female, drop = FALSE]))
      }
      
      nz_female <- if (length(idx_female) == 1) {
        as.numeric(X[, idx_female] > 0)
      } else {
        as.numeric(Matrix::rowSums(X[, idx_female, drop = FALSE] > 0))
      }
      
      female_sum <- if (is.null(female_sum)) s_female else female_sum + s_female
      nz_female_sum <- if (is.null(nz_female_sum)) nz_female else nz_female_sum + nz_female
      n_female_total <- n_female_total + length(idx_female)
    }
    
    # Accumulate male counts
    if (length(idx_male) > 0) {
      s_male <- if (length(idx_male) == 1) {
        as.numeric(X[, idx_male])
      } else {
        as.numeric(Matrix::rowSums(X[, idx_male, drop = FALSE]))
      }
      
      nz_male <- if (length(idx_male) == 1) {
        as.numeric(X[, idx_male] > 0)
      } else {
        as.numeric(Matrix::rowSums(X[, idx_male, drop = FALSE] > 0))
      }
      
      male_sum <- if (is.null(male_sum)) s_male else male_sum + s_male
      nz_male_sum <- if (is.null(nz_male_sum)) nz_male else nz_male_sum + nz_male
      n_male_total <- n_male_total + length(idx_male)
    }
    
    rm(X)
    gc()
  }
  
  # Skip if missing one sex
  if (is.null(female_sum) || is.null(male_sum) || 
      n_female_total == 0 || n_male_total == 0) {
    message("    Skipping (missing one sex)")
    next
  }
  
  names(female_sum) <- names(male_sum) <- genes_counts_ref
  names(nz_female_sum) <- names(nz_male_sum) <- genes_counts_ref
  
  # Calculate percentage of cells expressing each gene
  pct_male <- nz_male_sum / n_male_total
  pct_female <- nz_female_sum / n_female_total
  
  # Filter genes by minimum percentage
  keep_pct <- (pct_male >= params$min_pct_within_cluster) | 
    (pct_female >= params$min_pct_within_cluster)
  
  # Filter by specificity
  spec_vec <- SPECIFICITY[match(genes_counts_ref, rownames(SPECIFICITY)), cluster]
  spec_vec[is.na(spec_vec)] <- 0
  keep_spec <- spec_vec > params$specificity_cutoff
  
  # Combined filter
  keep_idx <- which(keep_pct & keep_spec)
  
  if (!length(keep_idx)) {
    message("    No genes passed filters")
    next
  }
  
  message("    Genes passing filters: ", length(keep_idx))
  message("    Cells - Female: ", n_female_total, ", Male: ", n_male_total)
  
  # Create pseudobulk matrix for edgeR
  # Using single pseudobulk per sex (2-column design)
  M2 <- cbind(
    Female = female_sum[keep_idx],
    Male = male_sum[keep_idx]
  )
  rownames(M2) <- genes_counts_ref[keep_idx]
  storage.mode(M2) <- "numeric"
  
  # Run edgeR with fallback dispersion (2 samples = insufficient for estimation)
  results <- edger_run(M2, factor(c("Female", "Male"), levels = c("Female", "Male")),
                       use_disp_fallback = TRUE, 
                       disp_fallback = params$disp_fallback)
  
  if (is.null(results) || !nrow(results)) {
    message("    No genes passed edgeR filtering")
    next
  }
  
  # Add percentage columns
  add_df <- tibble(
    gene = genes_counts_ref[keep_idx],
    pct_Male = as.numeric(pct_male[keep_idx]),
    pct_Female = as.numeric(pct_female[keep_idx]),
    specificity = as.numeric(spec_vec[keep_idx])
  )
  
  results <- results %>%
    left_join(add_df, by = "gene") %>%
    mutate(p_val_adj = ifelse(is.finite(p_val_adj), p_val_adj, 1))
  
  message("    DE genes: ", nrow(results))
  
  # Add to Excel workbook
  sheet_name <- gsub("[\\/:*?\\[\\]]", "_", cluster)
  sheet_name <- substr(sheet_name, 1, 31)
  addWorksheet(wb_clusters, sheet_name)
  writeData(wb_clusters, sheet_name, results)
  
  # Create volcano plot
  total_tested <- nrow(results)
  n_up <- sum(results$p_val_adj < params$p_adj_threshold & 
                results$avg_log2FC > params$lfc_threshold, na.rm = TRUE)
  n_down <- sum(results$p_val_adj < params$p_adj_threshold & 
                  results$avg_log2FC < -params$lfc_threshold, na.rm = TRUE)
  
  caption_text <- sprintf(
    "Total: %s | Up (Male>Female): %s | Down (Female>Male): %s",
    format(total_tested, big.mark = ","),
    format(n_up, big.mark = ","),
    format(n_down, big.mark = ",")
  )
  
  safe_cluster <- gsub("[^A-Za-z0-9_-]+", "_", cluster)
  
  make_volcano(
    results,
    title = paste0("Cluster: ", cluster, " (Male vs Female)"),
    outbase = file.path(plots_dir, paste0("volcano_", safe_cluster)),
    p_adj_thr = params$p_adj_threshold,
    lfc_thr = params$lfc_threshold,
    caption = caption_text
  )
  
  gc()
}

# Save per-cluster results
if (length(wb_clusters$worksheets) > 0) {
  excel_file <- file.path(output_dir, "DE_sex_per_cluster.xlsx")
  saveWorkbook(wb_clusters, excel_file, overwrite = TRUE)
  message("\nSaved per-cluster DE results: ", excel_file)
} else {
  message("\nNo per-cluster results to save")
}

# =============================================================================
# WHOLE-ADRENAL DIFFERENTIAL EXPRESSION
# =============================================================================

message("\n", rep("=", 80))
message("DIFFERENTIAL EXPRESSION: WHOLE ADRENAL (Male vs Female)")
message(rep("=", 80))

# Aggregate counts across all cells
whole_female_sum <- NULL
whole_male_sum <- NULL
whole_nz_female <- NULL
whole_nz_male <- NULL
n_female <- 0L
n_male <- 0L

for (layer in counts_layers) {
  X <- LayerData(combined[["Spatial"]], layer = layer)
  
  # Align genes
  if (!identical(rownames(X), genes_counts_ref)) {
    common <- intersect(genes_counts_ref, rownames(X))
    X2 <- Matrix::Matrix(0, nrow = length(genes_counts_ref), ncol = ncol(X), sparse = TRUE)
    rownames(X2) <- genes_counts_ref
    colnames(X2) <- colnames(X)
    if (length(common)) {
      X2[match(common, genes_counts_ref), ] <- X[common, , drop = FALSE]
    }
    X <- X2
    rm(X2)
  }
  
  sx <- sex_vec[colnames(X)]
  keep <- !is.na(sx)
  
  if (!any(keep)) {
    rm(X)
    next
  }
  
  sx <- sx[keep]
  X <- X[, keep, drop = FALSE]
  
  # Female cells
  cells_female <- which(sx == "Female")
  if (length(cells_female) > 0) {
    s_f <- as.numeric(Matrix::rowSums(X[, cells_female, drop = FALSE]))
    nz_f <- as.numeric(Matrix::rowSums(X[, cells_female, drop = FALSE] > 0))
    
    whole_female_sum <- if (is.null(whole_female_sum)) s_f else whole_female_sum + s_f
    whole_nz_female <- if (is.null(whole_nz_female)) nz_f else whole_nz_female + nz_f
    n_female <- n_female + length(cells_female)
  }
  
  # Male cells
  cells_male <- which(sx == "Male")
  if (length(cells_male) > 0) {
    s_m <- as.numeric(Matrix::rowSums(X[, cells_male, drop = FALSE]))
    nz_m <- as.numeric(Matrix::rowSums(X[, cells_male, drop = FALSE] > 0))
    
    whole_male_sum <- if (is.null(whole_male_sum)) s_m else whole_male_sum + s_m
    whole_nz_male <- if (is.null(whole_nz_male)) nz_m else whole_nz_male + nz_m
    n_male <- n_male + length(cells_male)
  }
  
  rm(X)
  gc()
}

message("Total cells - Female: ", n_female, ", Male: ", n_male)

# Calculate percentages
pct_female <- whole_nz_female / pmax(1, n_female)
pct_male <- whole_nz_male / pmax(1, n_male)
names(pct_female) <- names(pct_male) <- genes_counts_ref

# Filter by minimum percentage (at least one sex must have >= threshold)
min_pct_whole <- 0.01  # 1% for whole adrenal
keep_pct <- (pct_female >= min_pct_whole) | (pct_male >= min_pct_whole)
message("Genes passing min.pct filter: ", sum(keep_pct))

# Create pseudobulk matrix
names(whole_female_sum) <- names(whole_male_sum) <- genes_counts_ref
M2_whole <- cbind(
  Female = whole_female_sum[keep_pct],
  Male = whole_male_sum[keep_pct]
)

# Run edgeR
results_whole <- edger_run(
  M2_whole, 
  factor(c("Female", "Male"), levels = c("Female", "Male")),
  use_disp_fallback = TRUE,
  disp_fallback = params$disp_fallback
)

if (!is.null(results_whole) && nrow(results_whole) > 0) {
  # Add percentage columns
  results_whole <- results_whole %>%
    left_join(
      tibble(
        gene = genes_counts_ref[keep_pct],
        pct_Male = as.numeric(pct_male[keep_pct]),
        pct_Female = as.numeric(pct_female[keep_pct])
      ),
      by = "gene"
    )
  
  message("DE genes (whole adrenal): ", nrow(results_whole))
  
  # Save to Excel
  wb_whole <- createWorkbook()
  addWorksheet(wb_whole, "DE_whole_adrenal")
  writeData(wb_whole, "DE_whole_adrenal", results_whole)
  
  # Add summary sheet
  summary_whole <- tibble(
    metric = c("cells_Female", "cells_Male", "genes_total", 
               "genes_after_filter", "DE_method"),
    value = c(n_female, n_male, length(genes_counts_ref), 
              sum(keep_pct), "edgeR_fallback_disp")
  )
  addWorksheet(wb_whole, "Summary")
  writeData(wb_whole, "Summary", summary_whole)
  
  excel_whole <- file.path(output_dir, "DE_sex_whole_adrenal.xlsx")
  saveWorkbook(wb_whole, excel_whole, overwrite = TRUE)
  message("Saved whole-adrenal DE results: ", excel_whole)
  
  # Create volcano plot
  total_tested <- nrow(results_whole)
  n_up <- sum(results_whole$p_val_adj < params$p_adj_threshold & 
                results_whole$avg_log2FC > params$lfc_threshold, na.rm = TRUE)
  n_down <- sum(results_whole$p_val_adj < params$p_adj_threshold & 
                  results_whole$avg_log2FC < -params$lfc_threshold, na.rm = TRUE)
  
  caption_text <- sprintf(
    "Total: %s | Up (Male>Female): %s | Down (Female>Male): %s",
    format(total_tested, big.mark = ","),
    format(n_up, big.mark = ","),
    format(n_down, big.mark = ",")
  )
  
  make_volcano(
    results_whole,
    title = "Whole Adrenal: Male vs Female",
    outbase = file.path(plots_dir, "volcano_whole_adrenal"),
    p_adj_thr = params$p_adj_threshold,
    lfc_thr = params$lfc_threshold,
    caption = caption_text
  )
}

# =============================================================================
# SUMMARY
# =============================================================================

message("\n", rep("=", 80))
message("SEX-DIFFERENTIAL EXPRESSION ANALYSIS COMPLETE")
message(rep("=", 80))

message("\nGenerated files:")
message("  - DE_sex_per_cluster.xlsx")
message("  - DE_sex_whole_adrenal.xlsx")
message("  - Volcano plots in: ", plots_dir)

message("\nStatistical approach:")
message("  - Pseudobulk aggregation per sex")
message("  - edgeR quasi-likelihood with fixed dispersion (", 
        params$disp_fallback, ")")
message("  - Gene filtering: min.pct >= ", params$min_pct_within_cluster,
        ", specificity > ", params$specificity_cutoff)
message("  - Significance thresholds: padj < ", params$p_adj_threshold,
        ", |log2FC| > ", params$lfc_threshold)