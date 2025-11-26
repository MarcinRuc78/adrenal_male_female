# =============================================================================
# Visium HD Spatial Transcriptomics Analysis Pipeline
# Script 06: Cell-Cell Communication Analysis using CellChat
# =============================================================================
#
# Description:
#   This script performs cell-cell communication analysis using CellChat v2,
#   comparing signaling networks between male and female mouse adrenal glands.
#   The analysis leverages spatial coordinates from Visium HD cell segmentation
#   to incorporate distance-dependent communication probability.
#
# Input:
#   - Annotated Seurat object from Script 03
#   - Optional: CSV files with cell barcodes for subsetting
#
# Output:
#   - CellChat objects for male and female samples
#   - Diagnostic tables for cell distribution
#   - Communication network visualizations
#
# Key Features:
#   - Sex-stratified analysis with balanced downsampling
#   - Biological exception handling (e.g., X-zone preserved)
#   - Spatial distance-dependent communication modeling
#   - Optimized parameters for Visium HD cell segmentation data
#
# Methods:
#   - Secreted signaling pathway analysis (CellChatDB.mouse)
#   - Truncated mean expression calculation
#   - Distance-weighted communication probability
#   - Contact-dependent and paracrine signaling modes
#
# Dependencies:
#   - Seurat (>= 5.0)
#   - CellChat (>= 2.1.0, recommended 2.2.0)
#   - Matrix
#   - dplyr
#   - ggplot2
#   - future (parallel processing)
#   - RANN (nearest neighbor calculations)
#
# References:
#   - Jin S, et al. (2021). Inference and analysis of cell-cell communication
#     using CellChat. Nature Communications, 12(1), 1088.
#   - Jin S, et al. (2024). CellChat for systematic analysis of cell-cell
#     communication from single-cell transcriptomics. Nature Protocols.
#
# =============================================================================

# -----------------------------------------------------------------------------
# Load Required Packages
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(future)
  library(RANN)
})

set.seed(42)

# Check CellChat version
message("CellChat version: ", packageVersion("CellChat"))
if (packageVersion("CellChat") < "2.1.0") {
  warning("CellChat version < 2.1.0 detected. Update recommended: ",
          "devtools::install_github('jinworks/CellChat')")
}

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Output directory
output_dir <- "results"
cellchat_dir <- file.path(output_dir, "cellchat")
dir.create(cellchat_dir, showWarnings = FALSE, recursive = TRUE)

# Memory settings for large datasets
options(future.globals.maxSize = 48 * 1024^3)  # 48 GB

# -----------------------------------------------------------------------------
# Downsampling Configuration
# -----------------------------------------------------------------------------
# 
# Balanced downsampling ensures fair comparison between sexes by equalizing
# cell numbers per cluster. However, some clusters may show genuine biological
# sex differences (e.g., X-zone in adrenal glands) and should be excluded
# from downsampling to preserve these differences.

DOWNSAMPLE_CONFIG <- list(
  # Strategy options:
  # - "balanced_with_exceptions": Match to smallest/geometric mean, except excluded
  # - "match_female": Downsample male to match female counts
  # - "match_male": Downsample female to match male counts
  # - "match_smallest": Both sexes match smaller count
  # - "none": No downsampling
  strategy = "balanced_with_exceptions",
  
  # Clusters excluded from downsampling (biological sex differences)
  # X-zone (X/JMZ) shows sex-specific biology - larger in females
  # Add cluster names here that represent genuine biological differences
  exclude_clusters = c(),  # e.g., c("X/JMZ") if sex difference is biological
  
  # For non-excluded clusters: how to calculate target count
  # - "smallest": Use minimum of male/female counts
  # - "balanced": Use geometric mean sqrt(n_male * n_female)
  match_to = "smallest",
  
  # Minimum cells per cluster (prevents over-downsampling)
  min_cells_per_cluster = 50,
  
  # Random seed for reproducibility
  seed = 42
)

# -----------------------------------------------------------------------------
# Spatial Parameters for Visium HD
# -----------------------------------------------------------------------------
#
# These parameters define the spatial scale for cell-cell communication:
# - contact_range: Maximum distance for direct cell-cell contact signaling
# - interaction_range: Maximum distance for secreted/paracrine signaling
#
# For Visium HD with cell segmentation, distances are in pixels and
# should be calibrated based on actual cell sizes and tissue architecture.

SPATIAL_CONFIG <- list(
  # Multiplier for contact range (relative to median nearest neighbor distance)
  # Typical value: 2-3x cell diameter
  contact_multiplier = 2.5,
  
  
  # Multiplier for interaction range (paracrine signaling)
  # Typical value: 10-20x cell diameter (accounts for diffusion)
  interaction_multiplier = 15.0,
  
  # Nearest neighbor k for calibration
  k_main = 5,    # For base nearest neighbor distance
  k_scale = 50   # For 99th percentile calculation
)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# -----------------------------------------------------------------------------
# Layer Management for Seurat v5
# -----------------------------------------------------------------------------

#' Join Layers if Needed for Seurat v5 Compatibility
#'
#' Seurat v5 may have multiple layers (counts.1, counts.2, etc.) after merging.
#' This function joins them into single counts/data layers and ensures
#' normalized data is available.
#'
#' @param obj Seurat object
#' @param assay Assay name (default: "Spatial")
#' @return Seurat object with joined layers
join_layers_if_needed <- function(obj, assay = "Spatial") {
  stopifnot(assay %in% Assays(obj))
  
  layers <- Layers(obj[[assay]])
  counts_layers <- grep("^counts(\\.|$)", layers, value = TRUE)
  data_layers <- grep("^data(\\.|$)", layers, value = TRUE)
  
  # Join counts layers if multiple exist
  if (length(counts_layers) > 1 && !"counts" %in% layers) {
    obj <- JoinLayers(obj, assay = assay, layers = counts_layers, 
                      new.layer = "counts")
  }
  
  # Join data layers if multiple exist
  if (length(data_layers) > 1 && !"data" %in% Layers(obj[[assay]])) {
    obj <- JoinLayers(obj, assay = assay, layers = data_layers, 
                      new.layer = "data")
  }
  
  # Normalize if data layer doesn't exist
  if (!"data" %in% Layers(obj[[assay]])) {
    DefaultLayer(obj[[assay]]) <- "counts"
    obj <- NormalizeData(obj, assay = assay, 
                         normalization.method = "LogNormalize",
                         scale.factor = 1e4, verbose = FALSE)
  }
  
  DefaultLayer(obj[[assay]]) <- "data"
  return(obj)
}

# -----------------------------------------------------------------------------
# Spatial Coordinate Extraction
# -----------------------------------------------------------------------------

#' Extract Spatial Coordinates from Seurat Object
#'
#' Attempts to extract spatial coordinates from various possible locations
#' in the Seurat object, handling different naming conventions.
#'
#' @param obj Seurat object
#' @return Data frame with x, y columns and cell names as row names
#'
#' @details
#' Searches for coordinates in the following order:
#' 1. GetTissueCoordinates() function
#' 2. Metadata columns: centroid_x/y, x_centroid/y_centroid, etc.
#' 3. Standard columns: imagecol/imagerow, pxl_col/pxl_row, x/y
get_spatial_coords <- function(obj) {
  
  # Try GetTissueCoordinates first
  tc <- tryCatch(GetTissueCoordinates(obj), error = function(e) NULL)
  
  if (!is.null(tc) && all(c("x", "y") %in% colnames(tc))) {
    row_names <- if ("cell" %in% colnames(tc)) tc$cell else Cells(obj)
    coords <- data.frame(
      x = tc[, "x"],
      y = tc[, "y"],
      row.names = row_names,
      check.names = FALSE
    )
    message("Coordinates from GetTissueCoordinates: ", nrow(coords), " cells")
    return(coords[Cells(obj), , drop = FALSE])
  }
  
  # Search metadata for coordinate columns
  meta <- obj@meta.data
  col_names_lower <- tolower(colnames(meta))
  
  # Centroid-based coordinates (from cell segmentation)
  centroid_patterns <- list(
    c("centroid_x", "centroid_y"),
    c("x_centroid", "y_centroid"),
    c("center_x", "center_y"),
    c("x_location", "y_location")
  )
  
  for (pattern in centroid_patterns) {
    if (all(pattern %in% col_names_lower)) {
      x_col <- colnames(meta)[match(pattern[1], col_names_lower)]
      y_col <- colnames(meta)[match(pattern[2], col_names_lower)]
      message("Coordinates from centroid columns: ", x_col, ", ", y_col)
      return(data.frame(
        x = meta[[x_col]],
        y = meta[[y_col]],
        row.names = rownames(meta),
        check.names = FALSE
      )[Cells(obj), , drop = FALSE])
    }
  }
  
  # Standard coordinate columns
  standard_patterns <- list(
    c("imagecol", "imagerow"),
    c("pxl_col_in_fullres", "pxl_row_in_fullres"),
    c("x", "y"),
    c("col", "row")
  )
  
  for (pattern in standard_patterns) {
    if (all(pattern %in% col_names_lower)) {
      x_col <- colnames(meta)[match(pattern[1], col_names_lower)]
      y_col <- colnames(meta)[match(pattern[2], col_names_lower)]
      message("Coordinates from standard columns: ", x_col, ", ", y_col)
      return(data.frame(
        x = meta[[x_col]],
        y = meta[[y_col]],
        row.names = rownames(meta),
        check.names = FALSE
      )[Cells(obj), , drop = FALSE])
    }
  }
  
  stop("No spatial coordinates found in Seurat object")
}

#' Align Metadata and Coordinates with Seurat Object
#'
#' Ensures metadata, coordinates, and cell names are aligned,
#' removing cells with invalid coordinates.
#'
#' @param obj Seurat object
#' @return List with aligned obj, meta, and coords
align_meta_and_coords <- function(obj) {
  cell_names <- Cells(obj)
  
  # Create metadata frame
  meta <- data.frame(
    labels = as.character(Idents(obj)),
    samples = if ("sample" %in% colnames(obj@meta.data)) {
      as.character(obj$sample)
    } else {
      "sample1"
    },
    row.names = cell_names,
    check.names = FALSE
  )
  
  # Get coordinates
  coords <- get_spatial_coords(obj)
  coords$x <- suppressWarnings(as.numeric(coords$x))
  coords$y <- suppressWarnings(as.numeric(coords$y))
  
  # Remove cells with invalid coordinates
  valid_coords <- rowSums(is.finite(as.matrix(coords))) == 2
  n_invalid <- sum(!valid_coords)
  
  if (n_invalid > 0) {
    message("Removing ", n_invalid, " cells with invalid coordinates")
    keep_cells <- rownames(coords)[valid_coords]
    obj <- subset(obj, cells = keep_cells)
    meta <- meta[keep_cells, , drop = FALSE]
    coords <- coords[keep_cells, , drop = FALSE]
  }
  
  # Verify alignment
  stopifnot(all(Cells(obj) == rownames(meta)))
  stopifnot(all(Cells(obj) == rownames(coords)))
  
  message("Final: ", nrow(coords), " cells with valid coordinates")
  
  return(list(obj = obj, meta = meta, coords = coords))
}

# -----------------------------------------------------------------------------
# Distance Calibration for Spatial Analysis
# -----------------------------------------------------------------------------

#' Calibrate Spatial Distance Parameters
#'
#' Calculates optimal distance parameters for CellChat based on the
#' actual spatial distribution of cells in the tissue.
#'
#' @param coords Data frame with x, y coordinates
#' @param k_main Number of neighbors for median NN distance
#' @param k_scale Number of neighbors for 99th percentile
#' @return List with base_nn (median NN distance) and scale_distance
#'
#' @details
#' For Visium HD with cell segmentation:
#' - base_nn represents typical cell-cell distance (approximately cell diameter)
#' - scale_distance is the decay factor for distance-weighted probability
#'
#' The function automatically adjusts scale_distance to reasonable bounds
#' for paracrine signaling detection.
calibrate_spatial_distance <- function(coords, 
                                       k_main = SPATIAL_CONFIG$k_main, 
                                       k_scale = SPATIAL_CONFIG$k_scale) {
  n <- nrow(coords)
  stopifnot(n >= 3)
  
  coord_matrix <- as.matrix(coords)
  
  # Calculate median nearest neighbor distance
  k_use <- min(max(2, k_main), n)
  nn_result <- RANN::nn2(coord_matrix, k = k_use)
  
  stopifnot(ncol(nn_result$nn.dists) >= 2)
  
  # Get valid distances (exclude self, which is always 0)
  valid_distances <- nn_result$nn.dists[, 2]
  valid_distances <- valid_distances[is.finite(valid_distances) & valid_distances > 0]
  base_nn <- median(valid_distances)
  
  # Calculate 99th percentile for scale normalization
  k_scale_use <- min(max(5, k_scale), n)
  nn_scale <- RANN::nn2(coord_matrix, k = k_scale_use)
  distance_matrix <- nn_scale$nn.dists[, -1, drop = FALSE]  # Exclude self
  q99 <- stats::quantile(as.numeric(distance_matrix), 0.99, 
                         na.rm = TRUE, names = FALSE)
  
  # Calculate and bound scale factor
  raw_scale <- as.numeric(1 / q99)
  
  # Adjust scale for Visium HD cell segmentation
  # Too high: overly restrictive, misses paracrine signaling
  # Too low: includes spurious long-range interactions
  if (raw_scale > 10) {
    message("Auto-scale too high (", round(raw_scale, 2), "), capping at 5.0")
    scale_distance <- 5.0
  } else if (raw_scale < 0.5) {
    message("Auto-scale too low (", round(raw_scale, 2), "), setting to 2.0")
    scale_distance <- 2.0
  } else {
    scale_distance <- raw_scale
  }
  
  message("\nSpatial calibration results:")
  message("  Median NN distance: ", round(base_nn, 2), " px")
  message("  99th percentile distance: ", round(q99, 2), " px")
  message("  Scale factor: ", round(scale_distance, 3))
  
  return(list(base_nn = base_nn, scale_distance = scale_distance))
}

# -----------------------------------------------------------------------------
# Downsampling Functions
# -----------------------------------------------------------------------------

#' Perform Balanced Downsampling with Biological Exceptions
#'
#' Downsamples cells to balance sex representation while preserving
#' clusters with genuine biological sex differences.
#'
#' @param obj_male Seurat object for male sample
#' @param obj_female Seurat object for female sample
#' @param config Downsampling configuration list
#' @return List with downsampled male and female Seurat objects
#'
#' @details
#' The function:
#' 1. Identifies common clusters between sexes
#' 2. Excludes specified clusters from downsampling (biological exceptions)
#' 3. Calculates target cell counts based on strategy
#' 4. Randomly samples cells to reach targets
#' 5. Reports balance metrics before and after
downsample_by_strategy <- function(obj_male, obj_female, 
                                   config = DOWNSAMPLE_CONFIG) {
  set.seed(config$seed)
  
  # Get cluster counts
  clusters_m <- table(Idents(obj_male))
  clusters_f <- table(Idents(obj_female))
  
  # Find common clusters
  common_clusters <- intersect(names(clusters_m), names(clusters_f))
  
  message("\n", rep("=", 60))
  message("DOWNSAMPLING STRATEGY: ", config$strategy)
  message(rep("=", 60))
  
  # Handle "no downsampling" case
  if (config$strategy == "none") {
    message("No downsampling applied")
    return(list(male = obj_male, female = obj_female))
  }
  
  # Identify excluded clusters (biological differences)
  excluded <- if (!is.null(config$exclude_clusters)) {
    intersect(config$exclude_clusters, common_clusters)
  } else {
    character(0)
  }
  
  if (length(excluded) > 0) {
    message("\nEXCLUDED from downsampling (biological sex differences):")
    message("  ", paste(excluded, collapse = ", "))
    message("  These clusters retain original cell counts\n")
  }
  
  # Calculate targets per cluster
  targets <- list()
  
  for (cluster in common_clusters) {
    n_m <- as.numeric(clusters_m[cluster])
    n_f <- as.numeric(clusters_f[cluster])
    
    # Check if excluded
    if (cluster %in% excluded) {
      targets[[cluster]] <- list(male = n_m, female = n_f, excluded = TRUE)
      message(sprintf("  %s: EXCLUDED - Male %d | Female %d (preserved)",
                      cluster, n_m, n_f))
      next
    }
    
    # Calculate target for non-excluded clusters
    if (config$strategy == "balanced_with_exceptions") {
      target <- switch(
        config$match_to,
        "smallest" = min(n_m, n_f),
        "balanced" = round(sqrt(n_m * n_f)),
        min(n_m, n_f)  # default
      )
    } else {
      target <- switch(
        config$strategy,
        "match_female" = n_f,
        "match_male" = n_m,
        "match_smallest" = min(n_m, n_f),
        "balanced" = round(sqrt(n_m * n_f))
      )
    }
    
    # Enforce minimum
    target <- max(target, config$min_cells_per_cluster)
    
    # Don't upsample - cap at actual counts
    target_m <- min(target, n_m)
    target_f <- min(target, n_f)
    
    targets[[cluster]] <- list(male = target_m, female = target_f, excluded = FALSE)
    
    ratio_before <- n_m / n_f
    ratio_after <- target_m / target_f
    
    message(sprintf("  %s: Male %d->%d | Female %d->%d | Ratio %.2fx->%.2fx",
                    cluster, n_m, target_m, n_f, target_f, 
                    ratio_before, ratio_after))
  }
  
  # Perform downsampling - Male
  cells_keep_m <- c()
  for (cluster in common_clusters) {
    cells_cluster <- WhichCells(obj_male, idents = cluster)
    n_target <- targets[[cluster]]$male
    
    if (length(cells_cluster) > n_target) {
      cells_cluster <- sample(cells_cluster, n_target)
    }
    cells_keep_m <- c(cells_keep_m, cells_cluster)
  }
  
  obj_male_ds <- subset(obj_male, cells = cells_keep_m)
  Idents(obj_male_ds) <- droplevels(Idents(obj_male_ds))
  
  # Perform downsampling - Female
  cells_keep_f <- c()
  for (cluster in common_clusters) {
    cells_cluster <- WhichCells(obj_female, idents = cluster)
    n_target <- targets[[cluster]]$female
    
    if (length(cells_cluster) > n_target) {
      cells_cluster <- sample(cells_cluster, n_target)
    }
    cells_keep_f <- c(cells_keep_f, cells_cluster)
  }
  
  obj_female_ds <- subset(obj_female, cells = cells_keep_f)
  Idents(obj_female_ds) <- droplevels(Idents(obj_female_ds))
  
  # Summary
  pct_removed_m <- round((ncol(obj_male) - ncol(obj_male_ds)) / 
                           ncol(obj_male) * 100, 1)
  pct_removed_f <- round((ncol(obj_female) - ncol(obj_female_ds)) / 
                           ncol(obj_female) * 100, 1)
  
  message("\nDownsampling complete:")
  message("  Male:   ", ncol(obj_male), " -> ", ncol(obj_male_ds),
          " cells (", pct_removed_m, "% removed)")
  message("  Female: ", ncol(obj_female), " -> ", ncol(obj_female_ds),
          " cells (", pct_removed_f, "% removed)")
  
  # Final ratio check
  final_m <- table(Idents(obj_male_ds))
  final_f <- table(Idents(obj_female_ds))
  
  message("\nFinal Male/Female ratios:")
  for (cluster in common_clusters) {
    ratio <- as.numeric(final_m[cluster]) / as.numeric(final_f[cluster])
    
    if (cluster %in% excluded) {
      status <- "[PRESERVED]"
    } else if (ratio < 1.2) {
      status <- "BALANCED"
    } else if (ratio < 1.5) {
      status <- "OK"
    } else {
      status <- "IMBALANCED"
    }
    
    message(sprintf("  %s: %.2fx %s", cluster, ratio, status))
  }
  
  return(list(male = obj_male_ds, female = obj_female_ds))
}

#' Compare Cell Distribution Between Sexes
#'
#' Creates a detailed comparison table of cell counts per cluster.
#'
#' @param obj_m Male Seurat object
#' @param obj_f Female Seurat object
#' @param label Optional label for output
#' @param excluded_clusters Clusters to mark as excluded
#' @return Data frame with comparison statistics
compare_cell_distribution <- function(obj_m, obj_f, label = "", 
                                      excluded_clusters = NULL) {
  tab_m <- table(Idents(obj_m))
  tab_f <- table(Idents(obj_f))
  all_clusters <- union(names(tab_m), names(tab_f))
  
  df <- data.frame(
    Cluster = all_clusters,
    Male_n = as.numeric(tab_m[all_clusters]),
    Male_pct = round(as.numeric(tab_m[all_clusters]) / sum(tab_m) * 100, 1),
    Female_n = as.numeric(tab_f[all_clusters]),
    Female_pct = round(as.numeric(tab_f[all_clusters]) / sum(tab_f) * 100, 1)
  )
  
  # Handle NAs
  df$Male_n[is.na(df$Male_n)] <- 0
  df$Female_n[is.na(df$Female_n)] <- 0
  df$Male_pct[is.na(df$Male_pct)] <- 0
  df$Female_pct[is.na(df$Female_pct)] <- 0
  
  df$Ratio <- round(df$Male_n / df$Female_n, 2)
  df$Diff_pct <- round(df$Male_pct - df$Female_pct, 1)
  
  # Add status
  df$Status <- sapply(seq_len(nrow(df)), function(i) {
    cluster <- df$Cluster[i]
    ratio <- df$Ratio[i]
    
    if (!is.null(excluded_clusters) && cluster %in% excluded_clusters) {
      "Biological"
    } else if (ratio < 1.2) {
      "Balanced"
    } else if (ratio < 1.5) {
      "Acceptable"
    } else {
      "Imbalanced"
    }
  })
  
  df <- df[order(-df$Male_n), ]
  
  if (nzchar(label)) {
    message("\n", label)
    message(rep("=", nchar(label)))
  }
  
  print(df, row.names = FALSE)
  
  invisible(df)
}

# =============================================================================
# MAIN CELLCHAT BUILDING FUNCTION
# =============================================================================

#' Build CellChat Object for Visium HD Spatial Data
#'
#' Creates and processes a CellChat object from a Seurat object,
#' optimized for Visium HD with cell segmentation.
#'
#' @param seurat_obj Seurat object with spatial data
#' @param out_dir Output directory for saving results
#' @param label_suffix Label for this analysis (e.g., "male", "female")
#' @return CellChat object with computed communication probabilities
#'
#' @details
#' The function performs the following steps:
#' 1. Prepares Seurat object (joins layers, extracts coordinates)
#' 2. Calibrates spatial distance parameters
#' 3. Creates CellChat object with mouse secreted signaling database
#' 4. Identifies overexpressed genes and interactions
#' 5. Computes communication probabilities with spatial constraints
#' 6. Aggregates networks at pathway level
#'
#' For Visium HD cell segmentation data, the function uses cell centroids
#' rather than spot coordinates, enabling more accurate distance calculations.
build_cellchat_mouse <- function(seurat_obj, out_dir, label_suffix) {
  
  message("\n", rep("=", 60))
  message("Building CellChat: ", label_suffix)
  message(rep("=", 60))
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # -------------------------------------------------------------------------
  # Prepare Seurat object
  # -------------------------------------------------------------------------
  
  DefaultAssay(seurat_obj) <- "Spatial"
  seurat_obj <- join_layers_if_needed(seurat_obj, "Spatial")
  
  # Align metadata and coordinates
  aligned <- align_meta_and_coords(seurat_obj)
  seurat_prepared <- aligned$obj
  meta <- aligned$meta
  coords <- aligned$coords
  
  # -------------------------------------------------------------------------
  # Calibrate spatial parameters
  # -------------------------------------------------------------------------
  
  calibration <- calibrate_spatial_distance(coords)
  base_nn <- calibration$base_nn
  scale_distance <- calibration$scale_distance
  
  # Calculate range parameters based on median NN distance
  contact_range <- SPATIAL_CONFIG$contact_multiplier * base_nn
  interaction_range <- SPATIAL_CONFIG$interaction_multiplier * base_nn
  
  message("\nSpatial parameters for Visium HD (adrenal gland):")
  message("  Contact range: ", round(contact_range, 2), 
          " px (~", round(SPATIAL_CONFIG$contact_multiplier, 1), 
          "x cell diameter)")
  message("  Interaction range: ", round(interaction_range, 2), 
          " px (~", round(SPATIAL_CONFIG$interaction_multiplier, 1), 
          "x cell diameter)")
  message("  Scale distance: ", round(scale_distance, 3))
  
  # -------------------------------------------------------------------------
  # Create CellChat object
  # -------------------------------------------------------------------------
  
  message("\nCreating CellChat object...")
  message("  Mode: Cell-based (centroids from Visium HD segmentation)")
  
  cellchat <- createCellChat(
    object = seurat_prepared,
    meta = meta,
    group.by = "labels"
  )
  
  # Add spatial coordinates
  cellchat@images <- list(coordinates = coords)
  
  # Set spatial factors (CellChat >= 2.1.0)
  available_slots <- slotNames(cellchat)
  
  if ("spatial.factors" %in% available_slots) {
    tryCatch({
      cellchat@spatial.factors <- data.frame(ratio = 1, tol = contact_range)
      message("  Set @spatial.factors (CellChat >= 2.1)")
    }, error = function(e) {
      cellchat@images$spatial.factors <- data.frame(ratio = 1, tol = contact_range)
      message("  Set @images$spatial.factors (fallback)")
    })
  } else {
    cellchat@images$spatial.factors <- data.frame(ratio = 1, tol = contact_range)
    message("  Set @images$spatial.factors (CellChat < 2.1)")
  }
  
  message("  Cells: ", nrow(coords))
  message("  Clusters: ", length(unique(meta$labels)))
  
  # -------------------------------------------------------------------------
  # Load ligand-receptor database
  # -------------------------------------------------------------------------
  
  # Use mouse database with secreted signaling pathways
  # Secreted signaling is most relevant for spatial analysis as it
  # involves distance-dependent ligand diffusion
  CellChatDB <- CellChatDB.mouse
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", 
                             key = "annotation")
  cellchat@DB <- CellChatDB.use
  
  message("  Database: ", nrow(CellChatDB.use$interaction), 
          " secreted signaling interactions")
  
  # -------------------------------------------------------------------------
  # Prepare expression data
  # -------------------------------------------------------------------------
  
  DefaultAssay(seurat_prepared) <- "Spatial"
  
  # Ensure data layer exists
  if (!"data" %in% Layers(seurat_prepared[["Spatial"]])) {
    DefaultLayer(seurat_prepared[["Spatial"]]) <- "counts"
    seurat_prepared <- NormalizeData(seurat_prepared, assay = "Spatial",
                                     normalization.method = "LogNormalize",
                                     scale.factor = 1e4, verbose = FALSE)
  }
  DefaultLayer(seurat_prepared[["Spatial"]]) <- "data"
  
  # Remove duplicated genes if any
  gene_duplicated <- duplicated(rownames(seurat_prepared))
  if (any(gene_duplicated)) {
    message("  Removing ", sum(gene_duplicated), " duplicated genes")
    seurat_prepared <- seurat_prepared[!gene_duplicated, ]
  }
  
  # Load expression data into CellChat
  cellchat@data <- as(
    GetAssayData(seurat_prepared, assay = "Spatial", 
                 layer = DefaultLayer(seurat_prepared[["Spatial"]])),
    "dgCMatrix"
  )
  
  # Subset to signaling genes
  cellchat <- subsetData(cellchat)
  
  if (is.null(cellchat@data.signaling) || nrow(cellchat@data.signaling) == 0) {
    stop("subsetData returned empty matrix - no signaling genes found")
  }
  
  message("  Signaling genes in data: ", nrow(cellchat@data.signaling))
  
  # -------------------------------------------------------------------------
  # Identify overexpressed genes and interactions
  # -------------------------------------------------------------------------
  
  # Set up parallel processing (limit workers for stability)
  n_cores <- max(1, min(6, parallel::detectCores() - 2))
  future::plan("multisession", workers = n_cores)
  message("  Parallel processing: ", n_cores, " workers")
  
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  message("  Overexpressed genes and interactions identified")
  
  # -------------------------------------------------------------------------
  # Compute communication probabilities
  # -------------------------------------------------------------------------
  
  message("\nComputing communication probabilities...")
  message("  Method: truncatedMean (trim = 0.1)")
  message("  Distance-based: YES")
  message("  Contact-dependent: YES")
  
  # Try full spatial model first, with fallbacks
  cellchat <- tryCatch({
    computeCommunProb(
      cellchat,
      type = "truncatedMean",
      trim = 0.1,
      distance.use = TRUE,
      interaction.range = interaction_range,
      scale.distance = scale_distance,
      contact.dependent = TRUE,
      contact.range = contact_range
    )
  }, error = function(e) {
    message("  Full spatial model failed, trying simplified...")
    
    tryCatch({
      computeCommunProb(
        cellchat,
        type = "truncatedMean",
        trim = 0.1,
        distance.use = TRUE,
        scale.distance = scale_distance
      )
    }, error = function(e2) {
      message("  Distance model failed, using basic model...")
      computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1)
    })
  })
  
  message("  Communication probabilities computed")
  
  # -------------------------------------------------------------------------
  # Filter, compute pathways, aggregate
  # -------------------------------------------------------------------------
  
  # Filter interactions with too few cells
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  message("  Filtered communications (min.cells = 10)")
  
  # Compute pathway-level probabilities
  cellchat <- computeCommunProbPathway(cellchat)
  message("  Pathway-level probabilities computed")
  
  # Aggregate network
  cellchat <- aggregateNet(cellchat)
  message("  Network aggregated")
  
  # -------------------------------------------------------------------------
  # Summary statistics
  # -------------------------------------------------------------------------
  
  n_interactions <- sum(cellchat@net$count > 0)
  total_weight <- sum(cellchat@net$weight)
  n_pathways <- length(cellchat@netP)
  
  message("\nCellChat Results Summary:")
  message("  Total interactions detected: ", n_interactions)
  message("  Total communication weight: ", round(total_weight, 6))
  message("  Number of signaling pathways: ", n_pathways)
  
  # Check for specific pathways of interest
  if (n_pathways > 0) {
    pathway_names <- names(cellchat@netP)
    message("\nKey pathways detected:")
    
    # Pathways relevant for adrenal gland
    key_pathways <- c("WNT", "HH", "BMP", "NOTCH", "IGF", "VEGF", "PDGF")
    
    for (pw in key_pathways) {
      if (pw %in% pathway_names) {
        pw_weight <- sum(cellchat@netP[[pw]]$prob, na.rm = TRUE)
        message("  ", pw, ": DETECTED (weight = ", round(pw_weight, 6), ")")
      }
    }
  }
  
  # -------------------------------------------------------------------------
  # Save result
  # -------------------------------------------------------------------------
  
  rds_file <- file.path(out_dir, sprintf("cellchat_%s.rds", label_suffix))
  saveRDS(cellchat, file = rds_file)
  message("\nSaved: ", rds_file)
  
  # Reset parallel processing
  future::plan("sequential")
  
  invisible(cellchat)
}

# =============================================================================
# MAIN WORKFLOW
# =============================================================================

message("\n", rep("=", 80))
message("CELLCHAT ANALYSIS PIPELINE")
message(rep("=", 80))

# -----------------------------------------------------------------------------
# Load Data
# -----------------------------------------------------------------------------

message("\nLoading annotated Seurat object...")
combined <- readRDS(file.path(output_dir, "combined_annotated.rds"))

stopifnot(inherits(combined, "Seurat"))
stopifnot("Spatial" %in% Assays(combined))

message("Loaded object with ", ncol(combined), " cells")

# -----------------------------------------------------------------------------
# Split by Sex
# -----------------------------------------------------------------------------

message("\n", rep("=", 60))
message("SPLITTING BY SEX")
message(rep("=", 60))

# Find sex column
sex_col <- if ("sex" %in% colnames(combined@meta.data)) {
  "sex"
} else if ("sex2" %in% colnames(combined@meta.data)) {
  "sex2"
} else {
  stop("Sex column not found in metadata")
}

# Normalize sex labels
normalize_sex_label <- function(x) {
  x_lower <- tolower(trimws(as.character(x)))
  ifelse(x_lower %in% c("female", "f"), "Female",
         ifelse(x_lower %in% c("male", "m"), "Male", NA))
}

sex_normalized <- normalize_sex_label(combined@meta.data[[sex_col]])
names(sex_normalized) <- colnames(combined)

# Get cells by sex
cells_female <- names(sex_normalized)[sex_normalized == "Female" & !is.na(sex_normalized)]
cells_male <- names(sex_normalized)[sex_normalized == "Male" & !is.na(sex_normalized)]

message("Total cells: ", length(sex_normalized))
message("  Female: ", length(cells_female))
message("  Male: ", length(cells_male))
message("  NA: ", sum(is.na(sex_normalized)))

if (!length(cells_female) || !length(cells_male)) {
  stop("Missing cells for one sex - cannot proceed")
}

# Create sex-specific objects
female_obj <- subset(combined, cells = cells_female)
male_obj <- subset(combined, cells = cells_male)

Idents(female_obj) <- droplevels(Idents(female_obj))
Idents(male_obj) <- droplevels(Idents(male_obj))

# -----------------------------------------------------------------------------
# Remove Unwanted Clusters (Optional)
# -----------------------------------------------------------------------------

# Clusters to exclude from analysis (e.g., low quality, non-adrenal)
clusters_to_remove <- c("fibroblasts")  # Adjust as needed

if (length(clusters_to_remove) > 0) {
  present_to_remove <- intersect(
    clusters_to_remove,
    c(levels(Idents(female_obj)), levels(Idents(male_obj)))
  )
  
  if (length(present_to_remove) > 0) {
    message("\nRemoving clusters: ", paste(present_to_remove, collapse = ", "))
    
    female_obj <- subset(female_obj, idents = present_to_remove, invert = TRUE)
    male_obj <- subset(male_obj, idents = present_to_remove, invert = TRUE)
    
    Idents(female_obj) <- droplevels(Idents(female_obj))
    Idents(male_obj) <- droplevels(Idents(male_obj))
  }
}

message("\nAfter cluster filtering:")
message("  Male:   ", ncol(male_obj), " cells, ", 
        length(levels(Idents(male_obj))), " clusters")
message("  Female: ", ncol(female_obj), " cells, ", 
        length(levels(Idents(female_obj))), " clusters")

# -----------------------------------------------------------------------------
# Cell Distribution Before Downsampling
# -----------------------------------------------------------------------------

message("\n", rep("=", 60))
message("CELL DISTRIBUTION BEFORE DOWNSAMPLING")
message(rep("=", 60))

df_before <- compare_cell_distribution(
  male_obj, female_obj,
  label = "Cell counts per cluster",
  excluded_clusters = DOWNSAMPLE_CONFIG$exclude_clusters
)

# -----------------------------------------------------------------------------
# Perform Downsampling
# -----------------------------------------------------------------------------

downsampled <- downsample_by_strategy(male_obj, female_obj, DOWNSAMPLE_CONFIG)
male_obj_ds <- downsampled$male
female_obj_ds <- downsampled$female

# -----------------------------------------------------------------------------
# Cell Distribution After Downsampling
# -----------------------------------------------------------------------------

message("\n", rep("=", 60))
message("CELL DISTRIBUTION AFTER DOWNSAMPLING")
message(rep("=", 60))

df_after <- compare_cell_distribution(
  male_obj_ds, female_obj_ds,
  label = "Cell counts per cluster (downsampled)",
  excluded_clusters = DOWNSAMPLE_CONFIG$exclude_clusters
)

# Save diagnostic tables
diagnostics_dir <- file.path(cellchat_dir, "diagnostics")
dir.create(diagnostics_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(df_before,
          file.path(diagnostics_dir, "cluster_balance_BEFORE_downsampling.csv"),
          row.names = FALSE)
write.csv(df_after,
          file.path(diagnostics_dir, "cluster_balance_AFTER_downsampling.csv"),
          row.names = FALSE)

message("\nSaved diagnostic tables to: ", diagnostics_dir)

# -----------------------------------------------------------------------------
# Save Downsampled Seurat Objects
# -----------------------------------------------------------------------------

female_dir <- file.path(cellchat_dir, "female")
male_dir <- file.path(cellchat_dir, "male")
dir.create(female_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(male_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(female_obj_ds, file.path(female_dir, "female_seurat_downsampled.rds"))
saveRDS(male_obj_ds, file.path(male_dir, "male_seurat_downsampled.rds"))

message("Saved downsampled Seurat objects")

# -----------------------------------------------------------------------------
# Build CellChat Objects
# -----------------------------------------------------------------------------

message("\n", rep("=", 60))
message("BUILDING CELLCHAT OBJECTS")
message(rep("=", 60))

future::plan("sequential")

# Female CellChat
cellchat_female <- build_cellchat_mouse(
  female_obj_ds,
  female_dir,
  "female_downsampled"
)

# Male CellChat
cellchat_male <- build_cellchat_mouse(
  male_obj_ds,
  male_dir,
  "male_downsampled"
)

# -----------------------------------------------------------------------------
# Save Final CellChat Objects
# -----------------------------------------------------------------------------

saveRDS(cellchat_female, file.path(cellchat_dir, "cellchat_female_final.rds"))
saveRDS(cellchat_male, file.path(cellchat_dir, "cellchat_male_final.rds"))

# =============================================================================
# ANALYSIS COMPLETE
# =============================================================================

message("\n", rep("=", 80))
message("CELLCHAT ANALYSIS COMPLETED")
message(rep("=", 80))

message("\nAnalysis summary:")
message("  Type: Sex-stratified with balanced downsampling")
message("  Technology: Visium HD with cell segmentation")
message("  CellChat version: ", packageVersion("CellChat"))

message("\nOutput files:")
message("  - ", file.path(cellchat_dir, "cellchat_female_final.rds"))
message("  - ", file.path(cellchat_dir, "cellchat_male_final.rds"))
message("  - Diagnostic tables in: ", diagnostics_dir)

message("\nNext steps:")
message("  1. Run comparative analysis (Script 07)")
message("  2. Visualize communication networks")
message("  3. Identify sex-specific signaling pathways")