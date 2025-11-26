# =============================================================================
# Visium HD Spatial Transcriptomics Analysis Pipeline
# Script 01: Data Loading and Preprocessing
# =============================================================================
#
# Description:
#   This script loads 10x Genomics Visium HD spatial transcriptomics data from
#   multiple samples (male and female mouse adrenal glands), creates Seurat
#   objects, attaches spatial images, and performs initial quality control.
#
# Input:
#   - Visium HD output directories containing:
#     - filtered_feature_cell_matrix.h5 (cell x gene count matrix)
#     - spatial/ directory with tissue images and coordinate files
#
# Output:
#   - Individual Seurat objects for each sample with spatial images attached
#   - QC metrics (mitochondrial %, ribosomal %) added to metadata
#
# Dependencies:
#   - Seurat (>= 5.0)
#   - SeuratObject
#   - Matrix
#   - png
#

# =============================================================================

# -----------------------------------------------------------------------------
# Load Required Packages
# -----------------------------------------------------------------------------

library(Seurat)
library(SeuratObject)
library(Matrix)
library(png)

# Set random seed for reproducibility
set.seed(1234)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Define sample identifiers and sex mapping
# CCAA2905: Male sample, CCAA2105: Female sample
samples <- c("CCAA2905", "CCAA2105")
sample_to_sex <- c("CCAA2905" = "Male", "CCAA2105" = "Female")

# Define base data directory (adjust to your local path)
base_dir <- "/path/to/visium_hd_data"

# Output directory for results
output_dir <- "results"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

#' Canonize Barcode Strings
#'
#' Normalizes cell barcode strings to enable matching between different
#' barcode formats used by 10x Genomics pipelines.
#'
#' @param v Character vector of barcodes
#' @return Character vector of normalized barcodes
#'
#' @details
#' This function handles various barcode formats:
#' - Removes slice/image prefixes (e.g., "slice1_", "image0_")
#' - Extracts core cell identifiers (cellid_X or cell_X patterns)
#' - Removes trailing "-1" suffixes
canonize_barcode <- function(v) {
  v <- tolower(v)
  v <- sub("^(slice[0-9]+_|image[0-9]+_)", "", v)
  v <- ifelse(grepl("cellid_", v, fixed = TRUE),
              sub("^.*?(cellid_.*)$", "\\1", v, perl = TRUE), v)
  v <- ifelse(grepl("cell_", v, fixed = TRUE),
              sub("^.*?(cell_.*)$", "\\1", v, perl = TRUE), v)
  v <- sub("-1$", "", v)
  return(v)
}

#' Attach Spatial Image to Seurat Object with Barcode Alignment
#'
#' Aligns cell barcodes between a Seurat object and its spatial image,
#' then attaches the image to the object. This handles cases where
#' barcode formats differ between count matrix and spatial files.
#'
#' @param obj Seurat object
#' @param img VisiumV2 image object from Read10X_Image()
#' @param verbose Logical, print progress messages
#' @return Seurat object with attached spatial image
#'
#' @details
#' The function:
#' 1. Extracts barcodes from both object and image
#' 2. Normalizes barcodes using canonize_barcode()
#' 3. Identifies common cells between object and image
#' 4. Subsets and reorders both to match
#' 5. Renames cells in object to match image barcodes
#' 6. Attaches the aligned image
attach_image_align <- function(obj, img, verbose = TRUE) {
  stopifnot(inherits(obj, "Seurat"))
  DefaultAssay(obj) <- "Spatial"
  
  # Extract barcodes from object and image
  bc.obj <- colnames(obj)
  bc.img <- Cells(img)
  
  # Normalize barcodes for matching
  key.obj <- canonize_barcode(bc.obj)
  key.img <- canonize_barcode(bc.img)
  common_keys <- intersect(key.obj, key.img)
  
  if (!length(common_keys)) {
    stop("No common cells between object and image after normalization.")
  }
  
  # Map normalized keys back to original barcodes
  first_match <- function(keys, raw, keys_ref) raw[match(keys_ref, keys)]
  obj_keep <- first_match(key.obj, bc.obj, common_keys)
  img_keep <- first_match(key.img, bc.img, common_keys)
  ok <- !is.na(obj_keep) & !is.na(img_keep)
  obj_keep <- obj_keep[ok]
  img_keep <- img_keep[ok]
  
  if (!length(obj_keep)) stop("No common cells after filtering.")
  
  # Subset and align
  obj <- subset(obj, cells = obj_keep)
  img <- subset(img, cells = img_keep)
  ord <- match(img_keep, colnames(obj))
  obj <- obj[, ord]
  
  # Rename cells to match image barcodes and attach image
  obj <- RenameCells(obj, new.names = setNames(img_keep, colnames(obj)))
  obj[["image"]] <- img
  
  if (isTRUE(verbose)) {
    message("Attached image. Common cells: ", ncol(obj))
  }
  
  return(obj)
}

# -----------------------------------------------------------------------------
# Load and Process Individual Samples
# -----------------------------------------------------------------------------

message("\n", rep("=", 80))
message("STEP 1: LOADING VISIUM HD SAMPLES")
message(rep("=", 80))

objects_list <- vector("list", length(samples))

for (i in seq_along(samples)) {
  sample_id <- samples[i]
  message("\nProcessing sample: ", sample_id)
  
  # Define paths
  sample_dir <- file.path(base_dir, sample_id)
  segmented_dir <- file.path(sample_dir, "segmented_outputs")
  spatial_dir <- file.path(segmented_dir, "spatial")
  
  # Create spatial directory if needed
  if (!dir.exists(spatial_dir)) {
    dir.create(spatial_dir, recursive = TRUE)
  }
  
  # -------------------------------------------------------------------------
  # Load count matrix
  # -------------------------------------------------------------------------
  message("  - Loading count matrix...")
  
  counts_matrix <- Read10X_h5(
    file.path(segmented_dir, "filtered_feature_cell_matrix.h5")
  )
  
  # Create Seurat object with basic QC filtering
  # min.cells = 3: Keep genes expressed in at least 3 cells
  # min.features = 200: Keep cells with at least 200 genes detected
  obj <- CreateSeuratObject(
    counts = counts_matrix,
    assay = "Spatial",
    project = "VisiumHD_AdrenalGland",
    min.cells = 3,
    min.features = 200
  )
  
  # -------------------------------------------------------------------------
  # Copy and load spatial files
  # -------------------------------------------------------------------------
  # Copy spatial files from original location to segmented_outputs/spatial
  original_spatial <- file.path(sample_dir, "spatial")
  if (dir.exists(original_spatial)) {
    items <- list.files(original_spatial, all.files = TRUE, full.names = TRUE,
                        recursive = FALSE, no.. = TRUE)
    file.copy(from = items, to = spatial_dir, recursive = TRUE,
              overwrite = FALSE, copy.mode = TRUE, copy.date = TRUE)
  }
  
  # Load spatial image
  message("  - Loading spatial image...")
  img <- Read10X_Image(spatial_dir)
  DefaultAssay(img) <- "Spatial"
  
  # Attach image with barcode alignment
  obj <- attach_image_align(obj, img, verbose = TRUE)
  
  # -------------------------------------------------------------------------
  # Upgrade to high-resolution image if available
  # -------------------------------------------------------------------------
  # The tissue_hires_image.png provides better spatial visualization
  hires_path <- file.path(spatial_dir, "tissue_hires_image.png")
  if (file.exists(hires_path)) {
    img_key <- Images(obj)[1]
    old_img <- obj@images[[img_key]]@image
    old_width <- if (!is.null(dim(old_img))) dim(old_img)[2] else NA_real_
    
    new_img <- png::readPNG(hires_path)
    new_width <- if (!is.null(dim(new_img))) dim(new_img)[2] else NA_real_
    
    # Update image and adjust scale factors proportionally
    if (is.finite(old_width) && is.finite(new_width) && 
        old_width > 0 && new_width > 0) {
      obj@images[[img_key]]@image <- new_img
      scale_factors <- obj@images[[img_key]]@scale.factors
      if (!is.null(scale_factors[["lowres"]])) {
        obj@images[[img_key]]@scale.factors[["lowres"]] <- 
          as.numeric(scale_factors[["lowres"]]) * (new_width / old_width)
      }
      message("  - Upgraded to high-resolution image")
    }
  }
  
  # -------------------------------------------------------------------------
  # Add metadata
  # -------------------------------------------------------------------------
  obj$sample <- sample_id
  obj$sex <- unname(sample_to_sex[sample_id])
  obj$orig.ident <- sample_id
  
  # -------------------------------------------------------------------------
  # Calculate QC metrics
  # -------------------------------------------------------------------------
  # Mitochondrial gene percentage (indicator of cell stress/death)
  # Pattern matches both mouse (mt-) and human (MT-) mitochondrial genes
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^(mt-|MT-)")
  
  # Ribosomal gene percentage (indicator of translation activity)
  # RPS/RPL genes encode ribosomal proteins
  features <- rownames(obj)
  ribo_features <- grep("^RP[SL]", features, value = TRUE, ignore.case = TRUE)
  if (length(ribo_features) > 0) {
    obj[["percent.ribo"]] <- PercentageFeatureSet(obj, features = ribo_features)
  } else {
    obj[["percent.ribo"]] <- 0
  }
  
  # -------------------------------------------------------------------------
  # Print sample summary
  # -------------------------------------------------------------------------
  message("  - Cells: ", ncol(obj))
  message("  - Features (genes): ", nrow(obj))
  message("  - Sex: ", unique(obj$sex))
  message("  - Median MT%: ", round(median(obj$percent.mt), 2))
  message("  - Median Ribo%: ", round(median(obj$percent.ribo), 2))
  
  # Store in list
  objects_list[[i]] <- obj
  names(objects_list)[i] <- sample_id
  
  # Clean up memory
  rm(counts_matrix, img)
  gc()
}

message("\n", rep("=", 80))
message("All samples loaded successfully")
message(rep("=", 80))

# -----------------------------------------------------------------------------
# Save Individual Objects (Optional)
# -----------------------------------------------------------------------------

# Uncomment to save individual sample objects before merging
# for (sample_id in names(objects_list)) {
#   saveRDS(objects_list[[sample_id]], 
#           file.path(output_dir, paste0(sample_id, "_seurat.rds")))
# }

# -----------------------------------------------------------------------------
# Session Info
# -----------------------------------------------------------------------------

message("\nR Session Info:")
print(sessionInfo())