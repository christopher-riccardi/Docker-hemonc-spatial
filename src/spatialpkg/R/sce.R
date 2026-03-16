#' Spatial Transcriptomics Utility Functions
#'
#' Project: spatialpkg
#' Authors: Christopher Riccardi,
#' Institution: Children's Hospital Los Angeles (CHLA),
#'              in collaboration with University of Southern California (USC)
#'
#' Description:
#' This module provides helper functions for constructing and manipulating
#' spatial transcriptomics data objects used throughout analyses.
#' The code is designed to integrate seqFISH-derived measurements with
#' SingleCellExperiment- and Seurat-based workflows, enabling unified preprocessing,
#' metadata harmonization, and downstream spatial analyses.
#'
#' The functions defined here are intended for internal use within the
#' spatialpkg project library and are loaded dynamically via the
#' package-style loader (`loader()`).
#'
#' Notes:
#' - Designed for large-scale spatial transcriptomics datasets.
#' - Compatible with seqFISH, MERFISH, and related imaging-based modalities.
#' - Intended for research use.
#'
#' Repository:
#' 
#' https://github
#'
#' Last Updated: 2026-03-15
#'
NULL

CreateExperiment <- function(counts, metadata,
                             x_col = "x", y_col = "y",
                             type = c("SpatialExperiment", "SingleCellExperiment")) {
  # counts : numeric matrix with rows = cells, cols = genes (assumed validated)
  # metadata: data.frame with rownames set to the same cell identifiers (assumed validated)
  # x_col, y_col : names of spatial columns in metadata (used only for SpatialExperiment)
  # type : which class to construct
  
  type <- match.arg(type)
  
  # check required Bioconductor packages are available (fail early if missing)
  if (type == "SingleCellExperiment" && !requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Package 'SingleCellExperiment' is required but not installed")
  }
  if (type == "SpatialExperiment" && !requireNamespace("SpatialExperiment", quietly = TRUE)) {
    stop("Package 'SpatialExperiment' is required but not installed")
  }
  if (!requireNamespace("S4Vectors", quietly = TRUE)) {
    stop("Package 'S4Vectors' is required but not installed")
  }
  
  # --- prepare assay: SingleCellExperiment / SpatialExperiment expect genes x cells
  # assume counts rows = cells, cols = genes (as produced by your loader)
  # transpose once to produce genes x cells
  counts_assay <- t(counts)
  
  # --- prepare colData: convert metadata to S4Vectors::DataFrame using existing rownames
  colData <- S4Vectors::DataFrame(metadata, row.names = rownames(metadata))
  
  # --- construct requested object
  if (type == "SingleCellExperiment") {
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = counts_assay),
      colData = colData
    )
    
    # logging: shapes and example ids
    message("Constructed SingleCellExperiment")
    message("Assay dimensions (genes x cells): ", nrow(counts_assay), " x ", ncol(counts_assay))
    message("First 5 cell identifiers: ", paste(utils::head(colnames(counts_assay), 5), collapse = ", "))
    
    return(sce)
  }
  
  # SpatialExperiment branch
  # assume x_col and y_col exist in metadata and are numeric/coercible (validated upstream)
  coords <- as.matrix(metadata[, c(x_col, y_col), drop = FALSE])
  rownames(coords) <- rownames(metadata)
  
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts_assay),
    colData = colData,
    spatialCoords = coords
  )
  
  # logging: shapes, coords used, example ids
  message("Constructed SpatialExperiment")
  message("Assay dimensions (genes x cells): ", nrow(counts_assay), " x ", ncol(counts_assay))
  message("Spatial coordinate columns used: ", x_col, ", ", y_col)
  message("First 5 cell identifiers: ", paste(utils::head(colnames(counts_assay), 5), collapse = ", "))
  
  return(spe)
}
