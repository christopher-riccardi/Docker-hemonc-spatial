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

CreateSeurat <- function(counts, metadata, assay = "RNA", project = "SeuratProject") {
  # counts : numeric matrix with rows = cells, cols = genes (assumed validated)
  # metadata: data.frame with rownames set to the same cell identifiers (assumed validated)
  # assay   : name of assay to create in Seurat object
  # project : Seurat project name
  
  # require Seurat package
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required but not installed. Install it with install.packages('Seurat') or appropriate method.")
  }
  
  # prepare assay: Seurat expects genes x cells, so transpose counts
  counts_assay <- t(counts)
  
  # convert metadata to plain data.frame for Seurat meta.data
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  
  # create Seurat object with min.cells = 0 and min.features = 0 to retain all cells and genes
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = counts_assay,
    project = project,
    assay = assay,
    meta.data = metadata,
    min.cells = 0,
    min.features = 0
  )
  
  # logging: shapes and example ids (safe base indexing for first 5)
  n_cells <- if (!is.null(colnames(counts_assay))) ncol(counts_assay) else 0
  n_genes <- nrow(counts_assay)
  first_ids <- character(0)
  if (n_cells > 0) {
    first_ids <- colnames(counts_assay)[seq_len(min(5, n_cells))]
  }
  message("Constructed Seurat object (assay = '", assay, "', project = '", project, "')")
  message("Assay dimensions (genes x cells): ", n_genes, " x ", n_cells)
  if (length(first_ids) > 0) {
    message("First 5 cell identifiers: ", paste(first_ids, collapse = ", "))
  } else {
    message("First 5 cell identifiers: (none available)")
  }
  
  return(seurat_obj)
}