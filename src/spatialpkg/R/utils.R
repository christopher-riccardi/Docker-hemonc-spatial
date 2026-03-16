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

.filetype <- function(path) {
  f <- file(path)
  ext <- summary(f)$class
  close.connection(f)
  ext
}
