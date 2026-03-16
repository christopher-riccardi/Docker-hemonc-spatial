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

LoadSpatialExpression <- function(cell.coordinates, cell.x.gene) {
  # minimal runtime check for data.table availability
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Please install the 'data.table' package: install.packages('data.table')")
  }
  
  # helper: read input (file path OR data.frame/data.table)
  # uses data.table::fread for file paths (handles .gz automatically)
  dtf <- function(x) {
    if (is.character(x) && length(x) == 1) {
      if (!file.exists(x)) stop("File not found: ", x)
      return(data.table::fread(x, data.table = TRUE, stringsAsFactors = FALSE, showProgress = FALSE))
    }
    if (inherits(x, "data.table")) return(x)
    if (is.data.frame(x)) return(data.table::as.data.table(x))
    stop("Input must be a filepath or a data.frame/data.table")
  }
  
  # read metadata and counts (accepts both file paths and in-memory tables)
  meta_dt <- dtf(cell.coordinates)
  counts_dt <- dtf(cell.x.gene)
  
  diagnostics <- list()
  warnings_out <- character(0)
  
  # preserve any original rownames (fread does not preserve rownames)
  meta_rownames <- rownames(meta_dt)
  counts_rownames <- rownames(counts_dt)
  
  # basic sanity: ensure counts and metadata have same number of rows
  if (nrow(meta_dt) != nrow(counts_dt)) {
    stop("Metadata and counts have different numbers of rows: metadata = ",
         nrow(meta_dt), ", counts = ", nrow(counts_dt),
         ". They must have the same number of entries (one row per cell).")
  }
  
  # basic sanity
  n <- nrow(counts_dt)
  if (n == 0) stop("Counts table contains zero rows")
  
  # helper: determine whether a column is numeric-like
  # treat empty string as NA; consider a column numeric-like when all non-NA entries coerce to numeric
  col_numeric_like <- function(col) {
    ch <- as.character(col)
    ch[ch == ""] <- NA
    non_na <- ch[!is.na(ch)]
    if (length(non_na) == 0) return(TRUE)
    num <- suppressWarnings(as.numeric(non_na))
    !any(is.na(num))
  }
  
  # detect non-numeric columns in counts
  non_num <- vapply(counts_dt, function(c) !col_numeric_like(c), logical(1))
  n_non_num <- sum(non_num)
  
  # if more than one non-numeric column, ambiguous identifier situation -> stop
  if (n_non_num > 1) {
    stop("Counts matrix contains more than one non-numeric column: ",
         paste(names(counts_dt)[non_num], collapse = ", "))
  }
  
  # ---- CASE 1: exactly one non-numeric column -> treat it as identifier column ----
  if (n_non_num == 1) {
    id_col <- names(counts_dt)[which(non_num)]
    ids <- as.character(counts_dt[[id_col]])
    
    # numeric portion (remaining columns)
    nums_dt <- counts_dt[, !non_num, with = FALSE]
    if (ncol(nums_dt) == 0) stop("After removing identifier column, no numeric columns remain in counts")
    
    # coerce numeric portion to numeric matrix, validating coercion for each column
    m <- as.matrix(data.frame(lapply(nums_dt, function(v) {
      vv <- suppressWarnings(as.numeric(as.character(v)))
      if (any(is.na(vv) & !is.na(v))) stop("Non-numeric entries detected in count column: ", deparse(substitute(v)))
      vv
    }), check.names = FALSE, stringsAsFactors = FALSE))
    
    # set rownames on counts using detected ids
    rownames(m) <- ids
    
    diagnostics$id_source <- "counts_identifier_column"
    diagnostics$id_column <- id_col
    
    # attempt to find an exact matching metadata column (same values, same order)
    matched <- NULL
    for (mc in names(meta_dt)) {
      if (identical(as.character(meta_dt[[mc]]), ids)) {
        matched <- mc
        break
      }
    }
    
    # convert metadata to data.frame for returned object
    meta_df <- as.data.frame(meta_dt, stringsAsFactors = FALSE)
    
    if (!is.null(matched)) {
      # assign the matched column as rownames, then remove that column from metadata
      diagnostics$matched_metadata_col <- matched
      rownames(meta_df) <- ids
      meta_df[[matched]] <- NULL
      message("Case: counts contain a single string column '", id_col,
              "' used as cell identifiers; matched metadata column '", matched,
              "' (assigned to rownames and removed from metadata).")
    } else if (!is.null(meta_rownames) && identical(meta_rownames, ids)) {
      diagnostics$matched_metadata_col <- "metadata_rownames"
      rownames(meta_df) <- ids
      # rownames were already present, nothing to remove
      message("Case: counts contain a single string column '", id_col,
              "' used as cell identifiers; matched metadata rownames.")
    } else {
      stop("Identifier column found in counts ('", id_col, "') but no matching metadata column or metadata rownames found in the same order.")
    }
    
    # final verification: metadata and counts must have same number of rows after adjustments
    if (nrow(meta_df) != nrow(m)) {
      stop("After assigning identifiers, metadata and counts have different row counts: metadata = ",
           nrow(meta_df), ", counts = ", nrow(m))
    }
    
    # log shapes and example ids (use utils::head explicitly)
    message("Metadata dimensions: ", nrow(meta_df), " x ", ncol(meta_df))
    message("Counts dimensions: ", nrow(m), " x ", ncol(m))
    message("First 5 cell identifiers: ", paste(utils::head(rownames(m), 5), collapse = ", "))
    
    return(list(metadata = meta_df, counts = m, diagnostics = diagnostics, warnings = warnings_out))
  }
  
  # ---- CASE 2: zero non-numeric columns -> look for progressive integer column (1:n) ----
  int_col <- NULL
  for (nm in names(counts_dt)) {
    v <- suppressWarnings(as.numeric(as.character(counts_dt[[nm]])))
    if (any(is.na(v))) next
    if (all(v == floor(v)) && identical(as.integer(v), as.integer(seq_len(n)))) {
      int_col <- nm
      break
    }
  }
  
  if (!is.null(int_col)) {
    nums_dt <- counts_dt[, setdiff(names(counts_dt), int_col), with = FALSE]
    if (ncol(nums_dt) == 0) stop("After removing progressive id column, no numeric columns remain in counts")
    
    m <- as.matrix(data.frame(lapply(nums_dt, function(v) {
      suppressWarnings(as.numeric(as.character(v)))
    }), check.names = FALSE, stringsAsFactors = FALSE))
    
    ids <- as.character(seq_len(n))
    rownames(m) <- ids
    
    diagnostics$id_source <- "progressive_integer_column"
    diagnostics$id_column <- int_col
    
    meta_df <- as.data.frame(meta_dt, stringsAsFactors = FALSE)
    
    # try to find metadata column that matches 1:n
    matched <- NULL
    for (mc in names(meta_dt)) {
      meta_num <- suppressWarnings(as.numeric(as.character(meta_dt[[mc]])))
      if (!any(is.na(meta_num)) && identical(as.integer(meta_num), as.integer(seq_len(nrow(meta_dt))))) {
        matched <- mc
        break
      }
    }
    
    if (!is.null(matched)) {
      # assign metadata rownames to matched column and remove that column
      diagnostics$matched_metadata_col <- matched
      rownames(meta_df) <- as.character(meta_dt[[matched]])
      meta_df[[matched]] <- NULL
      message("Case: progressive integer column '", int_col, "' used as cell identifiers; matched metadata column '", matched,
              "' (assigned to rownames and removed from metadata).")
    } else {
      diagnostics$matched_metadata_col <- NA
      rownames(meta_df) <- ids
      warnings_out <- c(warnings_out, "Progressive identifier detected in counts but no matching metadata column found; using row numbers.")
      message("Case: progressive integer identifiers detected in counts; no matching metadata column found; using row numbers for metadata.")
    }
    
    # final verification
    if (nrow(meta_df) != nrow(m)) {
      stop("After assigning identifiers, metadata and counts have different row counts: metadata = ",
           nrow(meta_df), ", counts = ", nrow(m))
    }
    
    message("Metadata dimensions: ", nrow(meta_df), " x ", ncol(meta_df))
    message("Counts dimensions: ", nrow(m), " x ", ncol(m))
    message("First 5 cell identifiers: ", paste(utils::head(rownames(m), 5), collapse = ", "))
    
    return(list(metadata = meta_df, counts = m, diagnostics = diagnostics, warnings = warnings_out))
  }
  
  # ---- CASE 3: fallback ----
  # coerce all columns to numeric (error if any non-numeric values present),
  # then try to use metadata or counts rownames; otherwise generate sequential ids.
  
  m <- as.matrix(data.frame(lapply(counts_dt, function(v) {
    vv <- suppressWarnings(as.numeric(as.character(v)))
    if (any(is.na(vv) & !is.na(v))) stop("Non-numeric values detected in counts matrix; if an ID column exists it must be the only non-numeric column")
    vv
  }), check.names = FALSE, stringsAsFactors = FALSE))
  
  meta_df <- as.data.frame(meta_dt, stringsAsFactors = FALSE)
  
  if (!is.null(meta_rownames) && length(meta_rownames) == nrow(m)) {
    rownames(m) <- meta_rownames
    rownames(meta_df) <- meta_rownames
    diagnostics$id_source <- "metadata_rownames"
    message("Case: all count columns numeric; using metadata rownames as cell identifiers.")
  } else if (!is.null(counts_rownames) && length(counts_rownames) == nrow(m)) {
    rownames(m) <- counts_rownames
    rownames(meta_df) <- counts_rownames
    diagnostics$id_source <- "counts_rownames"
    message("Case: all count columns numeric; using counts' existing rownames as cell identifiers.")
  } else {
    ids <- as.character(seq_len(n))
    rownames(m) <- ids
    rownames(meta_df) <- ids
    diagnostics$id_source <- "generated_sequence"
    diagnostics$matched_metadata_col <- NA
    warnings_out <- c(warnings_out, "No explicit IDs detected; generated sequential identifiers. Verify ordering with metadata.")
    message("Case: no explicit identifiers found; generated sequential identifiers.")
  }
  
  # final verification: ensure same number of rows
  if (nrow(meta_df) != nrow(m)) {
    stop("After final alignment, metadata and counts have different row counts: metadata = ",
         nrow(meta_df), ", counts = ", nrow(m))
  }
  
  message("Metadata dimensions: ", nrow(meta_df), " x ", ncol(meta_df))
  message("Counts dimensions: ", nrow(m), " x ", ncol(m))
  message("First 5 cell identifiers: ", paste(utils::head(rownames(m), 5), collapse = ", "))
  
  return(list(metadata = meta_df, counts = m, diagnostics = diagnostics, warnings = warnings_out))
}







