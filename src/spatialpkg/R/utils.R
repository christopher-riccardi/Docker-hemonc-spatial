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

#' if in the future there will be more, they can be appended to this vector
.VALID_SINGLECELL_OBJECTS <- c("SingleCellExperiment", "SpatialExperiment", "Seurat")
.DEFAULT_POINT_SHAPE <- 16 #trick: shape 46 is a single pixel

.pull_metadata <- function(object){
  #' pull metadata from object. Object can only be one of:
  #' Seurat, SingleCellExperiment, SpatialExperiment
  object.class <- class(object)
  if(object.class == "SingleCellExperiment" || object.class == "SpatialExperiment"){
    as.data.frame(SingleCellExperiment::colData(object))
  } else if(object.class == "Seurat"){
    as.data.frame(object@meta.data)
  } else {
    stop("Unsupported object type.")
  }
}

.pull_counts <- function(object){
  #' pull counts from object. Object can only be one of:
  #' Seurat, SingleCellExperiment, SpatialExperiment
  #' returns gene x cell counts matrix
  object.class <- class(object)
  if(object.class == "SingleCellExperiment" || object.class == "SpatialExperiment"){
    SingleCellExperiment::counts(object)
  } else if(object.class == "Seurat"){
    object[["RNA"]]$counts
  } else {
    stop("Unsupported object type.")
  }
}

.min_max_scaling <- function(X, a=0, b=1){
  #' assume xmax - xmin is always > 0
  xmin <- min(X); xmax <- max(X)
  a + (X - xmin) * (b - a) / (xmax - xmin)
}


.plot_spatial_property <- function(dataframe,
                                   groups,
                                   pt_size, # point size for plotting (must be > 0)
                                   pt_alpha = 1,
                                   colors = NULL,
                                   conversion_factor=NULL,
                                   dark_theme = TRUE, 
                                   legend_name="Legend"){
  #' must have columns "x", "y", "property"
  #' the four coordinate extremes are computed so that no matter the groups chosen to plot,
  #' they always fall within the context of the whole tissue slide
  #' note that pt_size is a vector of values
  #' offset is mobile, it depends on the y coordinate 
  
  xmin <- min(dataframe$x); xmax <- max(dataframe$x)
  ymin <- min(dataframe$y); ymax <- max(dataframe$y)
  dataframe[["cell.pt.size"]] <- pt_size
  
  offset <- abs(ymax-ymin) / 12
  
  #' note: any additions to the data frame must be made before the following subsetting step
  if (!is.null(groups)) {
    dataframe <- dataframe[dataframe[["property"]] %in% groups, ]
  }
  
  px <- ggplot2::ggplot(dataframe, ggplot2::aes(x = x, y = y, color = property)) +
    ggplot2::geom_point(size = dataframe[["cell.pt.size"]], alpha = pt_alpha, shape = .DEFAULT_POINT_SHAPE) +
    ggplot2::theme_minimal() +
    #ggplot2::xlim(c(xmin - offset, xmax + offset)) +
    #ggplot2::ylim(c(ymin - offset, ymax + offset)) +
    #ggplot2::coord_fixed() +
    ggplot2::coord_fixed(xlim = c(xmin - offset, xmax + offset),
                         ylim = c(ymin - offset, ymax + offset),
                         expand = FALSE) +
    ggplot2::labs(color = legend_name) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(size = 3, alpha = 1)
      )
    )
  
  if (!is.null(colors)) {
    levs <- levels(factor(dataframe$property))
    default_cols <- scales::hue_pal()(length(levs))
    names(default_cols) <- levs
    default_cols[names(colors)] <- unlist(colors)
    px <- px + ggplot2::scale_color_manual(values = default_cols)
  }
  
  if (dark_theme) {
    px <- px + Seurat::DarkTheme()#ggplot2::theme_dark()
  }
  
  if (!is.null(conversion_factor)) {
    #' conversion_factor: microns per pixel
    scale_length_px <- 100 / conversion_factor #microns in pixels
    
    bar_x_start <- xmin# + offset
    bar_x_end <- bar_x_start + scale_length_px
    bar_y <- ymin - (offset * 0.25)
    label_y <- bar_y - (offset * 0.3)
    
    bar_col <- if (isTRUE(dark_theme)) "white" else "black"
    
    px <- px +
      ggplot2::annotate("segment",
                        x = bar_x_start, xend = bar_x_end,
                        y = bar_y, yend = bar_y,
                        linewidth = 0.8,
                        colour = bar_col) +
      ggplot2::annotate("text",
                        x = (bar_x_start + bar_x_end) / 2,
                        y = label_y,
                        label = "100 µm",
                        size = 2,
                        colour = bar_col)
  }
  return(px)
}



PlotSpatialCat <- function(object, 
                           category, 
                           groups = NULL, 
                           pt_size = 1, # point size for plotting (must be > 0)
                           pt_alpha = 1, 
                           colors = NULL,#NULL or named vector/list mapping property values to colors
                           order = TRUE, #order points so that fewer-populated category levels are plotted first
                           pt_size_prop=FALSE, #have the point size be proportional to the log2 of total counts
                           conversion_factor=NULL,
                           dark_theme = TRUE, 
                           legend_name = "Legend") {
  #' Show spatial tissue and project points for a categorical variable.  
  #' x and y columns must be present. The groups option allows to only visualize
  #' category levels of interested, mapped onto the entire tissue centroids.
  valid <- inherits(object, .VALID_SINGLECELL_OBJECTS)
  if (!valid) {
    stop("Unsupported object type.")
  }
  
  dataframe <- .pull_metadata(object = object)
  column.names <- names(dataframe)
  if (!all(c("x", "y", category) %in% column.names)) stop("metadata must contain columns 'x', 'y' and the selected category")
  
  dataframe <- dataframe[, c("x", "y", category)]
  colnames(dataframe)[3] <- "property"
  
  cell.pt.size <- rep(pt_size, nrow(dataframe))
  
  if(pt_size_prop){
    cell.total.counts <- Matrix::colSums(.pull_counts(object))
    cell.total.counts <- log2(1+cell.total.counts)
    #'scale so that the maximum of cell.total.counts equals pt_size
    #' minmax scaling
    pseudocount <- pt_size / 10
    cell.pt.size <- .min_max_scaling(cell.total.counts, 0+pseudocount, pt_size)
  }
  
  if(order){
    levels.frequency <- table(dataframe[["property"]])
    sorted.levels <- names(sort(levels.frequency, decreasing = T))
    dataframe <- dataframe[order(match(dataframe[["property"]], sorted.levels)), ]
  }
  
  px <- .plot_spatial_property(dataframe = dataframe,
                               groups = groups,
                               pt_size = cell.pt.size, # point size for plotting (must be > 0)
                               pt_alpha = pt_alpha,
                               colors = colors,
                               conversion_factor = conversion_factor,
                               dark_theme = dark_theme, 
                               legend_name = legend_name)
  print("Returning")
  return(px)
}


.plot_spatial_feature <- function(dataframe,
                                  pt_size,
                                  pt_alpha = 1,
                                  color = NULL, # palette name or vector of colors
                                  conversion_factor = NULL,
                                  dark_theme = TRUE,
                                  legend_name = NULL) {
  #' must have columns "x", "y", "feature"
  #' the four coordinate extremes are computed so that no matter the groups chosen to plot,
  #' they always fall within the context of the whole tissue slide
  #' note that pt_size is a vector of values
  #' offset is mobile, it depends on the y coordinate 
  
  if (!all(c("x", "y", "feature") %in% colnames(dataframe))) {
    stop("dataframe must contain columns: x, y and feature")
  }
  
  xmin <- min(dataframe$x); xmax <- max(dataframe$x)
  ymin <- min(dataframe$y); ymax <- max(dataframe$y)
  dataframe[["cell.pt.size"]] <- pt_size
  
  offset <- abs(ymax - ymin) / 12
  
  if("retain" %in% colnames(dataframe)){
    #meaning that a filter was applied upstream
    dataframe <- dataframe[dataframe[["retain"]]==T,]
  }
  
  px <- ggplot2::ggplot(dataframe, ggplot2::aes(x = x, y = y, color = feature)) +
    ggplot2::geom_point(size = dataframe[["cell.pt.size"]], alpha = pt_alpha, shape = .DEFAULT_POINT_SHAPE) +
    ggplot2::theme_minimal() +
    #ggplot2::xlim(c(xmin - offset, xmax + offset)) +
    #ggplot2::ylim(c(ymin - offset, ymax + offset)) +
    #ggplot2::coord_fixed() +
    ggplot2::coord_fixed(xlim = c(xmin - offset, xmax + offset),
                         ylim = c(ymin - offset, ymax + offset),
                         expand = FALSE) +
    ggplot2::labs(color = legend_name)# +
  #ggplot2::guides(
  #  color = ggplot2::guide_legend(
  #    override.aes = list(size = 3, alpha = 1)
  #  )
  #)
  
  # color mapping: user-supplied palette name or vector
  if (!is.null(color)) {
    if (is.character(color) && length(color) == 1) {
      pal_name <- tolower(color)
      if (!requireNamespace("viridis", quietly = TRUE)) {
        stop("Please install 'viridis' package or supply a vector of colors in 'color'")
      }
      cols <- switch(pal_name,
                     viridis = viridis::viridis(256),
                     magma   = viridis::magma(256),
                     inferno = viridis::inferno(256),
                     plasma  = viridis::plasma(256),
                     viridis::viridis(256))
      px <- px + ggplot2::scale_color_gradientn(colors = cols)
    } else if (is.character(color) && length(color) >= 2) {
      px <- px + ggplot2::scale_color_gradientn(colors = color)
    } else {
      stop("Invalid 'color' argument. Provide a palette name (e.g. 'viridis') or a vector of colors.")
    }
  } else {
    if (requireNamespace("viridis", quietly = TRUE)) {
      px <- px + ggplot2::scale_color_gradientn(colors = viridis::viridis(256))
    } else {
      px <- px + ggplot2::scale_color_gradient(low = "gray90", high = "red")
    }
  }
  
  px <- px + ggplot2::guides(
    color = ggplot2::guide_colorbar(), 
    size = ggplot2::guide_legend(override.aes = list(size = 3, alpha=1))
  )
  
  if (isTRUE(dark_theme)) {
    if (requireNamespace("Seurat", quietly = TRUE)) {
      px <- px + Seurat::DarkTheme()
    } else {
      px <- px + ggplot2::theme_dark()
    }
  }
  
  if (!is.null(conversion_factor)) {
    #' conversion_factor: microns per pixel
    scale_length_px <- 100 / conversion_factor #microns in pixels
    
    bar_x_start <- xmin# + offset
    bar_x_end <- bar_x_start + scale_length_px
    bar_y <- ymin - (offset * 0.25)
    label_y <- bar_y - (offset * 0.3)
    
    bar_col <- if (isTRUE(dark_theme)) "white" else "black"
    
    px <- px +
      ggplot2::annotate("segment",
                        x = bar_x_start, xend = bar_x_end,
                        y = bar_y, yend = bar_y,
                        linewidth = 0.8,
                        colour = bar_col) +
      ggplot2::annotate("text",
                        x = (bar_x_start + bar_x_end) / 2,
                        y = label_y,
                        label = "100 µm",
                        size = 2,
                        colour = bar_col)
  }
  
  return(px)
}




PlotSpatialFeature <- function(object,
                               feature,
                               filters = NULL,
                               pt_size = 1, # point size for plotting (must be > 0)
                               pt_alpha = 1,
                               log = FALSE,
                               color = NULL, # palette name or vector of colors
                               order = TRUE,
                               pt_size_prop = FALSE,
                               conversion_factor = NULL,
                               dark_theme = TRUE,
                               legend_name = NULL) {
  #' Show spatial tissue and project points for a categorical variable.  
  #' x and y columns must be present.
  #' filters example: filters = "feature > 10 & feature < 40"
  valid <- inherits(object, .VALID_SINGLECELL_OBJECTS)
  if (!valid) {
    stop("Unsupported object type.")
  }
  
  dataframe <- .pull_metadata(object = object)
  column.names <- names(dataframe)
  if (!all(c("x", "y") %in% column.names)) stop("metadata must contain columns 'x' and 'y'")
  
  counts <- NULL
  if (is.null(legend_name)) legend_name <- feature
  
  #' find feature: prefer metadata then counts
  if (feature %in% colnames(dataframe)) {
    dataframe[["feature"]] <- dataframe[[feature]]
  } else {
    counts <- .pull_counts(object)
    if (is.null(rownames(counts)) || !(feature %in% rownames(counts))) {
      stop("Feature '", feature, "' not found in metadata columns nor in counts rownames")
    }
    feature.values <- as.numeric(counts[feature,])
    dataframe[["feature"]] <- feature.values
  }
  
  #' feature has to be numeric. And if the log is requested but feature is less than 0, cannot transform
  if (!is.numeric(dataframe[["feature"]])) stop("Selected feature must be numeric")
  
  if (isTRUE(log)) {
    if (min(dataframe[["feature"]]) < 0){
      warning("Cannot log-transform feature, negative values present")
    }else{
      dataframe[["feature"]] <- log1p(dataframe[["feature"]])
    }
  }
  
  #' apply filters (expression or string) evaluated in dataframe environment
  if (!is.null(filters)) {
    fexpr <- if (is.character(filters)) parse(text = filters)[[1]] else substitute(filters)
    mask <- eval(fexpr, envir = dataframe, enclos = parent.frame())
    mask <- as.logical(mask)
    dataframe[["retain"]] <- mask
    if (sum(mask) == 0) stop("No cells remain after filtering")
  }
  
  cell.pt.size <- rep(pt_size, nrow(dataframe))
  
  if(pt_size_prop){
    cell.total.counts <- Matrix::colSums(.pull_counts(object))
    cell.total.counts <- log2(1+cell.total.counts)
    #'scale so that the maximum of cell.total.counts equals pt_size
    #' minmax scaling
    pseudocount <- pt_size / 10
    cell.pt.size <- .min_max_scaling(cell.total.counts, 0+pseudocount, pt_size)
  }
  
  #' set ordering so that feature with higher values are plotted last (i.e. on top)
  if (isTRUE(order)) {
    ord <- order(dataframe[["feature"]], decreasing = F)
    dataframe <- dataframe[ord, ]
    cell.pt.size <- cell.pt.size[ord]
  }
  
  px <- .plot_spatial_feature(dataframe = dataframe,
                              pt_size = cell.pt.size,
                              pt_alpha = pt_alpha,
                              color = color,
                              conversion_factor = conversion_factor,
                              dark_theme = dark_theme,
                              legend_name = legend_name)
  
  return(px)
}


ApproxConversionFactor <- function(object, zlim=1, assumed=10){
  #' an avg cell is X um wide. Distance between two centroids
  #' is approximately X um. Median distance between 
  #' high-quality cells and any surrounding ones should be X.
  #' assume X=10 um, X should tell the scaling factor for 10 um.
  #' This is an approximation that does not take composition into acct
  #' zlim is the Z score cutoff to pick reference cells. approx is the
  #' assumed cell centroid to cell centroid distance.

  object.class <- class(object)
  if(object.class == "SingleCellExperiment" || object.class == "SpatialExperiment"){
    total.cts <- log1p(SingleCellExperiment::colData(object)$sum)
    detected.trx <- log1p(SingleCellExperiment::colData(object)$detected)
  } else if(object.class == "Seurat"){
    total.cts <- log1p(object@meta.data$nCount_RNA)
    detected.trx <- log1p(object@meta.data$nFeature_RNA)
  } else {
    stop("Unsupported object type.")
  }
  
  x.total.cts <- mean(total.cts)
  sd.total.cts <- stats::sd(total.cts)
  z.total.cts <- (total.cts-x.total.cts)/sd.total.cts
  
  x.detected.trx <- mean(detected.trx)
  sd.detected.trx <- stats::sd(detected.trx)  
  z.detected.trx <- (detected.trx-x.detected.trx)/sd.detected.trx
  
  #' reference cells set composed of cells with good stats
  reference <- -zlim < z.total.cts & z.total.cts < zlim & 
    -zlim < z.detected.trx & z.detected.trx < zlim
  
  message("Using ", sum(reference), " cells as reference")
  
  coords <- as.data.frame(.pull_metadata(object)[c("x", "y")])
  ref.data <- cbind(coords, reference=reference)
  
  nn.res <- RANN::nn2(coords, k = 2) 
  
  closest.distances <- sapply(which(ref.data$reference), function(i) {
    nn_idx <- nn.res$nn.idx[i, 2]
    nn_dist <- nn.res$nn.dists[i, 2]
    return(nn_dist)
  })
  assumed/stats::median(closest.distances)
}




SalvadorChan <- function(y, min_frac = 0.25) {
  "
  This function applies the L method (Salvador and Chan, 2004) for the estimation 
  of the appropriate number of classes on an evaluation graph.
  Readapted from https://www.mathworks.com/matlabcentral/fileexchange/38771-l-method
  https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/38771/versions/2/screenshot.jpg
  but added a minimum segment length for robustness
  "
  x <- seq_along(y)
  
  n <- length(y)
  min_seg <- max(2, floor(min_frac * n))
  
  ks <- min_seg:(n - min_seg)
  
  total_sse <- sapply(ks, function(k) {
    fit_l <- stats::lm(y[1:k] ~ x[1:k])
    fit_r <- stats::lm(y[(k + 1):n] ~ x[(k + 1):n])
    
    sum(stats::resid(fit_l)^2) + sum(stats::resid(fit_r)^2)
  })
  
  k_star <- ks[which.min(total_sse)]
  
  list(
    elbow_index = k_star,
    elbow_x = x[k_star],
    total_sse = total_sse,
    ks = ks
  )
}

.summarize_pair <- function(pp, mark_i, mark_j, nsim, min.r=100) {
  
  "
  EDIT
  
  summarize the spatial association between two cell populations
  in a multitype point pattern using a Monte Carlo test based on the cross
  pair-correlation function pcfcross
  
  This way we test the H0 that the two cell types are spatially independent, 
  conditional on the observed cell locations (and this is done by 
  repeatedly randomizing cell-type labels and comparing the observed cross 
  pair-correlation function to its distribution under the null
  model)
  For each distance in r.vec the observed pair-correlation value is classified as 
  showing attraction above the global simulation envelope, repulsion (below) 
  or no detectable deviation from independence.
  then we aggregate these distance-wise classifications into a single
  signed summary statistic, defined as the fraction of distances exhibiting
  significant attraction minus the fraction exhibiting significant repulsion.
  
  Therefore the final summarized vale quantifies the dominant spatial trend between the
  two cell populations across the entire distance interval: positive values
  indicate overall attraction, negative values indicate overall repulsion, and
  values near zero indicate no consistent spatial association.
  "
  
  env <- spatstat.explore::envelope(
    pp,
    spatstat.explore::pcfcross.inhom,
    i = mark_i,
    j = mark_j,
    simulate = expression(spatstat.random::rlabel(pp)),
    nsim=nsim,
    global=F
  )
  
  df <- data.frame(env)
  df <- df[df$r > min.r,]
  df$signif <- TRUE
  df[which(df$obs>df$lo & df$obs < df$hi), "signif"] <- FALSE
  df$enrich <- ifelse(df$obs>1, 1, -1)
  sum(df[df$signif==T, "enrich"])/nrow(df)
}
