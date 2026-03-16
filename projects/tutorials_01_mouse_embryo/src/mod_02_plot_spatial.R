#' -------------------------------------------------------------------------- '#
#'
#' This script runs preliminary ops for loading data. 
#' Note that we force-setwd to the path in which the source script is contained
current.path <- rstudioapi::getActiveDocumentContext()$path
message(paste("This file is located at:", current.path, 
        ". Will setwd() to its hosting directory"))
setwd(dirname(current.path))

#' Additional custom code is hosted at ../../../src/ . we only source the loader
source(file.path("../../../src", "spatialpkg", "loader.R"))
#' then use loader to load code into a package-like environment and attach it
pkg_env <- load_spatialpkg(attach = TRUE)
#'
#' -------------------------------------------------------------------------- '#



#' Assume mod_01 was run and variables are accesible in environment
#' (sce, spe, obj)

PlotSpatialFeature(object=obj, 
                   feature = "nCount_RNA", 
                   pt_size = 0.5,
                   log = T, 
                   pt_size_prop = T, 
                   conversion_factor = 6, 
                   dark_theme = T,
                   legend_name = "Log of counts")

#' Plotting using sce object
PlotSpatialCat(object = sce, 
               category = "celltype_mapped_refined", 
               pt_size = 0.5,
               order = T,
               pt_size_prop = T,
               conversion_factor = 6,
               dark_theme = F) + NoLegend()

#' Identical result if one uses the Seurat object
PlotSpatialCat(object = obj, 
               category = "celltype_mapped_refined", 
               pt_size = 0.5,
               order = T,
               pt_size_prop = T,
               conversion_factor = 6,
               dark_theme = F) + NoLegend()

#' Restrict plotting to a subset of cell types
PlotSpatialCat(object = spe, 
               category = "celltype_mapped_refined", 
               groups = c("Cranial mesoderm", "Spinal cord", "Cardiomyocytes"),
               pt_size = 0.5,
               order = T,
               pt_size_prop = T,
               conversion_factor = 6,
               dark_theme = F,
               legend_name = "Selected types")
