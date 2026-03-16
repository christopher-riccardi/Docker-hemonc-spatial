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



#' When working with seqFISH, we may find the preprocessed data split into two
#' tables: cell information (such as centroid x,y coordinates) and a counts matrix
#' there are referred to as cell coordinates and cell by gene matrix, respectively. 
cell.coordinates <- "../data/Lohoff_et_al_mouse_embryo/embryo1/cell_coordinates.csv.gz"
cell.x.gene <- "../data/Lohoff_et_al_mouse_embryo/embryo1/cell_x_gene.csv.gz"

#' We wrote LoadSpatialExpression, a main input manager that checks inputs 
#' and reports basic information for the data. It returns a list of items: 
#' $metadata, $counts, $diagnostics and $warnings. 
inputs <- LoadSpatialExpression(cell.coordinates = cell.coordinates,
                                cell.x.gene = cell.x.gene)


#' if wanting to work in Bioconductor universe
sce <- CreateExperiment(counts = inputs$counts,
                        metadata = inputs$metadata, 
                        x_col = "x", 
                        y_col = "y",
                        type = "SingleCellExperiment")


spe <- CreateExperiment(counts = inputs$counts,
                        metadata = inputs$metadata, 
                        x_col = "x", 
                        y_col = "y",
                        type = "SpatialExperiment")


#' if wanting to use Seurat instead
obj <- CreateSeurat(counts = inputs$counts,
                    metadata = inputs$metadata, 
                    assay = "RNA", 
                    project = "Embryo")

