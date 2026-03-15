# install_packages.R -- installs the reduced package set requested for the hemonc-spatial image
# Sets a reproducible CRAN mirror and suppresses non-fatal remotes warnings that can interrupt automated installs.
options(repos = c(CRAN = "https://cran.rstudio.com"))
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")

# Determine number of physical CPU cores for parallel package builds; fall back to 2 if detection fails.
ncpu_try <- try(parallel::detectCores(logical = FALSE), silent = TRUE)
NCPUS <- if (inherits(ncpu_try, "try-error") || is.null(ncpu_try)) 2 else as.integer(ncpu_try)
if (is.na(NCPUS) || NCPUS < 1) NCPUS <- 2

# Helper that installs only packages not already present in the library.
# This preserves layer caching in Docker builds and avoids reinstalling unchanged packages.
install_if_missing <- function(pkgs, ...) {
  missing <- setdiff(pkgs, rownames(installed.packages()))
  if (length(missing)) install.packages(missing, Ncpus = NCPUS, ...)
}

# Ensure package manager tooling is available: BiocManager for Bioconductor packages and remotes for GitHub installs.
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", Ncpus = NCPUS)
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", Ncpus = NCPUS)

# Preinstall minimal helper packages that are commonly required at build time by other packages (Depends/Imports/Suggests).
# Preinstalling these reduces the chance of build-time failures caused by transitive suggests used during package checks.
preinstall <- c(
  "R.utils", "R6", "Rcpp", "BH", "digest", "glue",
  "pkgconfig", "backports", "rlang", "vctrs", "matrixStats",
  "S4Vectors", "BiocGenerics", "BiocParallel"
)
install_if_missing(preinstall)

# Install the explicitly requested CRAN packages that are missing from the library.
# Listing packages explicitly makes dependency resolution auditable and reproducible.
cran_pkgs <- c(
  "Seurat", "ggplot2", "dplyr", "tidyr", "purrr", "tibble",
  "RColorBrewer", "presto", "msigdbr", "cowplot", "openxlsx",
  "patchwork", "gganimate", "plotly", "RANN"
)
install_if_missing(cran_pkgs)

# Install the specified Bioconductor packages that are not already present.
# BiocManager ensures compatibility with the Bioconductor/R pairing and installs Bioconductor dependencies.
bioc_pkgs <- c(
  "edgeR", "limma", "DESeq2", "SpatialExperiment", "SingleCellExperiment",
  "SummarizedExperiment", "BiocParallel", "GenomicRanges", "scater",
  "scran", "GSVA", "fgsea", "glmGamPoi", "MAST"
)
to_install_bioc <- setdiff(bioc_pkgs, rownames(installed.packages()))
if (length(to_install_bioc)) {
  BiocManager::install(to_install_bioc, ask = FALSE, Ncpus = NCPUS, update = FALSE, installDependencies = TRUE)
}

# Attempt CRAN installation for harmony first; fall back to the canonical GitHub repository if CRAN install is unavailable.
# This covers cases where a required or more recent branch exists only on GitHub.
if (!"harmony" %in% rownames(installed.packages())) {
  try(install.packages("harmony", Ncpus = NCPUS), silent = TRUE)
  if (!"harmony" %in% rownames(installed.packages())) {
    remotes::install_github("immunogenomics/harmony", dependencies = TRUE, upgrade = "never",
                            INSTALL_opts = c("--no-multiarch", "--no-test-load"))
  }
}

# Install ggspavis from GitHub because it is not available on CRAN/Bioconductor.
# INSTALL_opts limit multi-architecture builds and skip test-load to speed installation during image creation.
if (!"ggspavis" %in% rownames(installed.packages())) {
  remotes::install_github("lmweber/ggspavis", dependencies = TRUE, upgrade = "never",
                          INSTALL_opts = c("--no-multiarch", "--no-test-load"))
}

# Ensure presence of fgsea and GSVA; explicitly install Bioconductor variants if missing to avoid namespace/version ambiguity.
if (!"fgsea" %in% rownames(installed.packages())) BiocManager::install("fgsea", ask = FALSE)
if (!"GSVA" %in% rownames(installed.packages())) BiocManager::install("GSVA", ask = FALSE)
