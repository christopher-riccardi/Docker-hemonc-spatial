#'Loader for spatialpkg
#'detects the caller file (when sourced from a script) and resolves package path

this.file <- function() {
  for (i in rev(seq_along(sys.frames()))) {
    if (!is.null(sys.frames()[[i]]$ofile)) {
      return(sys.frames()[[i]]$ofile)
    }
  }
}

current.file <- normalizePath(this.file())

load_spatialpkg <- function(pkg_dir = NULL, attach = TRUE, envir = NULL) {

  caller <- current.file
  pkg_dir <- dirname(caller)

  if (!dir.exists(pkg_dir)) stop("Package directory not found: ", pkg_dir)
  r_dir <- file.path(pkg_dir, "R")
  if (!dir.exists(r_dir)) stop("R/ directory not found inside package: ", r_dir)
  files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
  files <- sort(files)

  # create package environment
  pkg_env <- new.env(parent = baseenv())

  # record registry to avoid re-sourcing if called multiple times
  if (!exists('.spatialpkg_sourced', envir = pkg_env, inherits = FALSE)) {
    assign('.spatialpkg_sourced', character(), envir = pkg_env)
  }

  for (f in files) {
    nf <- normalizePath(f)
    if (nf %in% get('.spatialpkg_sourced', envir = pkg_env)) next
    assign('.spatialpkg_sourced', c(get('.spatialpkg_sourced', envir = pkg_env), nf), envir = pkg_env)
    sys.source(f, envir = pkg_env)
  }

  if (attach) {
    # attach the package environment so functions are available like a package
    name <- paste0('package:spatialpkg')
    if (name %in% search()) detach(name, character.only = TRUE)
    attach(pkg_env, name = name)
  }

  invisible(pkg_env)
}
