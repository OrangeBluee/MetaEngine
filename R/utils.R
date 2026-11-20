# MetaEngine utilities.

pm_require <- function(pkg,
                       bioc = FALSE,
                       install = FALSE) {

  if (requireNamespace(pkg, quietly = TRUE)) {
    return(invisible(TRUE))
  }

  if (!install) {
    stop("Required package missing: ", pkg, call. = FALSE)
  }

  message("MetaEngine: installing package '", pkg, "' ...")

  if (bioc) {

    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }

    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE
    )

  } else {

    install.packages(
      pkg,
      dependencies = TRUE
    )
  }

  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("MetaEngine could not install required package: ", pkg, call. = FALSE)
  }

  message("MetaEngine: package '", pkg, "' is ready.")
  invisible(TRUE)
}

pm_initialize <- function(install = FALSE) {
  
  message("Initializing MetaEngine...")
  
  # CRAN packages
  cran_packages <- c(
    "data.table",
    "dplyr",
    "ggplot2",
    "igraph",
    "readxl",
    "xml2"
  )
  
  # Bioconductor packages
  bioc_packages <- c(
    "limma",
    "preprocessCore"
  )
  
  for (p in cran_packages) {
    pm_require(
      pkg = p,
      bioc = FALSE, install = install
    )
  }
  
  for (p in bioc_packages) {
    pm_require(
      pkg = p,
      bioc = TRUE, install = install
    )
  }
  
  message("MetaEngine initialization complete.")

  # Check adduct-expansion dimensions during initialization.
  if (exists("pm_test_adduct_expansion", mode="function")) {
    pm_test_adduct_expansion()
  }

  invisible(TRUE)
}

pm_bh <- function(p) {
  z <- rep(NA_real_, length(p))
  ok <- is.finite(p)
  z[ok] <- p.adjust(p[ok], "BH")
  z
}

pm_read <- function(x) {
  if (is.data.frame(x)) return(x)
  if (!is.character(x) || length(x) != 1L) stop("Input must be a data.frame or file path.")
  if (!file.exists(x)) stop("File not found: ", x)
  if (tolower(tools::file_ext(x)) == "csv")
    read.csv(x, check.names=FALSE, stringsAsFactors=FALSE)
  else
    read.delim(x, check.names=FALSE, stringsAsFactors=FALSE)
}


`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

