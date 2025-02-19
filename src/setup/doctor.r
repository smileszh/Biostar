#!/usr/bin/env Rscript

# Bioconductor packages used in the script.
bio.p <- c("DESeq2", "gplots", "biomaRt",
           "tibble", "dplyr", "tools"
)

cat("# Doctor, Doctor! Give me the R news!\n")


check_packages <- function(packages) {
  errs = 0

  for (package in packages) {
    fmtname <- sprintf("%-10s", substr(package, 1, 10))

    cat("# Checking", fmtname, " ... ")
    # Try to load the package
    result <- tryCatch({
      suppressPackageStartupMessages(
        library(package, character.only = TRUE)
      )
      cat("OK\n")
    }, error = function(e) {
      cat("ERROR!\n")
      errs <<- errs + 1
    })

  }

  if (errs > 0) {
    cat("# Errors found", errs, "packages are missing\n")
  } else {
    cat("# You are doing well, Majesty!\n", errs)
  }

}

check_packages(bio.p)



