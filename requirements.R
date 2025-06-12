# https://stackoverflow.com/a/38928678

load_libs <- function(bundle = "wispack") {
  if (length(bundle) == 1L && bundle == "wispack") {
    packages <- c(
      "Rcpp",
      "RcppEigen",
      "methods",
      "dtwclust",
      "StanHeaders",
      "BH",
      "ggplot2",
      "grid",
      "gridExtra",
      "colorspace",
      "roxygen2",
      "shiny",
      "shinyjs"
    )
  }

  installed_check <- match(packages, utils::installed.packages()[, 1])

  to_be_installed <- packages[is.na(installed_check)]

  if (length(to_be_installed) > 0L) {
    utils::install.packages(to_be_installed,
      repos = "http://cran.wustl.edu"
    )
  } else {
    print("All required packages already installed")
  }

  for (package in packages) {
    suppressPackageStartupMessages(
      library(package, character.only = TRUE, quietly = TRUE)
    )
  }

}

load_libs("wispack")