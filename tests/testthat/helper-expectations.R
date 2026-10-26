# Shared helpers for the wispack test suite.

#' Assert that a numeric value falls within an inclusive [lo, hi] range.
#'
#' This is the testthat equivalent of the `if (x < lo || x > hi) stop(...)` idiom used throughout the tutorial_*
#' vignettes. Using one expectation everywhere keeps failure messages uniform and makes the tolerance explicit.
#' `bounds` is typically one of the `c(lo = ..., hi = ...)` entries in `.wisp_golden` (see R/golden_checks.R); pass
#' it as `bounds` or supply `lo`/`hi` directly.
#'
#' @param value Numeric value measured in the test.
#' @param lo Lower bound (inclusive). Ignored if `bounds` is supplied.
#' @param hi Upper bound (inclusive). Ignored if `bounds` is supplied.
#' @param bounds Optional length-2 vector `c(lo = ..., hi = ...)` supplying both bounds at once.
expect_between <- function(value, lo = bounds[["lo"]], hi = bounds[["hi"]], bounds = NULL) {
  label <- deparse(substitute(value))
  testthat::expect(
    isTRUE(value >= lo) && isTRUE(value <= hi),
    sprintf("%s = %s is not between %s and %s.", label, format(value), format(lo), format(hi))
  )
  invisible(value)
}

#' Assert that a numeric value is at least `min_val` (one-sided bound).
expect_at_least <- function(value, min_val) {
  label <- deparse(substitute(value))
  testthat::expect(
    isTRUE(value >= min_val),
    sprintf("%s = %s is not >= %s.", label, format(value), format(min_val))
  )
  invisible(value)
}

#' Skip the current test unless slow/full-scale tests have been explicitly requested.
#'
#' Full-scale tests replicate the tutorial_*.Rmd vignette checks at real vignette scale (real extdata, real
#' bootstrap/MCMC iteration counts). They are the thing that actually guards "does the vignette narrative still
#' match the code's output", but they are too slow to run on every `devtools::test()` / commit. They run when
#' WISPACK_RUN_SLOW_TESTS is set to a truthy value, or whenever NOT_CRAN is unset is irrelevant here since CRAN
#' should never see these either way.
skip_unless_slow_tests <- function() {
  testthat::skip_if_not(
    nzchar(Sys.getenv("WISPACK_RUN_SLOW_TESTS")) &&
      !identical(Sys.getenv("WISPACK_RUN_SLOW_TESTS"), "false") &&
      !identical(Sys.getenv("WISPACK_RUN_SLOW_TESTS"), "0"),
    "slow/full-scale tests skipped (set WISPACK_RUN_SLOW_TESTS=true to run)"
  )
}

#' Skip the current test if a named extdata file isn't present (defensive; the files ship with the package, but
#' this keeps the real-data tests from hard-failing in an environment that trimmed inst/extdata).
skip_unless_extdata <- function(filename) {
  path <- system.file("extdata", filename, package = "wispack")
  testthat::skip_if(identical(path, ""), sprintf("extdata file %s not found", filename))
  invisible(path)
}
