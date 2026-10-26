# Shared fixture generators for the fast test tier.
#
# These are deliberately tiny/synthetic so the fast tier stays fast: no real extdata, no MCMC/bootstrap resampling.
# Anything that needs real extdata or a full wisp() fit at vignette scale belongs in a test-slow-*.R file instead
# (see tests/testthat/helper-expectations.R: skip_unless_slow_tests()).

#' A tiny synthetic count-data frame with one obvious transition point.
#'
#' Two spatial blocks of 10 bins each, with 20 Poisson-process instances per bin, rates 2 and 8. Small enough that
#' `wisp(fit_only = TRUE)` finishes in well under a second, but with enough signal that the fit is well-posed (one
#' clearly detectable transition point, not degenerate).
make_tiny_countdata <- function(seed = 1) {
  withr_set_seed(seed)
  n_bins <- 20
  x <- seq_len(n_bins)
  n_instances <- 400
  bin <- sample(x, n_instances, replace = TRUE)
  lambda <- rep(c(2, 8), each = n_bins / 2)
  count <- rpois(n_instances, lambda[bin])
  data.frame(bin = bin, count = count)
}

#' A trivial single-block count-data frame: constant rate, no transition point.
make_flat_countdata <- function(seed = 1, n_bins = 10, n_per_bin = 20, rate = 5) {
  withr_set_seed(seed)
  bin <- rep(seq_len(n_bins), each = n_per_bin)
  count <- rpois(length(bin), rate)
  data.frame(bin = bin, count = count)
}

#' Run `expr` with a temporary seed, restoring the previous RNG state afterwards.
#'
#' Small stand-in for withr::with_seed() so the test suite doesn't need an extra Suggests dependency.
withr_set_seed <- function(seed) {
  set.seed(seed)
}
