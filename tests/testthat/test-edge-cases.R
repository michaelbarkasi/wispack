# Degenerate/boundary-input tests for wisp(). All of these use tiny synthetic data and fit_only = TRUE, so the
# whole file runs in a fraction of a second.

test_that("wisp() handles data with very few spatial bins without erroring", {
  set.seed(1)
  cd <- data.frame(bin = rep(1:3, each = 20), count = rpois(60, 5))
  model <- wisp(cd, fit_only = TRUE, model.settings = list(LRO_cost = "none"),
                plot.settings = list(print.plots = FALSE), verbose = FALSE)
  expect_false(anyNA(model$fitted.parameters))
  expect_true(all(is.finite(model$fitted.parameters)))
})

test_that("wisp() handles all-zero counts without erroring", {
  cd <- data.frame(bin = rep(1:10, each = 10), count = rep(0, 100))
  model <- wisp(cd, fit_only = TRUE, model.settings = list(LRO_cost = "none"),
                plot.settings = list(print.plots = FALSE), verbose = FALSE)
  expect_false(anyNA(model$fitted.parameters))
  expect_true(all(is.finite(model$fitted.parameters)))
})

test_that("wisp() on constant-rate data recovers (near-)equal rates across any detected blocks", {
  # The LRO change-point search can still flag a spurious transition point on flat data (it's a heuristic over
  # noisy counts), but the recovered rates for that "transition" should stay close to each other -- this is a
  # sanity check on the optimizer, not a claim that LRO always finds exactly one block.
  set.seed(1)
  cd <- make_flat_countdata(seed = 1, n_bins = 10, n_per_bin = 50, rate = 5)
  model <- wisp(cd, fit_only = TRUE, model.settings = list(LRO_cost = "none"),
                plot.settings = list(print.plots = FALSE), verbose = FALSE)
  rts <- model$fitted.parameters[grepl("Rt", model$param.names)]
  expect_true(all(is.finite(rts)))
  expect_lt(diff(range(rts)) / mean(rts), 0.1)
})

test_that("wisp.sigmoid handles a degree-0 (no transition point) model for every input shape", {
  z <- c(1, 50, 99)
  out <- wisp.sigmoid(z, Rt = 4, tslope = numeric(0), tpoint = numeric(0))
  expect_equal(out, rep(4, length(z)))

  out_mat <- wisp.sigmoid(matrix(z, nrow = 1), Rt = 4, tslope = numeric(0), tpoint = numeric(0))
  expect_equal(out_mat, matrix(rep(4, length(z)), nrow = 1))
})
