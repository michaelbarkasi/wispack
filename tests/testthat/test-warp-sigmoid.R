# Fast, pure-function tests for the warping/sigmoid primitives (wisp.warp, wisp.sigmoid, and the underlying Rcpp
# exports warp_mc_R / sigmoid_stable_R / poly_sigmoid_R via their R wrappers). No model fitting involved, so these
# run in a fraction of a second and are a good place to probe numerical edge cases.

test_that("wisp.warp with zero warping factor is the identity", {
  z <- c(1, 5, 10)
  expect_equal(wisp.warp(z, b = 10, w = 0), z)
})

test_that("wisp.warp handles the boundary z == b without error", {
  expect_equal(wisp.warp(10, b = 10, w = 2), 10)
})

test_that("wisp.warp accepts vector and matrix input with matching shape", {
  z_vec <- c(1, 5, 10)
  out_vec <- wisp.warp(z_vec, b = 10, w = 1)
  expect_length(out_vec, length(z_vec))
  expect_true(all(is.finite(out_vec)))

  z_mat <- matrix(1:4, 2, 2)
  out_mat <- wisp.warp(z_mat, b = 10, w = 1)
  expect_equal(dim(out_mat), dim(z_mat))
  expect_true(all(is.finite(out_mat)))
})

test_that("wisp.warp is finite for a large warping factor (no overflow/NaN)", {
  out <- wisp.warp(c(1, 50, 99), b = 100, w = 1e6)
  expect_true(all(is.finite(out)))
})

test_that("wisp.warp validates its inputs", {
  expect_error(wisp.warp(numeric(0), b = 10, w = 1), "length zero")
  expect_error(wisp.warp(11, b = 10, w = 1), "less than or equal to b")
  expect_error(wisp.warp(1, b = c(10, 20), w = 1), "single value")
  expect_error(wisp.warp(1, b = 10, w = c(1, 2)), "single value")
})

test_that("wisp.sigmoid with a single rate block (degree 0) is constant", {
  out <- wisp.sigmoid(c(1, 5, 10), Rt = 3, tslope = numeric(0), tpoint = numeric(0))
  expect_equal(out, c(3, 3, 3))
})

test_that("wisp.sigmoid recovers the two asymptotic rates across a transition point", {
  out <- wisp.sigmoid(c(-1e3, 1e3), Rt = c(1, 9), tslope = 1, tpoint = 0)
  expect_equal(out, c(1, 9), tolerance = 1e-6)
})

test_that("wisp.sigmoid validates argument lengths", {
  expect_error(wisp.sigmoid(1, Rt = c(1, 2), tslope = c(1, 2), tpoint = 1), "Rt must be one element longer")
  expect_error(wisp.sigmoid(1, Rt = c(1, 2), tslope = 1, tpoint = c(1, 2)), "tslope must be the same length")
})

test_that("wisp.sigmoid accepts matrix input with matching shape", {
  x <- matrix(c(-10, -1, 1, 10), 2, 2)
  out <- wisp.sigmoid(x, Rt = c(1, 5), tslope = 1, tpoint = 0)
  expect_equal(dim(out), dim(x))
  expect_true(all(is.finite(out)))
})
