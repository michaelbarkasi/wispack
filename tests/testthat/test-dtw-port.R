# Regression test for project_cp()'s DTW alignment, ported from R (which called dtw::dtw(...,
# step.pattern = dtw::symmetricP05)) to native C++ (src/model_fit.cpp's dtw_symmetricP05_align() and
# project_cp()) so the compiled model-fitting code no longer calls back into R, and project_cp need
# not be an exported R function (see NAMESPACE's note on why dtw was dropped from Imports).
#
# These expected values were computed once, before the port, by running the *original* R
# implementation (dtw::dtw(Xi, centroid, step.pattern = dtw::symmetricP05, distance.only = FALSE)) on
# each case below -- not re-derived from the new port, so this is a genuine check against the
# reference algorithm, not a tautology. The full validation (1000+ random/edge cases against the live
# dtw package) was run once during development and is not repeated here, since dtw is no longer a
# package dependency; these few hand-picked cases are the permanent, dependency-free regression guard.
#
# .project_cp_test_hook() is an internal, non-exported Rcpp wrapper around project_cp() (see
# src/model_fit.cpp), accessed here via wispack:::.

test_that("project_cp()'s C++ DTW port matches the original dtw::dtw(step.pattern = symmetricP05) oracle", {
  hook <- wispack:::.project_cp_test_hook

  expect_equal(
    hook(c(2, 4), c(1, 2, 3, 4, 5), matrix(c(1, 2, 3, 4, 5), ncol = 1)),
    matrix(c(2, 4), ncol = 1)
  )

  expect_equal(
    hook(c(1, 3, 5), c(5, 1, 9, 2, 8, 3), matrix(c(5, 1, 9, 2, 8, 3), ncol = 1)),
    matrix(c(1, 3, 5), ncol = 1)
  )

  expect_equal(
    hook(
      c(2, 5, 7),
      c(1, 1, 5, 5, 1, 1, 5, 5),
      cbind(c(1, 1, 1, 5, 5, 1, 1, 5), c(5, 1, 1, 5, 5, 1, 1, 1))
    ),
    matrix(c(2, 6, 8, 2, 6, 8), ncol = 2)
  )

  # NA in X_lik (project_cp() replaces it with 1 before alignment, matching the original R function)
  expect_equal(
    hook(c(3, 8), c(0, 0, 0, 10, 10, 10, 0, 0, 0, 0), matrix(c(0, 0, NA, 10, 10, 9, 0, 0, 1, 0), ncol = 1)),
    matrix(c(3, 10), ncol = 1)
  )
})

test_that("project_cp()'s C++ port returns a 1x1 matrix of -1 for degenerate inputs", {
  hook <- wispack:::.project_cp_test_hook

  expect_equal(hook(integer(0), c(1, 2, 3), matrix(1:3, ncol = 1)), matrix(-1L, 1, 1))
  expect_equal(hook(c(1L), numeric(0), matrix(1:3, ncol = 1)), matrix(-1L, 1, 1))
  expect_equal(hook(c(1L), c(1, 2, 3), matrix(numeric(0), nrow = 3, ncol = 0)), matrix(-1L, 1, 1))
})
