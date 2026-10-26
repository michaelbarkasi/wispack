# Tests for wisp()'s argument-checking contract (the internal check_list() helper applied to `variables`,
# `model.settings`, and `MCMC.settings`), plus the equivalent checks in attractor_simulation() and
# create.count.data.WSPmm(). None of this is exercised by the ported vignette checks, but it's exactly the kind
# of thing that silently rots: nobody notices an error message changed, or a validation path got skipped, until a
# user hits it. All of these use tiny synthetic data and fit_only = TRUE (or error before any fitting happens), so
# the whole file runs in well under a second.

tiny_countdata <- function() make_tiny_countdata(seed = 1)

test_that("wisp() validates the `variables` argument", {
  cd <- tiny_countdata()
  expect_error(wisp(cd, variables = "not a list", fit_only = TRUE, verbose = FALSE), "must be a list")
  expect_error(wisp(cd, variables = list(1, 2), fit_only = TRUE, verbose = FALSE), "named elements")
  expect_error(wisp(cd, variables = list(not_a_role = "bin"), fit_only = TRUE, verbose = FALSE), "can contain only")
  expect_error(wisp(cd, variables = list(bin = 5), fit_only = TRUE, verbose = FALSE), "must be of type character")
})

test_that("wisp() validates the `model.settings` argument", {
  cd <- tiny_countdata()
  expect_error(wisp(cd, model.settings = "not a list", fit_only = TRUE, verbose = FALSE), "must be a list")
  expect_error(wisp(cd, model.settings = list(not_a_setting = 1), fit_only = TRUE, verbose = FALSE), "can contain only")
  expect_error(wisp(cd, model.settings = list(ctol = "not numeric"), fit_only = TRUE, verbose = FALSE), "must be of type numeric")
})

test_that("wisp() validates the `MCMC.settings` argument (only reached when fit_only = FALSE)", {
  cd <- tiny_countdata()
  # fit_only = TRUE skips MCMC.settings validation entirely (MCMC never runs), so these need fit_only = FALSE.
  # The invalid-argument checks below error out of wisp() before any fitting happens, so they stay fast.
  expect_error(wisp(cd, MCMC.settings = list(not_a_setting = 1), bootstraps.num = 0, verbose = FALSE), "can contain only")
  expect_error(wisp(cd, MCMC.settings = list(MCMC.burnin = -1), bootstraps.num = 0, verbose = FALSE), "must be >= 0")
})

test_that("wisp() validates the `plot.settings` argument", {
  cd <- tiny_countdata()
  expect_error(wisp(cd, plot.settings = "not a list", fit_only = TRUE, verbose = FALSE), "must be a list")
})

test_that("attractor_simulation() validates its seed_data columns", {
  bad_seed <- data.frame(x = 1:5, y = 1:5)
  expect_error(attractor_simulation(bad_seed), "must contain columns")
})

test_that("create.count.data.WSPmm() validates bin.dim and required columns", {
  df <- data.frame(mouse = 1, cell_num = 1, trt = "a", x_bins = 1, y_bins = 1)
  expect_error(
    create.count.data.WSPmm(df, bin.dim = "z_bins", gene.list = "foo", fixed.effect.names = "trt"),
    "bin.dim must be one of"
  )
  expect_error(
    create.count.data.WSPmm(df[, setdiff(names(df), "mouse")], bin.dim = "x_bins", gene.list = "foo", fixed.effect.names = "trt"),
    "missing mouse, cell_num, or fixed effect column"
  )
})
