# Port of tutorial_corticallaminar.Rmd's silent checks 4-5 (the ones that need a real, default-settings wisp()
# fit -- MCMC with the default 1100 steps -- on real extdata, ~30s). Checks 1-3 (data-only) are in
# test-corticallaminar-data.R (fast tier).

g <- .wisp_golden$corticallaminar

test_that("wisp() default fit on the S1 laminar data matches the vignette narrative", {
  skip_unless_slow_tests()
  skip_unless_extdata("S1_laminar_countdata_demo.csv")
  skip_unless_extdata("S1_layer_boundary_bins.csv")

  countdata <- read.csv(system.file("extdata", "S1_laminar_countdata_demo.csv", package = "wispack"))
  data.variables <- list(species = "gene", ran = "mouse", timeseries = "age")
  layer.boundary.bins <- read.csv(system.file("extdata", "S1_layer_boundary_bins.csv", package = "wispack"))

  set.seed(123)
  model <- wisp(
    count.data     = countdata,
    variables      = data.variables,
    model.settings = list(LRO_cost = "none"),
    plot.settings  = list(print.plots = FALSE, dim.bounds = colMeans(layer.boundary.bins)),
    verbose        = FALSE
  )

  # check4
  expect_between(model$stats$parameter[1, 2], bounds = g$param1_estimate)

  # check5
  mask <- grepl("Rorb", model$stats$parameters$parameter) &
    grepl("right", model$stats$parameters$parameter) &
    grepl("beta", model$stats$parameters$parameter) &
    grepl("Rt", model$stats$parameters$parameter)
  test_val <- model$stats$parameter[mask, ][2, 3]
  expect_between(test_val, bounds = g$hemisphere_RORB)
})
