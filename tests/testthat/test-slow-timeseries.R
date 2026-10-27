# Port of tutorial_timeseries.Rmd's silent check. Fits a real wisp(fit_only = TRUE) model on real extdata
# (~4-5s). Gated behind WISPACK_RUN_SLOW_TESTS because it's "vignette-scale", even though fit_only keeps any one
# of these individually cheap -- see tests/testthat/helper-expectations.R.

g <- .wisp_golden$timeseries

test_that("wisp() on the ACx laminar timeseries data matches the vignette narrative", {
  skip_unless_slow_tests()
  skip_unless_extdata("ACx_laminar_countdata_demo.csv.gz")
  skip_unless_extdata("ACx_layer_boundary_bins.csv")

  countdata <- read.csv(system.file("extdata", "ACx_laminar_countdata_demo.csv.gz", package = "wispack"))
  data.variables <- list(species = "gene", ran = "mouse", timeseries = "age")
  layer.boundary.bins <- read.csv(system.file("extdata", "ACx_layer_boundary_bins.csv", package = "wispack"))

  model <- wisp(
    count.data     = countdata,
    variables      = data.variables,
    fit_only       = TRUE,
    model.settings = list(LRO_cost = "none", LROwindow_factor = 2.0),
    plot.settings  = list(print.plots = FALSE, dim.bounds = colMeans(layer.boundary.bins)),
    verbose        = FALSE
  )

  gb <- ggplot2::ggplot_build(model$plots$ratecount$plot_pred_context_cortex_fixEff_Nxph3)
  layer_data <- gb$data[[1]]
  test_val <- layer_data[, c("x", "y")][4, 2]
  expect_between(test_val, bounds = g$ratecount_Nxph3)
})
