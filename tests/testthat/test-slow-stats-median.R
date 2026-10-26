# Port of tutorial_stats.Rmd's check 5, the one check in that vignette which needs a real wisp() fit
# (use.median = TRUE, bootstraps.num = 0, but the default MCMC still runs) on real extdata. Checks 1-4 and 6-8
# only need the pre-fit fixture model and live in test-wisp-stats-fixture.R (fast tier).

g <- .wisp_golden$stats

test_that("wisp() with use.median = TRUE on the S1 laminar data matches the vignette narrative", {
  skip_unless_slow_tests()
  skip_unless_extdata("S1_laminar_countdata_demo.csv")

  countdata <- read.csv(system.file("extdata", "S1_laminar_countdata_demo.csv", package = "wispack"))
  data.variables <- list(species = "gene", ran = "mouse", timeseries = "age")

  model_median <- wisp(
    count.data     = countdata,
    variables      = data.variables,
    use.median     = TRUE,
    model.settings = list(LRO_cost = "none"),
    plot.settings  = list(print.plots = FALSE),
    bootstraps.num = 0,
    verbose        = FALSE
  )

  expect_between(model_median$stats$parameters[1, 2], bounds = g$median_fit_param1)
})
