# Mechanical exercise of the plotting/stats machinery (plot.ratecount, plot.decomposition, plot.timeseries,
# plot.parameters, plot.MCMC.walks, plot.MCMC.bs.comparison, sample.stats()) on a small synthetic model with
# more than one fitted parameter and a handful of bootstrap resamples -- see
# tests/testthat/helper-fixtures.R's make_small_multiparam_countdata(). This is decoupled from "does a specific
# RORB story reproduce" (that's the slow-tier real-data tests and the tutorial vignettes); it only checks that
# each function runs without error and returns the expected shape/structure, not specific values.
#
# A single-fitted-parameter model hits real, previously-fixed degenerate edge cases in this package, so the
# fixture deliberately has several (14 at the default sizing). Fitting the model once here and reusing it
# across test_that() blocks keeps the file fast (a few seconds total).

model <- wisp(
  count.data     = make_small_multiparam_countdata(),
  variables      = list(timeseries = "age"),
  model.settings = list(LRO_cost = "none"),
  plot.settings  = list(print.plots = FALSE),
  MCMC.settings  = list(MCMC.steps = 100, MCMC.burnin = 20),
  bootstraps.num = 10,
  max.fork       = 1,
  verbose        = FALSE
)

test_that("the synthetic fixture has more than one fitted parameter", {
  expect_gt(length(model$fitted.parameters), 1)
  expect_true(all(is.finite(model$fitted.parameters)))
})

test_that("plot.ratecount runs on the synthetic fixture", {
  rc <- plot.ratecount(model, species = "demoGene", verbose = FALSE)
  expect_true(is.list(rc))
  expect_true("plot_pred_context_demoContext_fixEff_demoGene" %in% names(rc))
  expect_s3_class(rc$plot_pred_context_demoContext_fixEff_demoGene, "ggplot")
})

test_that("plot.decomposition runs on the synthetic fixture", {
  decomp <- plot.decomposition(model, species = "demoGene")
  expect_true(is.list(decomp))
  expect_gt(length(decomp), 0)
})

test_that("plot.timeseries runs on the synthetic fixture", {
  ts <- plot.timeseries(model, species = "demoGene", verbose = FALSE)
  expect_true(is.list(ts))
  expect_true("plot_pred_context_demoContext_timeseries_demoGene" %in% names(ts))
  expect_s3_class(ts$plot_pred_context_demoContext_timeseries_demoGene, "ggplot")
})

test_that("plot.parameters runs on the synthetic fixture", {
  pp <- plot.parameters(model, species = "demoGene", verbose = FALSE)
  expect_true(is.list(pp))
  expect_true("plot_baseline_demoContext_demoGene" %in% names(pp))
  expect_true("plot_ranEff_demoContext_demoGene" %in% names(pp))
})

test_that("plot.MCMC.walks runs on the synthetic fixture", {
  mw <- plot.MCMC.walks(model, print.plots = FALSE, verbose = FALSE, low_samples = 5)
  expect_true(is.list(mw))
  expect_true(all(c("plot.walks.parameters_low", "plot.walks.parameters_high", "plot.walks.nll") %in% names(mw)))
})

test_that("plot.MCMC.bs.comparison runs on the synthetic fixture (needs >1 parameter to not crash)", {
  bc <- plot.MCMC.bs.comparison(model, verbose = FALSE)
  expect_true(is.list(bc))
  expect_true(all(c("plot_sample_correlations", "plot_comparison_Shaprio", "plot_comparison_density") %in% names(bc)))
  expect_s3_class(bc$plot_comparison_density, "ggplot")
})

test_that("sample.stats runs on the synthetic fixture", {
  ss <- sample.stats(model, verbose = FALSE)
  expect_true(is.data.frame(ss$stats$parameters))
  expect_gt(nrow(ss$stats$parameters), 1)
  expect_true(all(c("parameter", "estimate", "CI.low", "CI.high", "p.value") %in% names(ss$stats$parameters)))
})
