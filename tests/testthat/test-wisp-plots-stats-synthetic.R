# Mechanical exercise of the plotting/stats machinery (plot(model, type = "ratecount"[, decompose = TRUE]/
# "timeseries"/"parameters"/"MCMC.walks"/"MCMC.bs.comparison"), sample.stats(), and the print/summary/residuals/
# coef/confint S3 methods for class "wisp") on a small synthetic model with more than one fitted parameter and a
# handful of bootstrap resamples -- see tests/testthat/helper-fixtures.R's make_small_multiparam_countdata(). This
# is decoupled from "does a specific RORB story reproduce" (that's the slow-tier real-data tests and the tutorial
# vignettes); it only checks that each function runs without error and returns the expected shape/structure, not
# specific values.
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

test_that("plot(model, type = \"ratecount\") runs on the synthetic fixture", {
  rc <- plot(model, type = "ratecount", species = "demoGene", print.plots = FALSE, verbose = FALSE)
  expect_true(is.list(rc))
  expect_true("plot_pred_context_demoContext_fixEff_demoGene" %in% names(rc))
  expect_s3_class(rc$plot_pred_context_demoContext_fixEff_demoGene, "ggplot")
})

test_that("plot(model, type = \"ratecount\", decompose = TRUE) runs on the synthetic fixture", {
  decomp <- plot(model, type = "ratecount", species = "demoGene", decompose = TRUE, print.plots = FALSE)
  expect_true(is.list(decomp))
  expect_gt(length(decomp), 0)
})

test_that("plot(model, type = \"ratecount\", decompose = TRUE) errors clearly without exactly one species", {
  expect_error(
    plot(model, type = "ratecount", decompose = TRUE, print.plots = FALSE, verbose = FALSE),
    "exactly one"
  )
  expect_error(
    plot(model, type = "ratecount", species = c("demoGene", "other"), decompose = TRUE, print.plots = FALSE, verbose = FALSE),
    "exactly one"
  )
})

test_that("decompose = TRUE is rejected for any type other than \"ratecount\"", {
  expect_error(
    plot(model, type = "timeseries", species = "demoGene", decompose = TRUE, print.plots = FALSE, verbose = FALSE),
    "only valid when type"
  )
})

test_that("plot(model, type = \"timeseries\") runs on the synthetic fixture", {
  ts <- plot(model, type = "timeseries", species = "demoGene", print.plots = FALSE, verbose = FALSE)
  expect_true(is.list(ts))
  expect_true("plot_pred_context_demoContext_timeseries_demoGene" %in% names(ts))
  expect_s3_class(ts$plot_pred_context_demoContext_timeseries_demoGene, "ggplot")
})

test_that("plot(model, type = \"parameters\") runs on the synthetic fixture", {
  pp <- plot(model, type = "parameters", species = "demoGene", print.plots = FALSE, verbose = FALSE)
  expect_true(is.list(pp))
  expect_true("plot_baseline_demoContext_demoGene" %in% names(pp))
  expect_true("plot_ranEff_demoContext_demoGene" %in% names(pp))
})

test_that("plot(model, type = \"MCMC.walks\") runs on the synthetic fixture", {
  mw <- plot(model, type = "MCMC.walks", print.plots = FALSE, verbose = FALSE, low_samples = 5)
  expect_true(is.list(mw))
  expect_true(all(c("plot.walks.parameters_low", "plot.walks.parameters_high", "plot.walks.nll") %in% names(mw)))
})

test_that("plot(model, type = \"MCMC.bs.comparison\") runs on the synthetic fixture (needs >1 parameter to not crash)", {
  bc <- plot(model, type = "MCMC.bs.comparison", print.plots = FALSE, verbose = FALSE)
  expect_true(is.list(bc))
  expect_true(all(c("plot_sample_correlations", "plot_comparison_Shaprio", "plot_comparison_density") %in% names(bc)))
  expect_s3_class(bc$plot_comparison_density, "ggplot")
})

test_that("species is ignored with a warning for types with no species argument", {
  expect_warning(
    plot(model, type = "effect.dist", species = "demoGene", print.plots = FALSE, verbose = FALSE),
    "species"
  )
})

test_that("sample.stats runs on the synthetic fixture", {
  ss <- sample.stats(model, verbose = FALSE)
  expect_true(is.data.frame(ss$stats$parameters))
  expect_gt(nrow(ss$stats$parameters), 1)
  expect_true(all(c("parameter", "estimate", "CI.low", "CI.high", "p.value") %in% names(ss$stats$parameters)))
})

test_that("print(model) runs on the synthetic fixture", {
  expect_output(print(model), "Fitted wisp model")
})

test_that("summary(model)/print(summary(model)) run on the synthetic fixture", {
  s <- summary(model)
  expect_s3_class(s, "summary.wisp")
  expect_true(is.data.frame(s$parameters))
  expect_output(print(s), "Summary of fitted wisp model")
})

test_that("residuals(model) runs on the synthetic fixture", {
  r_raw <- residuals(model)
  r_log <- residuals(model, type = "log")
  expect_true(is.data.frame(r_raw))
  expect_true(is.data.frame(r_log))
  expect_true(all(c("group", "mean", "sd", "variance") %in% names(r_raw)))
})

test_that("coef(model) runs on the synthetic fixture", {
  cf <- coef(model)
  expect_true(is.numeric(cf))
  expect_false(is.null(names(cf)))
  expect_equal(unname(cf), unname(model$fitted.parameters))
})

test_that("confint(model) runs on the synthetic fixture", {
  ci <- confint(model)
  expect_true(is.matrix(ci))
  expect_true(all(c("CI.low", "CI.high") %in% colnames(ci)))
  expect_equal(nrow(ci), nrow(model$stats$parameters))

  ci_subset <- confint(model, parm = model$param.names[1])
  expect_equal(nrow(ci_subset), 1)

  ci_90 <- confint(model, level = 0.9)
  expect_true(is.matrix(ci_90))
  expect_equal(nrow(ci_90), nrow(model$stats$parameters))
})
