# Port of tutorial_stats.Rmd's silent checks that operate on the pre-fit fixture model (checks 1-4, 6-8). Only
# check 5 (a real wisp() fit with use.median = TRUE on real extdata) needs an actual model fit; that one lives in
# test-slow-stats-median.R instead. Everything here just loads inst/extdata/corticallaminar_model.rds and calls
# sample.stats()/plot.*() on it, so it's fast.
#
# The checks mutate `model` in sequence, exactly as the vignette does (sample.stats(Bonferroni = TRUE) is applied
# and never undone before the later checks), so this is written as one test with the checks in vignette order
# rather than as independent tests that would each need to redo the earlier mutations.

g <- .wisp_golden$stats

test_that("sample.stats / MCMC-walk / bootstrap-comparison plots on the fixture model match the vignette narrative", {
  path <- skip_unless_extdata("corticallaminar_model.rds")
  set.seed(123)
  model <- readRDS(path)

  # check1: default stats on the loaded model
  expect_between(model$stats$parameters[4, 2], bounds = g$default_param4)

  # check2: re-run sample.stats with alpha = 0.01, Bonferroni correction
  model <- sample.stats(
    model,
    alpha      = 0.01,
    Bonferroni = TRUE,
    verbose    = FALSE
  )
  expect_between(model$stats$parameters[4, 3], bounds = g$bonferroni_param4)

  # check3: MCMC walk plot for a low-expression parameter
  MCMC_walk_plots <- plot.MCMC.walks(
    model,
    print.plots = FALSE,
    verbose     = FALSE,
    low_samples = 10
  )
  test_val <- ggplot2::ggplot_build(MCMC_walk_plots$plot.walks.parameters_low)$data[[1]][4, 2]
  expect_between(test_val, bounds = g$mcmc_walk_low4)

  # check4: re-derive stats from the MCMC resamples instead of bootstraps
  model_MCMC <- model
  model_MCMC$sample.params.type <- "MCMC"
  model_MCMC <- sample.stats(model_MCMC, verbose = FALSE)
  expect_between(model_MCMC$stats$parameters[4, 4], bounds = g$mcmc_resample_param4)

  # checks 6-8: MCMC-vs-bootstrap comparison plots, back on the (Bonferroni-adjusted) bootstrap model
  bs_comp_plots <- plot.MCMC.bs.comparison(model, verbose = FALSE)

  test_val <- ggplot2::ggplot_build(bs_comp_plots$plot_sample_correlations)$data[[1]][4, 2]
  expect_between(test_val, bounds = g$bs_comparison_correlation)

  test_val <- ggplot2::ggplot_build(bs_comp_plots$plot_comparison_Shaprio)$data[[1]][4, 2]
  expect_between(test_val, bounds = g$bs_comparison_shapiro)

  test_val <- ggplot2::ggplot_build(bs_comp_plots$plot_comparison_density)$data[[1]][4, 2]
  expect_between(test_val, bounds = g$bs_comparison_density)
})
