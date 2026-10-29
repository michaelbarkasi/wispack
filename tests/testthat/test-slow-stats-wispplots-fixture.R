# Port of tutorial_stats.Rmd's checks 1-4, 6-7 and tutorial_wispplots.Rmd's checks 1-4, now that both
# vignettes fit the cortical-laminar RORB model (1000 bootstraps, max.fork = 20) on the fly instead of loading
# a shipped inst/extdata/corticallaminar_model.rds fixture (that file has been removed -- see
# R/golden_checks.R and the commit that dropped it). These checks used to live in the fast tier
# (test-wisp-stats-fixture.R / test-wisp-plots-fixture.R) by loading the pre-fit RDS; since there's no longer
# a pre-fit object to load for free, exercising them now means doing the real 1000-bootstrap fit, which only
# belongs in the slow tier.
#
# Fit once and reuse across both sets of checks (mirrors the vignettes' own call order, including the
# sample.stats()/plot(model, type = "MCMC.walks")/plot(model, type = "parameters") mutation sequence
# tutorial_stats.Rmd uses, and its median-model refit) rather than fitting twice, since a single real fit here
# is already the most expensive thing in this file. check5 (median_fit_param1) gets re-verified here as a side
# effect of replaying that sequence -- it's also covered independently by test-slow-stats-median.R, which only
# needs that one smaller fit and so stays separate rather than depending on this file's much larger bootstrap
# fit.
#
# tutorial_stats.Rmd's check8 (plot_comparison_density) is NOT ported here with an exact-value bound: it was
# originally, but that turned out to be irreproducible outside the literal knitr session that first produced
# it. plot(model, type = "MCMC.bs.comparison") selects which parameter to plot by ranking Shapiro-Wilk
# p-values computed from *random* resamples of each parameter's values (R/wspmm.R), so the result is sensitive
# to the exact R-level RNG stream position. Empirically, re-running this file's exact sequence (same seed, same
# call order, including the plot(model, type = "parameters")/median-refit calls that precede it) gave values
# differing by 4 orders of magnitude between a plain Rscript run and a testthat run -- not just "different
# knitr sessions disagree", but genuinely not portable across execution contexts. See R/golden_checks.R's note
# and tutorial_stats.Rmd's own check8 chunk (now a structural check, same treatment as
# tutorial_wispplots.Rmd's check4).

g  <- .wisp_golden$stats
gw <- .wisp_golden$wispplots

test_that("the real cortical-laminar bootstrap fit matches both tutorial_stats.Rmd's and tutorial_wispplots.Rmd's narratives", {
  skip_unless_slow_tests()
  skip_unless_extdata("S1_laminar_countdata_demo.csv.gz")
  skip_unless_extdata("S1_layer_boundary_bins.csv")

  layer.boundary.bins <- read.csv(system.file("extdata", "S1_layer_boundary_bins.csv", package = "wispack"))
  countdata <- read.csv(system.file("extdata", "S1_laminar_countdata_demo.csv.gz", package = "wispack"))
  data.variables <- list(species = "gene", ran = "mouse", timeseries = "age")

  set.seed(123)
  model <- wisp(
    count.data     = countdata,
    variables      = data.variables,
    model.settings = list(LRO_cost = "none"),
    plot.settings  = list(print.plots = FALSE, dim.bounds = colMeans(layer.boundary.bins)),
    bootstraps.num = 1000,
    max.fork       = 20,
    verbose        = FALSE
  )

  # stats check1: default stats on the freshly fit model
  expect_between(model$stats$parameters[4, 2], bounds = g$default_param4)

  # stats check2: re-run sample.stats with alpha = 0.01, Bonferroni correction
  model <- sample.stats(model, alpha = 0.01, Bonferroni = TRUE, verbose = FALSE)
  expect_between(model$stats$parameters[4, 3], bounds = g$bonferroni_param4)

  # stats check3: MCMC walk plot for a low-expression parameter
  MCMC_walk_plots <- plot(model, type = "MCMC.walks", print.plots = FALSE, verbose = FALSE, low_samples = 10)
  test_val <- ggplot2::ggplot_build(MCMC_walk_plots$plot.walks.parameters_low)$data[[1]][4, 2]
  expect_between(test_val, bounds = g$mcmc_walk_low4)

  # stats check4: re-derive stats from the MCMC resamples instead of bootstraps
  model_MCMC <- model
  model_MCMC$sample.params.type <- "MCMC"
  model_MCMC <- sample.stats(model_MCMC, verbose = FALSE)
  expect_between(model_MCMC$stats$parameters[4, 4], bounds = g$mcmc_resample_param4)

  # Not independently checked, but the vignette calls these between check4 and check5/6/7 -- replaying them
  # here keeps this file's call sequence matching the vignette's, for whatever it's worth given check8's
  # fragility is NOT fixed by this replay (see header note); kept mainly so plot(model, type = "parameters")
  # itself is still exercised against the real fit at least once in the suite.
  plots.param <- plot(model, type = "parameters", species = "Rorb", verbose = FALSE)
  invisible(plot(model_MCMC, type = "parameters", species = "Rorb", verbose = FALSE))
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

  # stats checks 6-7: MCMC-vs-bootstrap comparison plots, on the (Bonferroni-adjusted) bootstrap model.
  # (check8, plot_comparison_density, is checked structurally only -- see header note.)
  bs_comp_plots <- plot(model, type = "MCMC.bs.comparison", verbose = FALSE)

  test_val <- ggplot2::ggplot_build(bs_comp_plots$plot_sample_correlations)$data[[1]][4, 2]
  expect_between(test_val, bounds = g$bs_comparison_correlation)

  test_val <- ggplot2::ggplot_build(bs_comp_plots$plot_comparison_Shaprio)$data[[1]][4, 2]
  expect_between(test_val, bounds = g$bs_comparison_shapiro)

  d <- ggplot2::ggplot_build(bs_comp_plots$plot_comparison_density)$data[[1]]
  expect_true(all(c("x", "y", "colour") %in% names(d)))
  expect_true(is.finite(d[4, "y"]) && d[4, "y"] >= 0)

  # wispplots check1: plot(model, type = "ratecount")
  gb <- plot(model, type = "ratecount", dim.boundaries = colMeans(layer.boundary.bins), species = "Rorb", verbose = FALSE)
  test_val <- ggplot2::ggplot_build(gb$plot_pred_context_cortex_fixEff_Rorb)$data[[1]][4, 2]
  expect_between(test_val, bounds = gw$ratecount_Rorb)

  # wispplots check2: plot(model, type = "decomposition")
  plots.decomp <- plot(model, type = "decomposition", dim.boundaries = colMeans(layer.boundary.bins), species = "Rorb")
  layer_data <- ggplot2::ggplot_build(plots.decomp[[3]])$data[[1]]
  test_val <- layer_data[, c("x", "y")][43, 2]
  expect_between(test_val, bounds = gw$decomposition_panel3)

  # wispplots check3: plot(model, type = "timeseries")
  gb3 <- plot(model, type = "timeseries", species = "Rorb", verbose = FALSE)
  test_val <- gb3$plot_pred_context_cortex_timeseries_Rorb$data[2, 3]
  expect_between(test_val, bounds = gw$timeseries_Rorb)

  # wispplots check4: plot(model, type = "parameters")'s random-effect plot, checked structurally only -- see
  # tutorial_wispplots.Rmd's own check4 chunk and its comment for why. Reuses the plots.param computed
  # above (same plot(model, type = "parameters", species = "Rorb") call the vignette itself makes earlier).
  d <- plots.param$plot_ranEff_cortex_Rorb$data
  expect_true(all(c("means", "value", "parameter", "type", "species", "value_low", "value_high") %in% names(d)))
  expect_true(all(d$species == "Rorb"))
  expect_true(all(is.finite(d$value)))
  expect_true(all(is.finite(d$means)))
  expect_gt(nrow(d), 0)
})
