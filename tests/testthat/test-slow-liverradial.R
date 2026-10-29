# Port of tutorial_liverradial.Rmd's silent checks 1-7. All seven read off the same fitted model (1000
# bootstraps, max.fork = 10, MCMC disabled) -- this is the single most expensive test in the suite (minutes, not
# seconds), because unlike the other real-data fits ported here, this one actually bootstraps. It is only run
# when WISPACK_RUN_SLOW_TESTS is set.

g <- .wisp_golden$liverradial

test_that("wisp() on the Droin liver radial-zonation simulation matches the vignette narrative", {
  skip_unless_slow_tests()
  skip_unless_extdata("Droin_radial_count_data_sim.csv.gz")

  count_data <- read.csv(system.file("extdata", "Droin_radial_count_data_sim.csv.gz", package = "wispack"))
  data.variables <- list(
    count      = "count",
    bin        = "bin",
    context    = "liver",
    species    = "gene",
    ran        = "mouse",
    timeseries = "ZT"
  )
  model.settings <- list(LRO_cost = "none", LROcutoff = 1.5)
  plot.settings <- list(
    print.plots      = FALSE,
    title_size       = 14,
    log.scale        = TRUE,
    count_size       = 2.5,
    count_jitter     = 0.0,
    count.alpha.none = 0.8,
    count.alpha.ran  = 0.5,
    pred.alpha.ran   = 0.0
  )

  set.seed(123)
  radial.model <- wisp(
    count.data     = count_data,
    variables      = data.variables,
    bootstraps.num = 1e3,
    max.fork       = 10,
    model.settings = model.settings,
    MCMC.settings  = list(MCMC.steps = 0),
    plot.settings  = plot.settings
  )

  # check1
  expect_between(radial.model$stats$parameters[4, 2], bounds = g$param4_estimate)

  rc_plots <- plot(
    radial.model,
    type             = "ratecount",
    log.scale        = radial.model$plot.settings$log.scale,
    CI_style         = radial.model$plot.settings$CI_style,
    dim.boundaries   = radial.model$plot.settings$dim.bounds,
    count.alpha.none = radial.model$plot.settings$count.alpha.none,
    count.alpha.ran  = radial.model$plot.settings$count.alpha.ran,
    pred.alpha.none  = radial.model$plot.settings$pred.alpha.none,
    pred.alpha.ran   = radial.model$plot.settings$pred.alpha.ran,
    print.plots      = FALSE,
    verbose          = FALSE
  )
  ts_plots <- plot(
    radial.model,
    type             = "timeseries",
    splitting_factor = radial.model$plot.settings$splitting_factor,
    log.scale        = radial.model$plot.settings$log.scale,
    print.plots      = FALSE,
    verbose          = FALSE
  )

  # check2
  plt_glul <- rc_plots[["plot_pred.log_context_liver_fixEff_glul"]]
  test_val <- ggplot2::ggplot_build(plt_glul)$data[[1]][4, 2]
  expect_between(test_val, bounds = g$zonation_space_glul)

  # check3
  plt_glul <- ts_plots[["plot_pred.log_context_liver_timeseries_glul"]]
  test_val <- ggplot2::ggplot_build(plt_glul)$data[[1]][17, 2]
  expect_between(test_val, bounds = g$zonation_time_glul)

  # check4
  plt_arntl <- rc_plots[["plot_pred.log_context_liver_fixEff_arntl"]]
  test_val <- ggplot2::ggplot_build(plt_arntl)$data[[1]][4, 2]
  expect_between(test_val, bounds = g$rhythm_space_arntl)

  # check5
  plt_arntl <- ts_plots[["plot_pred.log_context_liver_timeseries_arntl"]]
  test_val <- ggplot2::ggplot_build(plt_arntl)$data[[1]][4, 2]
  expect_between(test_val, bounds = g$rhythm_time_arntl)

  # check6
  plt_elovl3 <- rc_plots[["plot_pred.log_context_liver_fixEff_elovl3"]]
  test_val <- ggplot2::ggplot_build(plt_elovl3)$data[[1]][4, 2]
  expect_between(test_val, bounds = g$both_space_elovl3)

  # check7
  plt_elovl3 <- ts_plots[["plot_pred.log_context_liver_timeseries_elovl3"]]
  test_val <- ggplot2::ggplot_build(plt_elovl3)$data[[1]][4, 2]
  expect_between(test_val, bounds = g$both_time_elovl3)
})
