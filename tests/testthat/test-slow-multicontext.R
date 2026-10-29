# Port of tutorial_multicontext.Rmd's silent checks 3-13 (the ones that need a real wisp(fit_only = TRUE) fit,
# ~10s, on real extdata). Checks 1-2 (data-only) are in test-multicontext-data.R (fast tier). All of checks 3-13
# read off the *same* fitted model, exactly as the vignette does, so this is written as one test with the fit
# done once.

g <- .wisp_golden$multicontext

test_that("wisp() on the multicontext (celltype) data matches the vignette narrative", {
  skip_unless_slow_tests()
  skip_unless_extdata("ACx_laminar_countdata_demo_with_celltype.csv.gz")
  skip_unless_extdata("ACx_layer_boundary_bins.csv")

  countdata <- read.csv(system.file("extdata", "ACx_laminar_countdata_demo_with_celltype.csv.gz", package = "wispack"))
  Glut_mask <- grepl("Glut", countdata$celltype_MMC)
  GABA_mask <- grepl("GABA", countdata$celltype_MMC)
  countdata$celltype_MMC[Glut_mask] <- "Glut"
  countdata$celltype_MMC[GABA_mask] <- "GABA"
  countdata <- countdata[Glut_mask | GABA_mask, ]

  layer.boundary.bins <- read.csv(system.file("extdata", "ACx_layer_boundary_bins.csv", package = "wispack"))
  data.variables <- list(context = "celltype_MMC", species = "gene", ran = "mouse", timeseries = "age")

  model <- wisp(
    count.data     = countdata,
    variables      = data.variables,
    fit_only       = TRUE,
    model.settings = list(LRO_cost = "none", LROwindow_factor = 1.75),
    plot.settings  = list(
      print.plots     = FALSE,
      dim.bounds      = colMeans(layer.boundary.bins),
      title_size      = 12,
      log.scale       = TRUE,
      count.alpha.ran = 0.0,
      pred.alpha.ran  = 0.0
    ),
    verbose = FALSE
  )

  # check3
  test_val <- model[["fitted.parameters"]][["beta_Rt_GABA_Bcl11b_right18_X_Tns/Blk2"]]
  expect_between(test_val, bounds = g$beta_Rt_GABA_Bcl11b_right18_Blk2)

  rc_plots <- plot(
    model,
    type             = "ratecount",
    log.scale        = model$plot.settings$log.scale,
    CI_style         = model$plot.settings$CI_style,
    dim.boundaries   = model$plot.settings$dim.bounds,
    count.alpha.none = model$plot.settings$count.alpha.none,
    count.alpha.ran  = model$plot.settings$count.alpha.ran,
    pred.alpha.none  = model$plot.settings$pred.alpha.none,
    pred.alpha.ran   = model$plot.settings$pred.alpha.ran,
    print.plots      = FALSE,
    verbose          = FALSE
  )

  # check4
  test_val <- ggplot2::ggplot_build(rc_plots[["plot_pred.log_context_Glut"]])$data[[1]][4, 2]
  expect_between(test_val, bounds = g$ratecount_Glut)

  # check5
  test_val <- ggplot2::ggplot_build(rc_plots[["plot_pred.log_context_Glut_fixEff_Bcl11b"]])$data[[1]][4, 2]
  expect_between(test_val, bounds = g$ratecount_Glut_Bcl11b)

  # check6
  test_val <- ggplot2::ggplot_build(rc_plots[["plot_pred.log_context_Glut_fixEff_Cux2"]])$data[[1]][4, 2]
  expect_between(test_val, bounds = g$ratecount_Glut_Cux2)

  # check7
  test_val <- ggplot2::ggplot_build(rc_plots[["plot_pred.log_context_Glut_fixEff_Rorb"]])$data[[1]][4, 2]
  expect_between(test_val, bounds = g$ratecount_Glut_Rorb)

  # check8
  test_val <- ggplot2::ggplot_build(rc_plots[["plot_pred.log_context_Glut_fixEff_Nxph3"]])$data[[1]][4, 2]
  expect_between(test_val, bounds = g$ratecount_Glut_Nxph3)

  ts_plots <- plot(
    model,
    type             = "timeseries",
    splitting_factor = model$plot.settings$splitting_factor,
    log.scale        = model$plot.settings$log.scale,
    print.plots      = FALSE,
    verbose          = FALSE
  )

  # check9
  test_val <- ggplot2::ggplot_build(ts_plots[["plot_pred.log_context_Glut_timeseries_Bcl11b"]])$data[[1]][4, 2]
  expect_between(test_val, bounds = g$timeseries_Glut_Bcl11b)

  # check10
  test_val <- ggplot2::ggplot_build(ts_plots[["plot_pred.log_context_Glut_timeseries_Cux2"]])$data[[1]][4, 2]
  expect_between(test_val, bounds = g$timeseries_Glut_Cux2)

  # check11
  test_val <- ggplot2::ggplot_build(ts_plots[["plot_pred.log_context_Glut_timeseries_Fezf2"]])$data[[1]][4, 2]
  expect_between(test_val, bounds = g$timeseries_Glut_Fezf2)

  # check12
  test_val <- ggplot2::ggplot_build(ts_plots[["plot_pred.log_context_Glut_timeseries_Nxph3"]])$data[[1]][4, 2]
  expect_between(test_val, bounds = g$timeseries_Glut_Nxph3)

  # check13
  test_val <- ggplot2::ggplot_build(ts_plots[["plot_pred.log_context_Glut_timeseries_Rorb"]])$data[[1]][4, 2]
  expect_between(test_val, bounds = g$timeseries_Glut_Rorb)
})
