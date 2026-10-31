# Golden-value reference table for the tutorial_* vignettes ############################################################
#
# Every numeric bound checked (silently) inside the tutorial_*.Rmd vignettes' check chunks is recorded here, once,
# keyed by vignette and check label. Both the vignettes' own check chunks (via wispack:::.check_between(), below)
# and the tests/testthat slow-tier tests (tests/testthat/test-slow-*.R, via expect_between()) read these bounds
# instead of hard-coding them a second time, so there is exactly one place that encodes "the narrative numbers are
# X" -- rendering a vignette and running the testthat suite both verify the narrative against the same source.
#
# This object (and the helper below) are internal (not exported); they exist purely to support the vignettes and
# the test suite.

#' Assert that a value falls within a golden-value bound, vignette-check style.
#'
#' The vignette equivalent of tests/testthat/helper-expectations.R's expect_between(): same bounds table, same
#' inclusive-range check, but raising a plain error (via stop()) instead of a testthat expectation, since this is
#' called from knitr chunks, not from test_that() blocks.
#'
#' @param value Numeric value measured in the vignette.
#' @param bounds A `c(lo = ..., hi = ...)` entry from `.wisp_golden`.
#' @param label Short label identifying the check, used in the error message if it fails.
#' @noRd
.check_between <- function(value, bounds, label) {
  if (value < bounds[["lo"]] || value > bounds[["hi"]]) {
    stop(sprintf("failed %s: %s is not between %s and %s", label, format(value), format(bounds[["lo"]]), format(bounds[["hi"]])))
  }
  invisible(value)
}

#' @noRd
.wisp_golden <- list(

  poisson = list(
    countdata_bin1    = c(lo = 31,        hi = 31),
    countdata_count2  = c(lo = 8,         hi = 8),
    Rt1               = c(lo = 4.5,       hi = 4.75),
    tslope2_min       = 10,  # one-sided: tslope2 must be >= this
    tpoint2           = c(lo = 49,        hi = 51),
    recovered_rate4   = c(lo = 7.3,       hi = 7.5)
  ),

  warping = list(
    warp_ratio_demo        = c(lo = 0.2702110,  hi = 0.2702112),
    warped_sigmoid_layer5  = c(lo = 27.8595328, hi = 27.8595330)
  ),

  wispplots = list(
    ratecount_Rorb        = c(lo = 9.2,      hi = 9.3),
    decomposition_panel3  = c(lo = 16.3,     hi = 16.4),
    timeseries_Rorb       = c(lo = 11.11764, hi = 11.11766)
    # NOTE: the vignette's check4 (plot(model, type = "parameters")'s ranEff plot, data row 2) is deliberately
    # not reproduced here. plot(model, type = "parameters") appears to draw its displayed random-effect points
    # from the global RNG stream at call time rather than from a fixed stored sample, so that bound only ever
    # held because of knitr's chunk-level RNG-state caching across the vignette's specific rendering order --
    # it isn't a property of the fitted model.
    # See the corresponding (structural, not value-based) test in tests/testthat/test-slow-stats-wispplots-fixture.R.
  ),

  stats = list(
    default_param4             = c(lo = 0.04,    hi = 0.05),
    bonferroni_param4          = c(lo = -0.1,    hi = -0.09),
    mcmc_walk_low4             = c(lo = -0.93,   hi = -0.9),
    mcmc_resample_param4       = c(lo = 0.05,    hi = 0.07),
    median_fit_param1          = c(lo = 3.4,     hi = 3.5),
    bs_comparison_correlation  = c(lo = 0.003,   hi = 0.004),
    bs_comparison_shapiro      = c(lo = 0.8,     hi = 0.9)
    # NOTE: the vignette's check8 (plot(model, type = "MCMC.bs.comparison")'s plot_comparison_density, data
    # row 4) is deliberately not reproduced here, for the same reason as the wispplots group's check4 below:
    # plot(model, type = "MCMC.bs.comparison") ranks parameters by a Shapiro-Wilk p-value computed from a *random* resample
    # of each parameter's bootstrap/MCMC values (see R/wspmm.R), so which column ends up selected for display
    # depends on the exact R-level RNG stream position at call time -- confirmed empirically to differ by
    # orders of magnitude between otherwise-identical re-runs (plain script vs. testthat vs. knitr all land on
    # different values), not just between knitr sessions. Checked structurally instead; see
    # tests/testthat/test-slow-stats-wispplots-fixture.R and tutorial_stats.Rmd's own check8 chunk.
  ),

  multicontext = list(
    n_celltypes_raw      = c(lo = 31, hi = 31),
    n_celltypes_filtered = c(lo = 2,  hi = 2),
    beta_Rt_GABA_Bcl11b_right18_Blk2 = c(lo = -0.17,     hi = -0.16),
    ratecount_Glut                   = c(lo = 2.890371,  hi = 2.890373),
    ratecount_Glut_Bcl11b            = c(lo = 3.1,       hi = 3.2),
    ratecount_Glut_Cux2              = c(lo = 2.3,       hi = 2.4),
    ratecount_Glut_Rorb              = c(lo = 1.9,       hi = 2.0),
    ratecount_Glut_Nxph3             = c(lo = 5.2,       hi = 5.3),
    timeseries_Glut_Bcl11b           = c(lo = 1.37,      hi = 1.4),
    timeseries_Glut_Cux2             = c(lo = 2.4,       hi = 2.5),
    timeseries_Glut_Fezf2            = c(lo = 1.2,       hi = 1.3),
    timeseries_Glut_Nxph3            = c(lo = 1.6,       hi = 1.7),
    timeseries_Glut_Rorb             = c(lo = 1.9,       hi = 2.0)
  ),

  timeseries = list(
    ratecount_Nxph3 = c(lo = 219, hi = 220)
  ),

  corticallaminar = list(
    n_genes          = c(lo = 6,     hi = 6),
    second_gene_name = "Fezf2",
    n_cells          = c(lo = 54116, hi = 54116),
    param1_estimate  = c(lo = 3.38,  hi = 3.4),
    hemisphere_RORB  = c(lo = -0.01, hi = 0.1)
  ),

  liverradial = list(
    param4_estimate   = c(lo = 0.13, hi = 0.16),
    zonation_space_glul   = c(lo = 4.4, hi = 4.5),
    zonation_time_glul    = c(lo = 6.9, hi = 7.0),
    rhythm_space_arntl    = c(lo = 2.4, hi = 2.6),
    rhythm_time_arntl     = c(lo = 1.9, hi = 2.0),
    both_space_elovl3     = c(lo = 5.4, hi = 5.5),
    both_time_elovl3      = c(lo = 4.6, hi = 4.8)
  ),

  benchmarks = list(
    seed_data_val          = c(lo = 2.180366,   hi = 2.180368),
    Vip_smoothed           = c(lo = 2.116861,   hi = 2.116863),
    Slc17a7_SVG            = c(lo = 2.063262,   hi = 2.063264),
    Tac2_effect            = c(lo = -0.4785649, hi = -0.4785647),
    Tac2_FSE               = c(lo = 3.691496,   hi = 3.691498),
    replicate_scalar1      = c(lo = -0.3201774, hi = -0.3201772),
    wisp_results_7_2       = c(lo = 0.101,      hi = 0.102),
    wisp_results_12_1      = c(lo = 1.4,        hi = 1.5),
    metrics_1_1            = c(lo = 0.227,      hi = 0.229),
    metrics_3_2            = c(lo = 0.017,      hi = 0.019),
    metrics_8_3            = c(lo = 0.427,      hi = 0.429),
    metrics_9_1            = c(lo = 1,          hi = 3.000),
    metrics_9_2            = c(lo = 26,         hi = 28.00),
    metrics_9_3            = c(lo = 426,        hi = 428.0)
  )

)
