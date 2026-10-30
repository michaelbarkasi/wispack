# Structural tests for wisp()'s return value and for the package's exported API surface. These catch an
# accidental breaking change to the public contract (a renamed/missing list element, a new name silently added to
# NAMESPACE) independently of whether any fitted number drifts -- i.e. they're a different kind of regression
# guard than the ported vignette checks.

test_that("wisp() returns the documented top-level structure", {
  cd <- make_tiny_countdata(seed = 1)
  model <- wisp(
    cd,
    fit_only       = TRUE,
    model.settings = list(LRO_cost = "none"),
    plot.settings  = list(print.plots = FALSE),
    verbose        = FALSE
  )

  expect_type(model, "list")
  expect_true(all(c(
    "model.component.list", "count.data.summed", "fitted.parameters", "param.names",
    "fix", "treatment", "grouping.variables", "param.idx0", "settings"
  ) %in% names(model)))
  # wisp() does not store built plot objects on its return value (see R/wspmm.R); regenerate
  # plots on demand via plot(model, type = "ratecount"/"timeseries"/"parameters"/etc.) instead.
  expect_false("plots" %in% names(model))

  # wisp() classes its return value "wisp" so that plot(model, type = ...) dispatches correctly.
  expect_s3_class(model, "wisp")

  expect_equal(model$model.component.list, c("Rt", "tslope", "tpoint"))
  expect_type(model$fitted.parameters, "double")
  expect_length(model$param.names, length(model$fitted.parameters))
  expect_false(anyNA(model$fitted.parameters))

  expect_true(all(c("context.lvls", "species.lvls", "ran.lvls") %in% names(model$grouping.variables)))
  expect_true(all(c("names", "lvls", "treat.lvl", "ref.lvl") %in% names(model$fix)))
  expect_true(all(c("names", "components") %in% names(model$treatment)))
})

test_that("wisp() is deterministic given the same seed and settings", {
  cd <- make_tiny_countdata(seed = 1)
  args <- list(
    count.data     = cd,
    fit_only       = TRUE,
    model.settings = list(LRO_cost = "none"),
    plot.settings  = list(print.plots = FALSE),
    verbose        = FALSE,
    ran.seed       = 42
  )
  model1 <- do.call(wisp, args)
  model2 <- do.call(wisp, args)

  expect_equal(model1$fitted.parameters, model2$fitted.parameters)
  expect_equal(model1$param.names, model2$param.names)
})

test_that("the package's exported API surface hasn't grown or shrunk unintentionally", {
  # wispack's NAMESPACE now uses a curated export() list (one export() per function carrying a roxygen
  # @export tag in R/*.R) instead of a blanket exportPattern, so internal helpers and Rcpp-module/S4
  # artifacts no longer leak into the public surface just by being defined at top level. This test
  # freezes the current surface so that change is visible and intentional. When you add a genuinely new,
  # intended export, update the vector below in the same commit.
  current <- sort(getNamespaceExports("wispack"))
  expected <- sort(c(
    ".__C__C++Object", ".__C__Rcpp_wspc",
    # S4 method-table registration artifacts re-exported by import(ggplot2)/import(grid)/import(gridExtra)
    # (those packages setMethod() some base generics on their own S4 classes; importing them re-exports
    # the method-table entries too). Not something wispack defines itself.
    ".__T__[:base", ".__T__[[<-:base", ".__T__[<-:base", ".__T__$:base", ".__T__$<-:base",
    "analyze_attractor_sim_benchmarks", "analyze.residuals",
    "attractor_simulation", "attractor_simulation_ground_truth",
    "cortical_coordinate_transform", "create.count.data.WSPmm", "create.count.data.WSPmm.allen",
    "demo.sigmoid.plots", "demo.warp.plots",
    "make_count_data", "make_count_data_csv",
    "model_attractor_simulation_wisp",
    "run_attractor_sim_benchmarks", "sample.stats",
    "wisp", "wisp_sigmoid", "wisp_warp"
    # coef.wisp, confint.wisp, plot.wisp, print.wisp, print.summary.wisp, residuals.wisp, and summary.wisp
    # are S3 methods: they dispatch via the S3method() entries in NAMESPACE and don't need (and no longer
    # have) a separate export() entry for the dotted method name itself.
  ))
  expect_equal(current, expected)
})
