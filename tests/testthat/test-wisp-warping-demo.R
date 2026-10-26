# Port of tutorial_warping.Rmd's silent checks. No wisp() model fit is involved (demo.warp.plots just evaluates
# the warp/sigmoid primitives and builds ggplot objects), so this runs in the fast tier.

g <- .wisp_golden$warping

test_that("the warp-ratio demo function matches the vignette narrative", {
  set.seed(123)
  warp_ratio <- function(X, rho, b) 1 - exp(rho^2)^(-X / b)

  rhos <- c(0, 0.75, 1.5, 3)
  bs <- c(100, 1e3)
  max_X <- 100

  demo_data <- data.frame(
    X     = rep(c(1:max_X), length(bs) * length(rhos)),
    rho   = factor(rep(rep(rhos, each = max_X), length(bs))),
    b     = factor(rep(bs, each = max_X * length(rhos))),
    ratio = c(
      vapply(bs,
        function(j) c(
          vapply(
            rhos,
            function(r) warp_ratio(c(1:max_X), r, j),
            numeric(max_X)
          )
        ),
        numeric(max_X * length(rhos))
      )
    )
  )

  expect_between(demo_data[214, 4], bounds = g$warp_ratio_demo)
})

test_that("demo.warp.plots matches the vignette narrative", {
  set.seed(123)
  warp_plots <- demo.warp.plots(
    w            = 2,
    point_pos    = 60,
    point_neg    = 40,
    Rt           = c(6, 3, 0.2, 6) * 4.65,
    tslope       = c(0.4, 0.75, 1),
    tpoint       = c(15, 38, 80),
    w_factors    = c(0.6, -0.9, 0.5),
    return_plots = TRUE
  )

  expect_s3_class(warp_plots$demo_plot_warpfunction, "ggplot")
  expect_s3_class(warp_plots$demo_plot_warpedsigmoid, "ggplot")

  gb <- ggplot2::ggplot_build(warp_plots[["demo_plot_warpedsigmoid"]])
  layer_data <- gb$data[[1]]
  expect_between(layer_data[, c("x", "y")][5, 2], bounds = g$warped_sigmoid_layer5)
})
