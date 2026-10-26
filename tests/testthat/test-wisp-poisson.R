# Port of tutorial_Poisson.Rmd's silent checks. This is a full-scale replication of the vignette's own simulation
# (same seed, same sample size) -- it just happens to also be fast (fit_only = TRUE on a 1-D, 100-bin problem fits
# in well under a second), so it runs in the default/fast tier rather than being gated behind
# WISPACK_RUN_SLOW_TESTS. See R/golden_checks.R for the shared bounds and tests/testthat/test-slow-*.R for the
# checks that *do* need gating because of real extdata or real MCMC/bootstrap costs.

g <- .wisp_golden$poisson

test_that("simulated Poisson-process count data matches the vignette narrative", {
  set.seed(123)
  n_bins <- 100
  x <- seq(1, n_bins)
  n_instances <- 10000
  bin <- sample(x, n_instances, replace = TRUE)
  lambda <- rep(c(1, 5, 2, 7), each = length(x) / 4)
  var_gamma <- 2
  process_lambda <- rgamma(
    n_instances,
    shape = lambda[bin]^2 / var_gamma,
    rate  = lambda[bin]  / var_gamma
  )
  count <- rpois(n_instances, process_lambda)
  countdata <- data.frame(bin, count)

  expect_between(countdata[1, 1], bounds = g$countdata_bin1)
  expect_between(countdata[2, 2], bounds = g$countdata_count2)

  model <- wisp(
    countdata,
    fit_only       = TRUE,
    model.settings = list(LRO_cost = "none"),
    plot.settings  = list(print.plots = FALSE),
    verbose        = FALSE
  )

  param_names <- names(model$fitted.parameters)
  est_param <- round(model$fitted.parameters, 3)

  expect_between(est_param[grepl("Rt", param_names)][1], bounds = g$Rt1)
  expect_at_least(est_param[grepl("tslope", param_names)][2], g$tslope2_min)
  expect_between(est_param[grepl("tpoint", param_names)][2], bounds = g$tpoint2)

  # Recover the original block-4 rate, as the vignette does, and check it matches the simulated value of 7.
  n_bins_local <- 100
  tokens_per_bin <- rep(NA, n_bins_local)
  for (b in seq_len(n_bins_local)) tokens_per_bin[b] <- sum(countdata$bin == b)
  block <- 4
  start_bin <- (block - 1) * 25 + 1
  end_bin   <- block * 25
  mean_tokens <- mean(tokens_per_bin[start_bin:end_bin])
  rate_param  <- model$fitted.parameters[block]
  recovered_rate <- (exp(rate_param) - 1) / mean_tokens

  expect_between(recovered_rate, bounds = g$recovered_rate4)
})
