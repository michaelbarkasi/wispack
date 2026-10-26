# Port of tutorial_benchmarks.Rmd's check 7, the one check in that vignette which fits actual wisp models
# (model_attractor_simulation_wisp with 100 bootstrap resamples, max_fork = 10; ~15s). Checks 1-6 and 8
# (data loading, attractor_simulation(), and aggregating a precomputed benchmark CSV) don't fit any model and
# live in test-attractor-benchmarks-fast.R (fast tier).

g <- .wisp_golden$benchmarks

test_that("model_attractor_simulation_wisp on the Allen seed patch matches the vignette narrative", {
  skip_unless_slow_tests()
  skip_unless_extdata("Allen_data.csv")

  countdata <- read.csv(system.file("extdata", "Allen_data.csv", package = "wispack"))
  seed_patch <- countdata[
    countdata$slice_num == 33 &
      countdata$coord_y >= 2 & countdata$coord_y <= 4 &
      countdata$coord_x >= 1 & countdata$coord_x <= 3, ]

  set.seed(12399)
  test_sim <- attractor_simulation(
    seed_data                = seed_patch,
    n_bins                   = 100,
    n_replicates             = 4,
    replicate_spatial_scalar = 0.05,
    min_effect_size          = 0.05,
    print_plots              = FALSE
  )

  wisp_results <- model_attractor_simulation_wisp(
    sim      = test_sim,
    sim_num  = 1,
    bs_num   = 100,
    max_fork = 10
  )

  expect_between(wisp_results[7, 2], bounds = g$wisp_results_7_2)
  expect_between(wisp_results[12, 1], bounds = g$wisp_results_12_1)
})
