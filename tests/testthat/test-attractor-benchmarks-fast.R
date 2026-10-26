# Port of tutorial_benchmarks.Rmd's silent checks that don't require fitting a wisp model (checks 1-6, 8). These
# cover reading the shipped seed data and calling attractor_simulation() (fast: ~0.4s even on the real seed patch)
# plus aggregating a precomputed benchmark-results CSV with analyze_attractor_sim_benchmarks(). Check 7, which
# calls model_attractor_simulation_wisp() (100 bootstrap resamples), is the one part of this vignette that
# actually fits a model; it lives in test-slow-benchmarks.R.

g <- .wisp_golden$benchmarks

test_that("Allen MERFISH seed data matches the vignette narrative", {
  path <- skip_unless_extdata("Allen_data.csv")
  countdata <- read.csv(path)
  expect_between(countdata[1, 5], bounds = g$seed_data_val)
})

test_that("attractor_simulation on the seed patch matches the vignette narrative", {
  path <- skip_unless_extdata("Allen_data.csv")
  countdata <- read.csv(path)
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
    print_plots              = TRUE
  )

  test_val <- ggplot2::ggplot_build(test_sim$plots[["Vip"]][[1]])$data[[1]][4, 2]
  expect_between(test_val, bounds = g$Vip_smoothed)

  test_val <- ggplot2::ggplot_build(test_sim$plots[["Slc17a7"]][[1]])$data[[1]][7, 2]
  expect_between(test_val, bounds = g$Slc17a7_SVG)

  expect_between(test_sim$effect[3], bounds = g$Tac2_effect)

  test_val <- ggplot2::ggplot_build(test_sim$plots[["Tac2"]][[1]])$data[[1]][7, 2]
  expect_between(test_val, bounds = g$Tac2_FSE)

  expect_between(test_sim$replicate_rate_scalars[1], bounds = g$replicate_scalar1)
})

test_that("ground truth extracted from an attractor simulation matches its own stipulated values", {
  # Correctness test, independent of the regression bounds above: the ground-truth helper should report TRUE
  # exactly at the indices the simulation itself recorded as SVGs/FSEs.
  path <- skip_unless_extdata("Allen_data.csv")
  countdata <- read.csv(path)
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

  truth <- attractor_simulation_ground_truth(test_sim)

  expect_equal(names(truth$SVGs), test_sim$genes)
  expect_equal(sort(names(which(truth$SVGs))), sort(test_sim$SVGs))
  expect_equal(sort(names(which(truth$FSEs))), sort(test_sim$FSEs))
  expect_equal(truth$fse_true, test_sim$effect)
  expect_equal(truth$ran_true, test_sim$replicate_rate_scalars)
})

test_that("analyze_attractor_sim_benchmarks on the precomputed results matches the vignette narrative", {
  path <- skip_unless_extdata("benchmark_results_final250.csv")
  benchmark_results <- read.csv(path)
  benchmark_metrics <- analyze_attractor_sim_benchmarks(benchmark_results)

  expect_between(benchmark_metrics[1, 1], bounds = g$metrics_1_1)
  expect_between(benchmark_metrics[3, 2], bounds = g$metrics_3_2)
  expect_between(benchmark_metrics[8, 3], bounds = g$metrics_8_3)
  expect_between(benchmark_metrics[9, 1], bounds = g$metrics_9_1)
  expect_between(benchmark_metrics[9, 2], bounds = g$metrics_9_2)
  expect_between(benchmark_metrics[9, 3], bounds = g$metrics_9_3)
})
