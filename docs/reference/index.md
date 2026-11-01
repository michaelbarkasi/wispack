# Package index

## Fitting & using a wisp model

- [`wisp()`](https://michaelbarkasi.github.io/wispack/reference/wisp.md)
  : Fit wisp to count data
- [`print(`*`<wisp>`*`)`](https://michaelbarkasi.github.io/wispack/reference/print.wisp.md)
  : Print a fitted wisp model
- [`summary(`*`<wisp>`*`)`](https://michaelbarkasi.github.io/wispack/reference/summary.wisp.md)
  : Summarize a fitted wisp model
- [`print(`*`<summary.wisp>`*`)`](https://michaelbarkasi.github.io/wispack/reference/print.summary.wisp.md)
  : Print a wisp model summary
- [`coef(`*`<wisp>`*`)`](https://michaelbarkasi.github.io/wispack/reference/coef.wisp.md)
  : Fitted parameter estimates from a wisp model
- [`confint(`*`<wisp>`*`)`](https://michaelbarkasi.github.io/wispack/reference/confint.wisp.md)
  : Confidence intervals for wisp model parameters
- [`residuals(`*`<wisp>`*`)`](https://michaelbarkasi.github.io/wispack/reference/residuals.wisp.md)
  : Residuals from a fitted wisp model
- [`plot(`*`<wisp>`*`)`](https://michaelbarkasi.github.io/wispack/reference/plot.wisp.md)
  : Plot a fitted wisp model

## Core model functions

- [`wisp_warp()`](https://michaelbarkasi.github.io/wispack/reference/wisp_warp.md)
  : Wisp warping function, implemented in R
- [`demo_warp_plots()`](https://michaelbarkasi.github.io/wispack/reference/demo_warp_plots.md)
  : Make plot demonstrating how the wisp function is warped by the
  warping factors
- [`wisp_sigmoid()`](https://michaelbarkasi.github.io/wispack/reference/wisp_sigmoid.md)
  : Wisp sigmoid function, implemented in R
- [`demo_sigmoid_plots()`](https://michaelbarkasi.github.io/wispack/reference/demo_sigmoid_plots.md)
  : Make plot demonstrating how the wisp function is determined by the
  Rt, slope, and tpoint parameters

## Analysis functions

- [`analyze_residuals()`](https://michaelbarkasi.github.io/wispack/reference/analyze_residuals.md)
  : Analyze residuals from wisp fit
- [`sample_stats()`](https://michaelbarkasi.github.io/wispack/reference/sample_stats.md)
  : Estimate p-values and confidence intervals from resampled parameters

## Benchmarking functions

- [`attractor_simulation_ground_truth()`](https://michaelbarkasi.github.io/wispack/reference/attractor_simulation_ground_truth.md)
  :

  Function to extract ground-truth from attractor-based simulation made
  with the function `attractor_simulation`

- [`attractor_simulation()`](https://michaelbarkasi.github.io/wispack/reference/attractor_simulation.md)
  : Function to simulate data using attractor-based spatial variability

- [`model_attractor_simulation_wisp()`](https://michaelbarkasi.github.io/wispack/reference/model_attractor_simulation_wisp.md)
  : Function to model simulation with wisp

- [`run_attractor_sim_benchmarks()`](https://michaelbarkasi.github.io/wispack/reference/run_attractor_sim_benchmarks.md)
  : Function to run benchmarks on attractor-based simulations

- [`analyze_attractor_sim_benchmarks()`](https://michaelbarkasi.github.io/wispack/reference/analyze_attractor_sim_benchmarks.md)
  : Function to analyze benchmark results from attractor-based
  simulations

## MERFISH preprocessing functions

- [`make_count_data()`](https://michaelbarkasi.github.io/wispack/reference/make_count_data.md)
  : Load and parse data, HDF5
- [`make_count_data_csv()`](https://michaelbarkasi.github.io/wispack/reference/make_count_data_csv.md)
  : Load and parse data, csv
- [`cortical_coordinate_transform()`](https://michaelbarkasi.github.io/wispack/reference/cortical_coordinate_transform.md)
  : Function to transform coordinates for each mouse and extract layer
  boundary estimates
- [`create_count_data_WSPmm()`](https://michaelbarkasi.github.io/wispack/reference/create_count_data_WSPmm.md)
  : Function to convert to WSPmm format
- [`create_count_data_WSPmm_allen()`](https://michaelbarkasi.github.io/wispack/reference/create_count_data_WSPmm_allen.md)
  : Function to convert Allen-formatted MERFISH data to WSPmm format
