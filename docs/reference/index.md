# Package index

## Main

- [`wisp()`](https://michaelbarkasi.github.io/wispack/reference/wisp.md)
  : Fit wisp to count data

## Model functions

- [`wisp.sigmoid()`](https://michaelbarkasi.github.io/wispack/reference/wisp.sigmoid.md)
  : Wisp sigmoid function, implemented in R
- [`wisp.warp()`](https://michaelbarkasi.github.io/wispack/reference/wisp.warp.md)
  : Wisp warping function, implemented in R

## Analysis functions

- [`analyze.residuals()`](https://michaelbarkasi.github.io/wispack/reference/analyze.residuals.md)
  : Analyze residuals from wisp fit
- [`pvalues.samples()`](https://michaelbarkasi.github.io/wispack/reference/pvalues.samples.md)
  : Compute p-values using ecdf from parameter resamples
- [`sample.stats()`](https://michaelbarkasi.github.io/wispack/reference/sample.stats.md)
  : Estimate p-values and confidence intervals from resampled parameters

## Plotting functions

- [`plot.MCMC.bs.comparison()`](https://michaelbarkasi.github.io/wispack/reference/plot-MCMC-bs-comparison.md)
  : Visually compare normality and autocorrelation of bootstrap and MCMC
  parameter estimates
- [`plot.MCMC.walks()`](https://michaelbarkasi.github.io/wispack/reference/plot-MCMC-walks.md)
  : Plot sampling of random walks and negative log likelihood from MCMC
  simulations
- [`plot.decomposition()`](https://michaelbarkasi.github.io/wispack/reference/plot-decomposition.md)
  : Plot individual components of wisp fit for a single species level
- [`plot.effect.dist()`](https://michaelbarkasi.github.io/wispack/reference/plot-effect-dist.md)
  : Plot parameter distributions from WISP results as histograms
- [`plot.parameters()`](https://michaelbarkasi.github.io/wispack/reference/plot-parameters.md)
  : Plot of wisp parameters
- [`plot.ratecount()`](https://michaelbarkasi.github.io/wispack/reference/plot-ratecount.md)
  : Plot fitted model and data
- [`plot.species.summary()`](https://michaelbarkasi.github.io/wispack/reference/plot-species-summary.md)
  : Print rate-count, residual, and parameter plots for one species
  level together.
- [`plot.timeseries()`](https://michaelbarkasi.github.io/wispack/reference/plot-timeseries.md)
  : Plot timeseries effects from wisp results
- [`demo.sigmoid.plots()`](https://michaelbarkasi.github.io/wispack/reference/demo.sigmoid.plots.md)
  : Make plot demonstrating how the wisp function is determined by the
  Rt, slope, and tpoint parameters
- [`demo.warp.plots()`](https://michaelbarkasi.github.io/wispack/reference/demo.warp.plots.md)
  : Make plot demonstrating how the wisp function is warped by the
  warping factors

## MERFISH preprocessing functions

- [`make_count_data()`](https://michaelbarkasi.github.io/wispack/reference/make_count_data.md)
  : Load and parse data, HDF5
- [`make_count_data_csv()`](https://michaelbarkasi.github.io/wispack/reference/make_count_data_csv.md)
  : Load and parse data, csv
- [`cortical_coordinate_transform()`](https://michaelbarkasi.github.io/wispack/reference/cortical_coordinate_transform.md)
  : Function to transform coordinates for each mouse and extract layer
  boundary estimates
- [`create.count.data.WSPmm()`](https://michaelbarkasi.github.io/wispack/reference/create.count.data.WSPmm.md)
  : Function to convert to WSPmm format
