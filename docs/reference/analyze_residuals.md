# Analyze residuals from wisp fit

Computes raw and log-scale residuals (predicted minus observed count,
per row of `wisp.results$count.data.summed`) for a fitted
[`wisp`](https://michaelbarkasi.github.io/wispack/reference/wisp.md)
model, then summarizes and visualizes them by group – overall, by
random-effect level, by treatment level, and by context/species
combination. For each group, builds a histogram of the log-residuals and
a Q-Q plot comparing their ordered values against both a standard-normal
and a gamma-convolved-normal theoretical quantile (using the model's
fitted gamma dispersion), plus the mean, standard deviation, and
variance of the residuals.

## Usage

``` r
analyze_residuals(wisp.results, verbose = TRUE)
```

## Arguments

- wisp.results:

  List, output of the wisp function.

- verbose:

  Logical, if TRUE, prints information during the statistical analysis.

## Value

A list with three elements: `stats` and `stats.log`, data frames (one
row per group) giving the `group` name plus `mean`, `sd`, and `variance`
of the raw and log residuals, respectively; and `plots`, a flat list of
the histogram, Q-Q, and side-by-side combined plots for each group,
named `"<group>_hist"`, `"<group>_qq"`, and `"<group>_resid"`.

## Examples

``` r
set.seed(1)
n_bins <- 20
bin <- sample(seq_len(n_bins), 400, replace = TRUE)
lambda <- rep(c(2, 8), each = n_bins / 2)
countdata <- data.frame(bin = bin, count = rpois(400, lambda[bin]))

# analyze_residuals() only needs predicted vs. observed counts, so it works
# even on a fit_only = TRUE model (no resampling required).
model <- wisp(
  countdata,
  fit_only       = TRUE,
  model.settings = list(LRO_cost = "none"),
  plot.settings  = list(print.plots = FALSE),
  verbose        = FALSE
)
resid_summary <- analyze_residuals(model, verbose = FALSE)
resid_summary$stats
#>             group     mean       sd variance
#> 1             all 5.503603 22.31537 497.9756
#> 2     ran_lvl_ran 5.503603 22.31537 497.9756
#> 3   treatment_ref 5.503603 22.31537 497.9756
#> 4 context_species 5.503603 22.31537 497.9756
```
