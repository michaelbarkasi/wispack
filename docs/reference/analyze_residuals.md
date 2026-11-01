# Analyze residuals from wisp fit

This function takes a vector `mu.B` of sampled values of a variable X,
and a single observation `mu.obs` of the same variable, and computes the
p-value of `mu.obs` using the empirical cumulative distribution function
(ecdf) of `mu.B`.

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

A list containing data frames `stats` and `stats.log` holding summary
stats for residuals and log-residuals, respectively, plus a list of
plots visualizing these stats.

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
