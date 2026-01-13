# Visually compare normality and autocorrelation of bootstrap and MCMC parameter estimates

Function allowing for visual comparison of the parameter estimates from
bootstrapping vs MCMC simulation. Shows density for ten representative
parameters from bootstrapping and MCMC walk, the distributions of
Shapiro-Walk p-values for bootstrapping vs MCMC walks, and the density
of autocorrelation results for bootstrapping vs MCMC walks.

## Usage

``` r
# S3 method for class 'MCMC.bs.comparison'
plot(wisp.results, print.plots = TRUE, verbose = TRUE)
```

## Arguments

- wisp.results:

  List, output of the wisp function.

- verbose:

  Logical, if TRUE, prints information during plotting.

## Value

List of ggplot objects containing plots of representative parameter
distributions, Shapiro-Wilk normality test results, and autocorrelation
plots for bootstrap and MCMC parameter estimates.
