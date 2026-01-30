# Visually compare normality and autocorrelation of bootstrap and MCMC parameter estimates

Function allowing for visual comparison of the parameter estimates from
bootstrapping vs MCMC simulation. Shows density for ten representative
parameters from bootstrapping and MCMC walk, the distributions of
Shapiro-Walk p-values for bootstrapping vs MCMC walks, and the density
of autocorrelation results for bootstrapping vs MCMC walks.

## Usage

``` r
plot.MCMC.bs.comparison(
 wisp.results, 
 print.plots = FALSE, 
 verbose = TRUE
)
```

## Arguments

- wisp.results:

  List, output of the wisp function.

- print.plots:

  Logical, if TRUE, prints the plots to the current graphics device.

- verbose:

  Logical, if TRUE, prints information during plotting.

- ran.seed:

  Integer, random seed for reproducibility. Defaults to 1234. If NULL,
  no seed set.

## Value

List of ggplot objects containing plots of representative parameter
distributions, Shapiro-Wilk normality test results, and autocorrelation
plots for bootstrap and MCMC parameter estimates.
