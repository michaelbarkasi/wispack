# Analyze residuals from wisp fit

This function takes a vector `mu.B` of sampled values of a variable X,
and a single observation `mu.obs` of the same variable, and computes the
p-value of `mu.obs` using the empirical cumulative distribution function
(ecdf) of `mu.B`.

## Usage

``` r
analyze.residuals(wisp.results, verbose = TRUE)
```

## Arguments

- wisp.results:

  List, output of the wisp function.

- verbose:

  Logical, if TRUE, prints information during the statistical analysis.

## Value

Numeric p-value.
