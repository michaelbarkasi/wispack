# Compute p-values using ecdf from parameter resamples

This function takes a vector `mu.B` of sampled values of a variable X,
and a single observation `mu.obs` of the same variable, and computes the
p-value of `mu.obs` using the empirical cumulative distribution function
(ecdf) of `mu.B`.

## Usage

``` r
pvalues.samples(mu.B, mu.obs)
```

## Arguments

- mu.B:

  Numeric vector of sampled values of a variable X, e.g., bootstrapped
  or MCMC estimates, used to make an empirical cumulative distribution
  function (ecdf).

- mu.obs:

  Numeric value of the observed variable X, e.g., the mean of `mu.B` or
  an actual observation.

## Value

Numeric p-value.
