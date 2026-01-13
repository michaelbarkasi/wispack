# Plot sampling of random walks and negative log likelihood from MCMC simulations

Function to make nicely formatted histograms of fitted parameters from
WISP results.

## Usage

``` r
plot.MCMC.walks(wisp.results, print.plots = TRUE, verbose = TRUE, low_samples = 10)
```

## Arguments

- wisp.results:

  List, output of the wisp function.

- print.plots:

  Logical, if TRUE, prints the plots to the current device.

- verbose:

  Logical, if TRUE, prints information during plotting.

- low_samples:

  Integer, number of low-value parameters to plot. Defaults to 10.

## Value

List of ggplot objects containing plots of walks for low-value
parameters, high-value parameters, and normalized negative log
likelihood.
