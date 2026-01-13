# Plot parameter distributions from WISP results as histograms

Function to make nicely formatted histograms of fitted parameters from
WISP results.

## Usage

``` r
# S3 method for class 'effect.dist'
plot(wisp.results, print.plots = TRUE, verbose = TRUE)
```

## Arguments

- wisp.results:

  List, output of the wisp function.

- verbose:

  Logical, if TRUE, prints information during plotting.

## Value

List of ggplot objects containing histograms of fitted parameters.
