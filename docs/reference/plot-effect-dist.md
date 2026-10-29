# Plot parameter distributions from WISP results as histograms

Function to make nicely formatted histograms of fitted parameters from
WISP results.

## Usage

``` r
plot.effect.dist(
 wisp.results, 
 print.plots = FALSE, 
 verbose     = TRUE
)
```

## Arguments

- wisp.results:

  List, output of the wisp function.

- print.plots:

  Logical, if TRUE, prints the plots to the current device.

- verbose:

  Logical, if TRUE, prints information during plotting.

## Value

List of ggplot objects containing histograms of fitted parameters.

## Note

Despite the dotted name, this is not an S3 method:
[`wisp()`](https://michaelbarkasi.github.io/wispack/reference/wisp.md)
results are plain, unclassed lists, so the generic
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) will never
dispatch here. Call `plot.effect.dist(...)` directly.
