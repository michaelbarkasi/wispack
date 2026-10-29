# Plot individual components of wisp fit for a single species level

Extension of `plot.ratecount` which plots the individual pieces of the
rate-count plot for a single species separately. Returns individual
plots for data points only, fit lines only, and data points plus fit
lines for individual random levels.

## Usage

``` r
plot.decomposition(
 wisp.results,
 species,
 log.scale      = FALSE,
 dim.boundaries = c(),
 y.lim          = NULL
)
```

## Arguments

- wisp.results:

  List, output of the wisp function.

- species:

  Character string, the species level to plot. Must be provided, and
  only one at a time.

- log.scale:

  Logical, if TRUE, plots on a log scale. Defaults to FALSE.

- dim.boundaries:

  Numeric vector, independent block boundaries to plot for comparison.
  If empty, the argument is ignored.

- y.lim:

  Numeric vector of length 2, limits for the y-axis of the plots. If NA,
  defaults to automatic limits.

## Value

List of ggplot objects containing the decomposed rate-count plots.

## Note

Despite the dotted name, this is not an S3 method:
[`wisp()`](https://michaelbarkasi.github.io/wispack/reference/wisp.md)
results are plain, unclassed lists, so the generic
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) will never
dispatch here. Call `plot.decomposition(...)` directly.
