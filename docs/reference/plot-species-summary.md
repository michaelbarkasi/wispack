# Print rate-count, residual, and parameter plots for one species level together.

Function to summarize all important information for an individual
species level on one plot.

## Usage

``` r
plot.species.summary(
 wisp.results, 
 contexts = c(), 
 species  = c(), 
 verbose  = TRUE
)
```

## Arguments

- wisp.results:

  List, output of the wisp function.

- contexts:

  Character vector, optional, specifies which context levels to
  summarize. Defaults to all.

- species:

  Character vector, optional, specifies which species levels to
  summarize. Defaults to all.

- verbose:

  Logical, if TRUE, prints information during plotting.

## Value

Nothing, plots are printed to the current graphics device.

## Note

Despite the dotted name, this is not an S3 method:
[`wisp()`](https://michaelbarkasi.github.io/wispack/reference/wisp.md)
results are plain, unclassed lists, so the generic
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) will never
dispatch here. Call `plot.species.summary(...)` directly.

[`wisp()`](https://michaelbarkasi.github.io/wispack/reference/wisp.md)
does not store built plot objects on its return value, so this function
always rebuilds the rate-count, residual, and parameter plots it needs
fresh from `wisp.results` (via
[`plot.ratecount()`](https://michaelbarkasi.github.io/wispack/reference/plot-ratecount.md),
[`analyze.residuals()`](https://michaelbarkasi.github.io/wispack/reference/analyze.residuals.md),
and
[`plot.parameters()`](https://michaelbarkasi.github.io/wispack/reference/plot-parameters.md))
rather than reading them back from a stored `$plots` field.
