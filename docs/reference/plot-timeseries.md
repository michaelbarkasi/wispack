# Plot timeseries effects from wisp results

Function to plot timeseries effects from wisp results, separated by a
splitting factor (e.g., hemisphere).

## Usage

``` r
plot.timeseries(
 wisp.results, 
 splitting_factor = c(), 
 tertiary_effects = NULL,
 log.scale        = FALSE,
 print.plots      = FALSE,
 species          = c(),
 verbose          = TRUE
)
```

## Arguments

- wisp.results:

  List, output of the wisp function.

- splitting_factor:

  Character string, the name of the fixed effect by which to split the
  timeseries. If empty, c(), the first non-timeseries fixed effect is
  used. If "none", will not split.

- tertiary_effects:

  Character vector, or NULL (default). Any fixed effects besides
  `timeseries` and the splitting factor are pruned down to a single
  level before rates are reconstructed (otherwise treatments belonging
  to different levels of those extra fixed effects would be mixed
  together). If NULL, each such fixed effect is pruned down to its
  reference level. Otherwise, for each additional fixed effect, if one
  of the levels named in `tertiary_effects` matches one of its levels,
  that level is used instead of the reference level.

- log.scale:

  Logical, if TRUE, plots rates and counts on a log scale; if FALSE,
  plots raw rates and counts.

- print.plots:

  Logical, if TRUE, prints plots; if FALSE, only returns plots in list
  without printing any.

- species:

  Character vector, list of species levels to plot. If c(), all species
  levels are plotted.

- verbose:

  Logical, if TRUE, prints updates about the plotting process.

## Value

List of ggplot objects for timeseries plots.

## Note

Despite the dotted name, this is not an S3 method:
[`wisp()`](https://michaelbarkasi.github.io/wispack/reference/wisp.md)
results are plain, unclassed lists, so the generic
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) will never
dispatch here. Call `plot.timeseries(...)` directly.
