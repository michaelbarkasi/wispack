# Plot timeseries effects from wisp results

Function to plot timeseries effects from wisp results, separated by a
splitting factor (e.g., hemisphere).

## Usage

``` r
plot.timeseries(
 wisp.results, 
 splitting_factor = NULL, 
 log.scale = FALSE,
 print.plots = FALSE,
 species = c(),
 verbose = TRUE
)
```

## Arguments

- wisp.results:

  List, output of the wisp function.

- splitting_factor:

  Character string, the name of the fixed effect by which to split the
  timeseries. If NULL, the first non-timeseries fixed effect is used. If
  "none", will not split.

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
