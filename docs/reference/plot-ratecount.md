# Plot fitted model and data

This function takes wisp model results (a fitted line) and observed data
(counts) and plots them together for visual comparison. It can also
include independent block boundaries for comparison if provided.

## Usage

``` r
plot.ratecount(
 wisp.results, 
 log.scale = FALSE,
 dim.boundaries = c(), 
 print.plots = FALSE, 
 y.lim = NA, 
 count.alpha.none = NA, 
 count.alpha.ran = NA, 
 pred.alpha.none = NA, 
 pred.alpha.ran = NA, 
 rans.to.print = c(), 
 species = c(), 
 verbose = TRUE
)
```

## Arguments

- wisp.results:

  List, output of the wisp function.

- log.scale:

  Logical, if TRUE, plots rates and counts on a log scale; if FALSE,
  plots raw rates and counts.

- CI_style:

  Logical, if TRUE, plots predicted lines with confidence interval style
  shading; if FALSE, plots predicted lines as solid lines with random
  levels and data points

- dim.boundaries:

  Numeric vector, independent block boundaries to plot for comparison.
  If empty, the argument is ignored.

- print.plots:

  Logical, if TRUE, prints plots; if FALSE, only returns plots in list
  without printing any.

- y.lim:

  Numeric vector of length 2, limits for the y-axis of the plots. If NA,
  defaults to automatic limits.

- count.alpha.none:

  Numeric, transparency for count points when random level is "none". If
  left NA, defaults to 0.25.

- count.alpha.ran:

  Numeric, transparency for count points when random level is not
  "none". If left NA, defaults to 0.25.

- pred.alpha.none:

  Numeric, transparency for predicted lines when random level is "none".
  If left NA, defaults to 1.0.

- pred.alpha.ran:

  Numeric, transparency for predicted lines when random level is not
  "none". If left NA, defaults to 0.9.

- rans.to.print:

  Character vector, list of random levels to include on each species
  plot. If c(), all random levels are included.

- species:

  Character vector, list of species levels to plot. If c(), all species
  levels are plotted.

- verbose:

  Logical, if TRUE, prints updates about the plotting process.

## Value

List of ggplot objects for rate-count plots.
