# Plot timeseries effects from wisp results

Function to plot timeseries effects from wisp results, separated by a
splitting factor (e.g., hemisphere).

## Usage

``` r
plot.timeseries(wisp.results, splitting_factor = NULL, pred.type = "pred")
```

## Arguments

- wisp.results:

  List, output of the wisp function.

- splitting_factor:

  Character string, the name of the fixed effect by which to split the
  timeseries. If NULL, the first non-timeseries fixed effect is used. If
  "none", will not split.

- pred.type:

  Character string, the name of the predicted rate column in the count
  data (e.g., "pred.log" or "pred").

- count.type:

  Character string, the name of the observed count column in the count
  data (e.g., "count.log" or "count").

- print.all:

  Logical, if TRUE, prints all plots; if FALSE, only returns plots in
  list without printing any.

- verbose:

  Logical, if TRUE, prints updates about the plotting process.

## Value

List of ggplot objects for timeseries plots.
