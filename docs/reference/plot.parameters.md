# Plot of wisp parameters

Function to make nicely formatted violin (or bar) plots of the fitted
wisp parameters, including confidence intervals if stat information is
available.

## Usage

``` r
# S3 method for class 'parameters'
plot(
  wisp.results,
  species.lvl = NULL,
  violin = TRUE,
  print.plots = FALSE,
  species.classes = NULL,
  verbose = TRUE
)
```

## Arguments

- wisp.results:

  List, output of the wisp function.

- species.lvl:

  Character string, the species level to be plotted. If NULL, all
  species levels are plotted.

- violin:

  Logical, if TRUE, plots violin plots for each parameter; if FALSE,
  uses bar plots.

- print.plots:

  Logical, if TRUE, prints the plots to the console; if FALSE, only
  returns a list of plots without printing.

- species.classes:

  List, a list of character vectors specifying which species levels to
  include together in plots. If NULL, all species levels are included in
  a single plot.

- verbose:

  Logical, if TRUE, prints updates about the plotting process.

## Value

List of ggplot objects for parameter plots.
