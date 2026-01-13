# Print rate-count, residual, and parameter plots for one child level together.

Function to summarize all important information for an individual child
level on one plot.

## Usage

``` r
# S3 method for class 'child.summary'
plot(wisp.results, these.parents = NULL, these.childs = NULL, verbose = TRUE)
```

## Arguments

- wisp.results:

  List, output of the wisp function.

- these.parents:

  Character vector, optional, specifies which parent levels to
  summarize. Defaults to all.

- these.childs:

  Character vector, optional, specifies which child levels to summarize.
  Defaults to all.

- verbose:

  Logical, if TRUE, prints information during plotting.

## Value

Nothing, plots are printed to the current graphics device.
