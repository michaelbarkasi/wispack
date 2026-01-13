# Print rate-count, residual, and parameter plots for one species level together.

Function to summarize all important information for an individual
species level on one plot.

## Usage

``` r
plot.species.summary(wisp.results, these.contexts = NULL, these.speciess = NULL, verbose = TRUE)
```

## Arguments

- wisp.results:

  List, output of the wisp function.

- these.contexts:

  Character vector, optional, specifies which context levels to
  summarize. Defaults to all.

- these.speciess:

  Character vector, optional, specifies which species levels to
  summarize. Defaults to all.

- verbose:

  Logical, if TRUE, prints information during plotting.

## Value

Nothing, plots are printed to the current graphics device.
