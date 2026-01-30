# Print rate-count, residual, and parameter plots for one species level together.

Function to summarize all important information for an individual
species level on one plot.

## Usage

``` r
plot.species.summary(
 wisp.results, 
 contexts = c(), 
 species = c(), 
 verbose = TRUE
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
