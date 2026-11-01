# Function to extract ground-truth from attractor-based simulation made with the function `attractor_simulation`

This function extracts the ground-truth values for functional spatial
effect values, replicate random effect values, and indicators for
spatially variable genes (SVGs) and genes with functional spatial
effects (FSEs) from a simulation object created by the
`attractor_simulation` function.

## Usage

``` r
attractor_simulation_ground_truth(sim)
```

## Arguments

- sim:

  A list object returned by the `attractor_simulation` function,
  containing the simulation data and parameters.

## Value

A list containing the following components: `fse_true`, `ran_true`,
`SVGs`, and `FSEs`.

## Examples

``` r
set.seed(1)
toy_seed <- data.frame(
  gene    = rep(c("geneA", "geneB"), each = 20),
  coord_x = runif(40, 0, 10),
  coord_y = runif(40, 0, 10),
  count   = rpois(40, 5)
)
sim <- attractor_simulation(toy_seed, n_bins = 5, n_replicates = 2)
truth <- attractor_simulation_ground_truth(sim)
truth$SVGs
#> geneA geneB 
#>  TRUE  TRUE 
```
