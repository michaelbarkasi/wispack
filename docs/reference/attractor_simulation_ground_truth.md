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
