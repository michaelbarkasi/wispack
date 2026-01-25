# Function to simulate data using attractor-based spatial variability

This function takes a seed dataset and simulates spatially variable
genes (SVGs) with functional spatial effects (FSEs) based on attractor
points. It generates multiple replicates for reference and treatment
conditions, applying random spatial transformations and rate scalars.

## Usage

``` r
attractor_simulation(
  seed_data,
  n_bins = 100,
  n_replicates = 4,
  replicate_spatial_scalar = 0.05,
  min_effect_size = 0.05,
  print_plots = FALSE
)
```

## Arguments

- seed_data:

  A data frame containing the seed dataset with at least the columns:
  `gene`, `coord_x`, `coord_y`, and `count`.

- n_bins:

  An integer specifying the number of spatial bins to divide the data
  coordinates into. Default is 100.

- n_replicates:

  An integer specifying the number of replicates to generate for each
  treatment condition. Default is 4.

- replicate_spatial_scalar:

  A numeric value controlling the amount of spatial variation introduced
  in each replicate for each simulation. Default is 0.05.

- min_effect_size:

  A numeric value specifying the minimum effect size for functional
  spatial effects (FSEs) in each simulation. Default is 0.05. Max
  positive effect size is always 4, i.e., a 4x change in rate. There is
  no max to the min effect size, i.e., an effect can drop the rate to
  zero.

- print_plots:

  A logical value indicating whether to return plots of the simulation
  process for each gene. Default is FALSE.

## Value

A list containing the following components: `genes`, `SVGs_idx`, `SVGs`,
`FSEs_idx`, `FSEs`, `attractor`, `effect`, `replicate_rate_scalars`,
`data`, and `plots`.
