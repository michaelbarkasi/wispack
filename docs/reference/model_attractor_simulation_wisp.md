# Function to model simulation with wisp

This function fits a wisp model to an attractor-based simulation dataset
and extracts the estimated functional spatial effects, replicate random
effects, and p-values for functional spatial effects. It then compiles
these results along with the ground-truth values from the simulation for
comparison.

## Usage

``` r
model_attractor_simulation_wisp(sim, sim_num, bs_num = 1000, max_fork = 1)
```

## Arguments

- sim:

  A list object returned by the `attractor_simulation` function,
  containing the simulation data and parameters.

- sim_num:

  An integer specifying the simulation number for setting the random
  seed and labeling results.

## Value

A data frame containing the estimated values (`est`), ground-truth
values (`true`), parameter types (`param`), parameter-associated gene or
replicate (`id`), modeling method used (`method`), and simulation number
(`sim`).
