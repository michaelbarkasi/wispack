# Function to run benchmarks on attractor-based simulations

This function benchmarks a provided list of modeling methods on
attractor-based simulation datasets. It simulates data multiple times,
fits each modeling method to the simulated data, and records the results
and computation times for each step.

## Usage

``` r
run_attractor_sim_benchmarks(
  seed_data,
  n_sims = 100,
  n_bins = 100,
  n_replicates = 4,
  replicate_spatial_scalar = 0.05,
  min_effect_size = 0.05,
  modeling_functions = list(wisp = model_attractor_simulation_wisp),
  modeling_function_args = list(wisp = list(bs_num = 1000, max_fork = 1))
)
```

## Arguments

- seed_data:

  A data frame containing the seed dataset with columns: `gene`,
  `coord_x`, `coord_y`, and `count`. Will be given to
  `attractor_simulation` to generate simulations.

- n_sims:

  An integer specifying the number of simulations to run. Default is
  100.

- n_bins:

  An integer specifying the number of spatial bins to divide the data
  coordinates into for each simulation. Default is 100.

- n_replicates:

  An integer specifying the number of replicates to generate for each
  treatment condition in each simulation. Default is 4.

- replicate_spatial_scalar:

  A numeric value controlling the amount of spatial variation introduced
  in each replicate for each simulation. Default is 0.05.

- min_effect_size:

  A numeric value specifying the minimum effect size for functional
  spatial effects (FSEs) in each simulation. Default is 0.05. Max
  positive effect size is always 4, i.e., a 4x change in rate. There is
  no max to the min effect size, i.e., an effect can drop the rate to
  zero.

- modeling_functions:

  A named list of modeling functions to benchmark. Each function should
  take at least the arguments `sim` (the simulation object produced by
  `attractor_simulation`) and `sim_num` (the simulation number). Default
  is a list containing only the `model_attractor_simulation_wisp`
  function. Functions provided in this list should return a dataframe
  with columns `est`, `true`, `param`, `id`, `method`, and `sim`. The
  `true` column should contain, for each simulation and each applicable
  parameter, the ground-truth value returned by the
  `attractor_simulation_ground_truth` function. The `est` column should
  contain the corresponding estimated value from the modeling function.
  The `param` column should contain the name of the parameter (one of
  "rate_effect", "random_effect", "FSE", or "SVG"). The `id` column
  should contain the name of the gene or replicate associated with the
  parameter. The `method` column should contain the name of the modeling
  method used. The `sim` column should contain the simulation number.

- modeling_function_args:

  A named list of lists, where each sub-list contains additional
  arguments to pass to the corresponding modeling function in
  `modeling_functions`. Default is a list containing arguments for the
  `model_attractor_simulation_wisp` function.

## Value

A list containing two components: `results`, a data frame compiling the
results from all simulations and modeling methods; and `times`, a data
frame recording the computation times for data simulation and each
modeling method for each simulation.
