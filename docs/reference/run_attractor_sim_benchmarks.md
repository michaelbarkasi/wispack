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

## Examples

``` r
# A toy seed large enough per bin/replicate for the bootstrap resamples
# inside model_attractor_simulation_wisp() to converge reliably.
set.seed(1)
toy_seed <- data.frame(
  gene    = rep(c("geneA", "geneB"), each = 100),
  coord_x = runif(200, 0, 10),
  coord_y = runif(200, 0, 10),
  count   = rpois(200, 20)
)

# n_sims = 1 and a small bs_num keep this fast; a real benchmark would use
# many more simulations and bootstrap resamples.
bm <- run_attractor_sim_benchmarks(
  toy_seed,
  n_sims                 = 1,
  n_bins                 = 10,
  n_replicates           = 3,
  modeling_function_args = list(wisp = list(bs_num = 8, max_fork = 1))
)
#> 
#> 
#> ----------------------
#> Running simulation  1 / 1
#>  data sim time: 0.002 s
#> 
#> Modeling with  wisp ...
#>   wisp time: 4.265 s
bm$results
#>            est       true         param         id method sim
#> 1 -0.098954514  0.0000000   rate_effect      geneA   wisp   1
#> 2  0.007648286  0.0000000   rate_effect      geneB   wisp   1
#> 3 -0.257198208 -0.3297745 random_effect replicate1   wisp   1
#> 4  0.189563699 -0.2851778 random_effect replicate2   wisp   1
#> 5 -0.362414149 -0.3631447 random_effect replicate3   wisp   1
#> 6  0.401671078 -0.1792952 random_effect replicate4   wisp   1
#> 7  0.200000000  0.0000000           FSE      geneA   wisp   1
#> 8  0.172839506  0.0000000           FSE      geneB   wisp   1
```
