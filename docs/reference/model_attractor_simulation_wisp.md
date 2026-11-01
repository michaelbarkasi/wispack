# Function to model simulation with wisp

This function fits a wisp model to an attractor-based simulation dataset
and extracts the estimated functional spatial effects, replicate random
effects, and p-values for functional spatial effects. It then compiles
these results along with the ground-truth values from the simulation for
comparison.

## Usage

``` r
model_attractor_simulation_wisp(
  sim,
  sim_num,
  LRO_cost = "none",
  bs_num = 1000,
  max_fork = 1
)
```

## Arguments

- sim:

  A list object returned by the `attractor_simulation` function,
  containing the simulation data and parameters.

- sim_num:

  An integer specifying the simulation number for setting the random
  seed and labeling results.

- LRO_cost:

  Passed to the `wisp` function. Character string giving cost function
  to use when evaluating likelihood-ratio outlier (LRO) parameters for
  testing for expression-rate transition points in the data ("AIC",
  "BIC", "NLL"; anything else defaults to the hand-specified LROcutoff
  and LROwindow_factor).

- bs_num:

  An integer specifying the number of bootstrap resamples to run.

- max_fork:

  An integer specifying the number of forks to run when computing
  bootstraps in parallel. For Windows, must be 1.

## Value

A data frame containing the estimated values (`est`), ground-truth
values (`true`), parameter types (`param`), parameter-associated gene or
replicate (`id`), modeling method used (`method`), and simulation number
(`sim`).

## Examples

``` r
# A larger toy seed than attractor_simulation()'s own example: fitting a wisp
# model needs enough data per bin/replicate for the bootstrap resamples to
# converge reliably.
set.seed(1)
toy_seed <- data.frame(
  gene    = rep(c("geneA", "geneB"), each = 100),
  coord_x = runif(200, 0, 10),
  coord_y = runif(200, 0, 10),
  count   = rpois(200, 20)
)
sim <- attractor_simulation(toy_seed, n_bins = 10, n_replicates = 3)

# Small bs_num keeps this fast; a real benchmark would use many more resamples.
results <- model_attractor_simulation_wisp(sim, sim_num = 1, bs_num = 8, max_fork = 1)
results
#>           est       true         param         id method sim
#> 1  0.17075617  0.0000000   rate_effect      geneA   wisp   1
#> 2  0.02472505  0.0000000   rate_effect      geneB   wisp   1
#> 3 -0.09512236 -0.3297745 random_effect replicate1   wisp   1
#> 4  0.23173355 -0.2851778 random_effect replicate2   wisp   1
#> 5 -0.19984913 -0.3631447 random_effect replicate3   wisp   1
#> 6  0.48148148  0.0000000           FSE      geneA   wisp   1
#> 7  1.05555556  0.0000000           FSE      geneB   wisp   1
```
