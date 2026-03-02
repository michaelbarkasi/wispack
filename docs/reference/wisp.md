# Fit wisp to count data

This function takes a data frame of wisp variables (as columns) and fits
a wisp model to it. Statistical analyses and plots are generated from
the fitted model.

## Usage

``` r
wisp(
  count.data,
  variables = list(),
  fit_only = FALSE,
  use.median = FALSE,
  bootstraps.num = 0,
  converged.resamples.only = TRUE,
  max.fork = 1,
  verbose = TRUE,
  model.settings = list(),
  MCMC.settings = list(),
  plot.settings = list(),
  ran.seed = 1234
)
```

## Arguments

- count.data:

  Data.frame, data to be modeled, with columns for model variables
  (count, bin, context, species, ran, fixedeffects), or equivalent
  variables as specified in the `variables` argument.

- variables:

  List, names of the columns in `count.data` that correspond to the
  model variables. The list should contain only (but not necessarily
  all) named elements: `count`, `bin`, `context`, `species`, `ran`, and
  `fixedeffects`.

- fit_only:

  Logical, if TRUE, only fits the model to the full data set using
  L-BFGS and returns the fitted model without running MCMC or
  bootstrapping; if FALSE, runs MCMC and/or bootstrapping to estimate
  parameter uncertainty.

- use.median:

  Logical, if TRUE, the median of the resamples is used as the final
  parameter estimates; if FALSE, the initial fit by L-BFGS is used.

- bootstraps.num:

  Integer, number of bootstrap resamples to perform. If 0, only MCMC is
  run.

- converged.resamples.only:

  Logical, if TRUE, only resamples with a converged fit are used for
  statistical analysis; if FALSE, all resamples are used. Applies only
  to bootstraps.

- max.fork:

  Integer, maximum number of parallel processes to use for
  bootstrapping.

- verbose:

  Logical, if TRUE, prints information during the fitting process.

- model.settings:

  List, settings for the C++ model, including `buffer_factor`, `ctol`,
  `max_penalty_at_distance_factor`, `LROcutoff`, `LROwindow_factor`,
  `rise_threshold_factor`, `max_evals`, `rng_seed`, `warp_precision`,
  `round_none`, and `trtKO`. Default values are provided.

- MCMC.settings:

  List, settings for the MCMC simulation, including `MCMC.burnin`,
  `MCMC.steps`, `MCMC.step.size`, `MCMC.prior`, and
  `MCMC.neighbor.filter`. Default values are provided.

- plot.settings:

  List, settings for plots to make, including `print.plots`,
  `dim.bounds`, `log.scale`, `splitting_factor`, `CI_style`,
  `label_size`, `title_size`, `axis_size`, `legend_size`, `count_size`,
  `count_jitter`, `count.alpha.ran`, `count.alpha.none`,
  `pred.alpha.ran`, and `pred.alpha.none`. Default values are provided.

- ran.seed:

  Integer, random seed for reproducibility. If NULL, no seed is set.
  Default is 1234.

## Value

List giving the results of the fitted model, including:
`model.component.list`, `count.data.summed`, `fitted.parameters`,
`gamma.disperson`, `param.names`, `fix`, `treatment`,
`grouping.variables`, `param.idx0`, `settings`, `sample.params`,
`sample.params.bs`, `sample.params.MCMC`, `diagnostics.bs`,
`diagnostics.MCMC`, `stats`, and `plots`
