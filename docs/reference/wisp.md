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
  all) named elements: `count`, `bin`, `context`, `species`, `ran`,
  `timeseries`, and `fixedeffects`.

- fit_only:

  Logical, if TRUE, only fits the model to the full data set using
  L-BFGS and returns the fitted model without running resampling (MCMC
  or bootstrapping); if FALSE, runs MCMC and/or bootstrapping to
  estimate parameter uncertainty.

- use.median:

  Logical, if TRUE, the median of the resamples is used as the final
  parameter estimates; if FALSE, the fit by L-BFGS on the full, original
  data is used.

- bootstraps.num:

  Integer, number of bootstrap resamples to perform. If 0, only MCMC is
  run.

- converged.resamples.only:

  Logical, if TRUE, only resamples with a converged fit are used for
  statistical analysis; if FALSE, all resamples are used. Applies only
  to bootstraps.

- max.fork:

  Integer, maximum number of parallel processes to use for
  bootstrapping. Applies only on Mac and Linux. Must be 1 for windows.

- verbose:

  Logical, if TRUE, prints information during the fitting process.

- model.settings:

  List, settings for the C++ model, including:

  - `buffer_factor` = 0.05: The minimum distance between t-points must
    be this proportion of `max_bin`.

  - `ctol` = 1e-6: Convergence tolerance for L-BFGS.

  - `max_penalty_at_distance_factor` = 0.01: Maximum penalty applied to
    negative log-likelihood when at a non-trivial distance from
    parameter bounds.

  - `LRO_cost` = "AIC": Character string giving cost function to use
    when evaluating likelihood-ratio outlier (LRO) parameters for
    testing for expression-rate transition points in the data ("AIC",
    "BIC", "NLL"; anything else defaults to the hand-specified LROcutoff
    and LROwindow_factor).

  - `LROcutoff` = 2.0: Cutoff used in the `LROcp_find` function, a
    multiple of standard deviation. Specifically, each point is
    evaluated as the ratio of the likelihood of the data on a model with
    a transition point to the likelihood of the data on a model without
    a change point. A point is taken to be a transition point only if
    this ratio is more than this number of standard deviations above the
    mean across all points.

  - `LROwindow_factor` = 1.25: controls size of window used in
    `LROcp_find` algorithm (window = LROwindow_factor X bin_num X
    buffer_factor)

  - `rise_threshold_factor` = 0.8: Amount of detected rise as fraction
    of total required to end run when making initial estimates of
    transition-point slopes.

  - `max_evals` = 1000: Maximum number of evaluations for optimization,
    when using L-BFGS.

  - `rng_seed` = 42: Seed for random number generator used internal to
    C++.

  - `warp_precision` = 1e-7: Decimal precision to retain when selecting
    really big number as pseudo infinity for unbound warping.

  - `round_none` = TRUE: Logical specifying whether to round
    extrapolated counts for "none" (no random effect) to nearest
    integer.

  - `trtKO` = c("none"): List of effect names to remove from the model,
    e.g., model without the Factor1 x Factor2 interaction. Must match
    one of the treatment levels of the full model.

  - `max_bin` = 0.0: Largest bin number; if left zero, will infer from
    count data

- MCMC.settings:

  List, settings for the MCMC simulation, including:

  - `MCMC.burnin` = 1e2: Number of random-walk steps to take to find the
    optimal parameters. If zero, will start the walk from the L-BFGS
    parameters.

  - `MCMC.steps` = 1e3: Number of random-walk steps to take when
    resampling parameters.

  - `MCMC.step.size` = 0.5: Standard deviation of the normal
    distribution controlling random-walk step size.

  - `MCMC.prior` = 1.0: Prior to compute the Bayesian update controlling
    step acceptance.

  - `MCMC.neighbor.filter` = 1: Will only save every this-many steps as
    a resample.

- plot.settings:

  List, settings for plots to make, including:

  - `print.plots` = TRUE: Logical indicating whether to print plots as
    they are produced during the model fit.

  - `dim.bounds` = Vector of bins; if non-empty, rate-count plots will
    plot vertical dot-dashed lines as these coordinates to indicate a
    prior-known boundary along the spatial dimension. Note that these
    values are not used in any way during the model fit itself.

  - `log.scale` = False: Logical indicating whether plots should show
    the log of expression-rate values (plus 1).

  - `splitting_factor` = c(): Name of a fixed-effect factor. Must be
    either empty or length 1. If provided, rate-count and time-series
    plots will determine the color hue of treatment levels by the value
    of this fixed-effect in that level. Within splitting factor levels,
    the other fixed-effect factor level determines color brightness.
    Note that this splitting-factor color scheme only has an effect if
    there are exactly two fixed-effect levels (one of which can be a
    time series, but the splitting factor itself cannot be a time
    series).

  - `CI_style` = TRUE: Logical indicating whether rate-count plots
    should show (TRUE:) predicted treatment-level expression rates free
    of individual effects plus shading showing 95% confidence invervals,
    or (FALSE:) observed data as points, plus fitted treatment-level
    expression rates for both individual samples and predictions free of
    individual effects.

  - `splitting_factor_colors` = c(120, 240): Vector of two integers
    (range 0:360) giving the hue of each level of the splitting factor.

  - `label_size`, `title_size`, `axis_size`, `legend_size`,
    `count_size`: Controls for plot element sizes.

  - `count_jitter` = 0.5: Scale of jitter to apply when plotting data
    for a given input.

  - `count.alpha.ran`, `count.alpha.none`, `pred.alpha.ran`, and
    `pred.alpha.none`: Controls for plot element alpha-levels.

- ran.seed:

  = 1234: Integer, random seed for reproducibility. Used within R. If
  NULL, no seed is set. Default is 1234.

## Value

List giving the results of the fitted model, including:

- `model.component.list` = "Rt", "tslope", "tpoint": Character vector
  giving the main model components; always "Rt" = rate, "tslope" =
  transition-point slope, and "tpoint" = transition-point location.

- `count.data.summed`: Data frame of the fitted data. There is one row
  per combination of spatial bin, fixed effects, and random effects.
  Contains columns for `bin`, `count` (sum of spots with a given species
  label), `pred` (predicted rate), `count.log` log of the count, plus
  one, `pred.log` (log of predicted rate plus one), `context`,
  `species`, `ran` (random-effect level), and `treatment` (treatment
  level).

- `fitted.parameters`: Vector of all fitted model parameters.

- `gamma.disperson`: Matrix giving the estimated gamma dispersion of the
  gamma-convolved Poisson model.

- `param.names`: Names of the model parameters, matching the order of
  `fitted.parameters`.

- `fix`: List containing four elements: `names` (names of fixed-effect
  factors), `lvls` (list, each element the names of the levels of a
  fixed-effect factor), `treat.lvl` (as before, but excluding reference
  levels), and `ref.lvl` (as before, but containing only the
  reference-level).

- `treatment`: List with two elements: `names` (names of the treatment
  levels) and `components` (list, each element the names of the
  fixed-effect factor levels constituting that treatment level).

- `grouping.variables`: List containing character vectors
  `context.lvls`, `species.lvls`, and `ran.lvls` holding, respectively,
  the grouping levels of the context, species, and random-effect
  grouping variabes.

- `param.idx0`: List containing zero-indexed index values for
  `fitted.parameters` for different combinations of treatment levels and
  grouping factors.

- `settings`: A copy of the original model settings provided to the
  `wisp` function.

- `sample.params.type`: A character string, "bs" or "MCMC", specifies
  which type of resamples were used to compute stats.

- `sample.params.bs`: A data frame giving the resampled parameters from
  bootstrapping (if any).

- `sample.params.MCMC`: A data frame giving the resampled parameters
  from the MCMC resampling (if any).

- `diagnostics.bs`: Diagnostic information for each bootstrap resample.

- `diagnostics.MCMC`: Diagnostic information for each MCMC random-walk
  step (i.e., resample).

- `stats`: List containing a main data frame `parameters` with the
  results (p-values, estimate, confidence intervals, etc) of the
  statistical analysis of the model parameters, plus data frames holding
  model residuals.

- `plots`: List containing whatever plots where generated during the
  model fit.
