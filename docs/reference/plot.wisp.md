# Plot a fitted wisp model

S3 method for plotting the results returned by
[`wisp`](https://michaelbarkasi.github.io/wispack/reference/wisp.md).
This is the single public entry point for all wisp plotting: pick which
plot to make with `type`, and pass that plot's own extra arguments
through `...`. (This replaces the pre-2.6.0 API, where each plot had its
own dotted, `plot.*`-named function that looked like an S3 method but
never was, since
[`wisp()`](https://michaelbarkasi.github.io/wispack/reference/wisp.md)
results used to be plain, unclassed lists; see
`vignette("release_notes", package = "wispack")` for the old-name to
`type` mapping.)

## Usage

``` r
# S3 method for class 'wisp'
plot(
  x,
  type,
  print.plots = TRUE,
  verbose = TRUE,
  species = "all",
  decompose = FALSE,
  ...
)
```

## Arguments

- x:

  List, output of the
  [`wisp`](https://michaelbarkasi.github.io/wispack/reference/wisp.md)
  function (an object of class `"wisp"`).

- type:

  Character string, which plot to make. One of:

  - `"ratecount"`: plot the fitted model and observed data (counts)
    together for visual comparison, optionally with independent block
    boundaries. Set `decompose = TRUE` to instead plot the individual
    pieces of this plot, for a single species, separately. See the
    "ratecount" section below for its extra arguments.

  - `"timeseries"`: plot timeseries effects, separated by a splitting
    factor (e.g., hemisphere). See the "timeseries" section below.

  - `"parameters"`: violin (or bar) plots of the fitted wisp parameters,
    including confidence intervals if available. See the "parameters"
    section below.

  - `"effect.dist"`: histograms of the fitted parameters. See the
    "effect.dist" section below.

  - `"species.summary"`: rate-count, residual, and parameter plots for
    one species level, combined on one plot. See the "species.summary"
    section below.

  - `"MCMC.walks"`: sampling of random walks and negative log likelihood
    from MCMC simulations. See the "MCMC.walks" section below.

  - `"MCMC.bs.comparison"`: visually compare normality and
    autocorrelation of bootstrap and MCMC parameter estimates. See the
    "MCMC.bs.comparison" section below.

- print.plots:

  Logical, defaults to TRUE. If TRUE, prints the plot(s) to the current
  graphics device as they are built; if FALSE, only returns them
  (without printing). Forwarded to every `type` except
  `"species.summary"`, which has no print-control of its own and always
  prints its combined plot (see the "species.summary" section below).

- verbose:

  Logical, defaults to TRUE. If TRUE, prints updates about the plotting
  process. Forwarded to every `type`.

- species:

  Character vector, or the string `"all"` (the default), giving which
  species level(s) to include. `"all"` is translated internally to
  [`c()`](https://rdrr.io/r/base/c.html), this package's own "no filter,
  include everything" convention. Only forwarded to
  `type %in% c("ratecount", "timeseries", "parameters", "species.summary")`,
  since only those underlying plotting functions have a `species`
  argument; passing a non-default `species` for any other `type` raises
  a warning and is otherwise ignored. When `decompose = TRUE`, `species`
  must resolve to exactly one species level (i.e. not `"all"`, and not a
  multi-element vector), or a clear error is raised.

- decompose:

  Logical, defaults to FALSE. Only meaningful when `type = "ratecount"`;
  if TRUE, plots the individual pieces of the rate-count plot (data
  points only, fit lines only, and data points plus fit lines for
  individual random levels) for a single species separately, instead of
  the ordinary combined rate-count plot. This replaces the pre-2.6.0
  `type = "decomposition"` (see
  `vignette("release_notes", package = "wispack")`). Setting
  `decompose = TRUE` for any other `type` is an error.

- ...:

  Additional arguments passed on to the plotting function selected by
  `type`; see the sections below for each type's own arguments.

## Value

The return value depends on `type`; see the sections below. Most types
return a list of ggplot objects, except `type = "species.summary"`,
which returns nothing (plots are printed to the current graphics device)
and `type = "MCMC.walks"`, which returns `NULL` early if the model has
no MCMC sampling.

## ratecount

Takes wisp model results (a fitted line) and observed data (counts) and
plots them together for visual comparison. Can also include independent
block boundaries for comparison if provided. With `decompose = TRUE`,
instead returns the individual pieces of this plot (data points only,
fit lines only, and data points plus fit lines for individual random
levels) for a single species, plotted separately; `species` must then
resolve to exactly one species level, or a clear error is raised.

- `log.scale` Logical, if TRUE, plots rates and counts on a log scale;
  if FALSE, plots raw rates and counts.

- `CI_style` Logical, if TRUE, plots predicted lines with confidence
  interval style shading; if FALSE, plots predicted lines as solid lines
  with random levels and data points. Ignored when `decompose = TRUE`.

- `dim.boundaries` Numeric vector, independent block boundaries to plot
  for comparison. If empty, the argument is ignored.

- `y.lim` Numeric vector of length 2, limits for the y-axis of the
  plots. If NA (plain ratecount) or NULL (`decompose = TRUE`), defaults
  to automatic limits.

- `count.alpha.none` Numeric, transparency for count points when random
  level is "none". If left NA, defaults to 0.25. Ignored when
  `decompose = TRUE`.

- `count.alpha.ran` Numeric, transparency for count points when random
  level is not "none". If left NA, defaults to 0.25. Ignored when
  `decompose = TRUE`.

- `pred.alpha.none` Numeric, transparency for predicted lines when
  random level is "none". If left NA, defaults to 1.0. Ignored when
  `decompose = TRUE`.

- `pred.alpha.ran` Numeric, transparency for predicted lines when random
  level is not "none". If left NA, defaults to 0.9. Ignored when
  `decompose = TRUE`.

- `rans.to.print` Character vector, list of random levels to include on
  each species plot. If c(), all random levels are included. Ignored
  when `decompose = TRUE`.

Returns a list of ggplot objects for rate-count plots (or, with
`decompose = TRUE`, for the decomposed rate-count plots).

## timeseries

Plots timeseries effects from wisp results, separated by a splitting
factor (e.g., hemisphere).

- `splitting_factor` Character string, the name of the fixed effect by
  which to split the timeseries. If empty, c(), the first non-timeseries
  fixed effect is used. If "none", will not split.

- `tertiary_effects` Character vector, or NULL (default). Any fixed
  effects besides `timeseries` and the splitting factor are pruned down
  to a single level before rates are reconstructed (otherwise treatments
  belonging to different levels of those extra fixed effects would be
  mixed together). If NULL, each such fixed effect is pruned down to its
  reference level. Otherwise, for each additional fixed effect, if one
  of the levels named in `tertiary_effects` matches one of its levels,
  that level is used instead of the reference level.

- `log.scale` Logical, if TRUE, plots rates and counts on a log scale;
  if FALSE, plots raw rates and counts.

Returns a list of ggplot objects for timeseries plots.

## parameters

Makes nicely formatted violin (or bar) plots of the fitted wisp
parameters, including confidence intervals if stat information is
available.

- `violin` Logical, if TRUE, plots violin plots for each parameter; if
  FALSE, uses bar plots.

- `species.classes` List, a list of character vectors specifying which
  species levels to include together in plots. If NULL, all species
  levels are plotted separately.

- `mc_type` Character string, type of model component to plot ("rate",
  "tpoint", "tslope"). If NULL, all model components are plotted.

Returns a list of ggplot objects for parameter plots.

## effect.dist

Makes nicely formatted histograms of fitted parameters from wisp
results. Has no `species` argument; `species` is ignored (with a
warning) if set to a non-default value for this type. Returns a list of
ggplot objects containing histograms of fitted parameters.

## species.summary

Summarizes all important information (rate-count, residual, and
parameter plots) for one or more individual species levels, each on one
combined plot. Note that
[`wisp()`](https://michaelbarkasi.github.io/wispack/reference/wisp.md)
does not store built plot objects on its return value, so this always
rebuilds the rate-count, residual, and parameter plots it needs fresh
from `x` rather than reading them back from a stored `$plots` field.
Unlike every other `type`, this one has no `print.plots` argument of its
own and always prints its combined plots to the current graphics device;
the dispatcher's `print.plots` argument has no effect here.

- `contexts` Character vector, optional, specifies which context levels
  to summarize. Defaults to all.

Returns nothing; plots are printed to the current graphics device.

## MCMC.walks

Makes nicely formatted plots of the random walks and (normalized)
negative log likelihood sampled during MCMC simulation. Returns `NULL`
early if the model has no MCMC sampling (i.e. `x$sample.params.MCMC` is
empty/NULL). Has no `species` argument; `species` is ignored (with a
warning) if set to a non-default value for this type.

- `low_samples` Integer, number of low-value parameters to plot.
  Defaults to 10.

Returns a list of ggplot objects containing plots of walks for low-value
parameters, high-value parameters, and normalized negative log
likelihood (or `NULL`, see above).

## MCMC.bs.comparison

Allows for visual comparison of the parameter estimates from
bootstrapping vs MCMC simulation. Shows density for ten representative
parameters from bootstrapping and MCMC walk, the distributions of
Shapiro-Wilk p-values for bootstrapping vs MCMC walks, and the density
of autocorrelation results for bootstrapping vs MCMC walks. Has no
`species` argument; `species` is ignored (with a warning) if set to a
non-default value for this type.

- `ran.seed` Integer, random seed for reproducibility. Defaults to 1234.
  If NULL, no seed set.

Returns a list of ggplot objects containing plots of representative
parameter distributions, Shapiro-Wilk normality test results, and
autocorrelation plots for bootstrap and MCMC parameter estimates.

## Examples

``` r
set.seed(1)
n_bins <- 20
bin <- sample(seq_len(n_bins), 400, replace = TRUE)
lambda <- rep(c(2, 8), each = n_bins / 2)
countdata <- data.frame(bin = bin, count = rpois(400, lambda[bin]))
model <- wisp(
  countdata,
  fit_only       = TRUE,
  model.settings = list(LRO_cost = "none"),
  plot.settings  = list(print.plots = FALSE),
  verbose        = FALSE
)
rc <- plot(model, type = "ratecount", print.plots = FALSE, verbose = FALSE)
names(rc)
#> [1] "plot_pred_context_context_fixEff_species"
```
