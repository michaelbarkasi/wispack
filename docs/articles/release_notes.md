# Release Notes

## Planned changes for future releases

- Experiment with extending the treatment knock-out (trtKO)
  functionality into a likelihood-ratio test for hypothesis testing.
  Update the tutorial on customizing statistical analyses, if needed.
- Autodetect fixed-effect type, e.g., fixed-effect columns of R type
  “character” or “factor” are treated as two-level categorical factors,
  and numeric columns are treated as discrete additive effects (e.g.,
  like a time series). Will need to update the RORB tutorial and any
  other tutorials which discuss time-series modeling.
- Ensure any columns flagged as a time series are also treated as
  fixed-effects.
- Add a check during the fit to see if any transition points can be
  replaced with a slope.
- Add variable check to the loading of plot.settings in wisp. Currently
  only runs check_list without doing a more substantive check. Compare
  with the more detailed checks of MCMC.settings.
- Allow MCMC.prior to be a vector of different prior distributions for
  each parameter. When completed, update the tutorial on customizing
  statistical analyses.
- Ensure nested loops for matrix operations have the major axis (column)
  on outside.
- Extracellular (out-of-soma) transcripts do not follow a Poisson
  distribution, and so any attempt to model them with wisp needs to
  account for this.

## Version 2.6.1 (November 2, 2026)

- plot.wisp’s type argument now defaults to “ratecount”, so plot(model)
  no longer requires naming a type.
- Fixed plot.wisp silently ignoring the plot.settings a model was
  originally fit with when later re-plotting it: type = “ratecount” and
  type = “timeseries” now default their own log.scale, CI_style,
  dim.boundaries, count.alpha.ran/count.alpha.none,
  pred.alpha.ran/pred.alpha.none, and splitting_factor arguments to the
  model’s own plot.settings (as configured in the original wisp() call),
  rather than to a fixed value – matching how type = “species.summary”
  already behaved, and how every type’s cosmetic sizing/color settings
  already worked. These arguments can still be overridden on the plot()
  call itself. **Note:** a call that relied on the old behavior (a
  non-default plot.settings being ignored on manual re-plotting) will
  now see that setting applied by default; pass the argument explicitly
  to plot() to keep the old output.
- Fixed analyze_residuals’s description, which had been copy-pasted from
  an unrelated internal helper (pvalues.samples) and incorrectly claimed
  the function computes p-values from an empirical CDF; it does not – it
  computes and summarizes raw/log residuals by group. Also expanded its
  documented return value with the group and plot-naming structure.
- Expanded the wisp example, and the coef.wisp/confint.wisp
  documentation, to explain and demonstrate how to recover the original
  per-instance rate from a fitted Rt parameter (exponentiate, subtract
  one, divide by the mean number of data rows per bin – see the “Poisson
  processes and sigmoids” tutorial), and to distinguish it from the
  tslope/tpoint parameters (already on their natural scale) and the
  wfactor_rate/wfactor_slope/wfactor_point random-effect warping-factor
  parameters (present only when the random-effect factor has more than
  one level; a different kind of quantity from the like-named base
  parameters, not a rescaled version of one).
- README.md is now generated from a new README.Rmd (via
  devtools::build_readme()), carrying the same self-updating “date
  (vX.Y)” build stamp as the other articles, instead of a hand-typed
  string that had already drifted out of date.
- Fixed the peak-memory table under “Result 2” of the forking article,
  which was missing its 20-bootstraps/20-workers row – the data point
  the following paragraph’s “holding workers fixed at 20, doubling
  bootstraps from 20 to 40” comparison actually depends on.

## Version 2.6.0 (October 29, 2026)

- **Breaking API change**: the eight dotted plotting functions
  (plot.ratecount, plot.timeseries, plot.parameters, plot.effect.dist,
  plot.species.summary, plot.MCMC.walks, plot.decomposition,
  plot.MCMC.bs.comparison) have been removed. They looked like S3
  methods but never were one, since wisp() used to return a plain,
  unclassed list, so plot() could never dispatch to them; users had to
  call each by its full dotted name. wisp() now classes its return value
  “wisp”, and all plotting goes through one real, exported S3 method,
  plot.wisp, selected with a type argument. All other arguments to each
  old function are unchanged and are simply passed through …. The
  old-name to new-call mapping is:
  - plot.ratecount(model, …) → plot(model, type = “ratecount”, …)
  - plot.timeseries(model, …) → plot(model, type = “timeseries”, …)
  - plot.parameters(model, …) → plot(model, type = “parameters”, …)
  - plot.effect.dist(model, …) → plot(model, type = “effect.dist”, …)
  - plot.species.summary(model, …) → plot(model, type =
    “species.summary”, …)
  - plot.MCMC.walks(model, …) → plot(model, type = “MCMC.walks”, …)
  - plot.decomposition(model, …) → plot(model, type = “decomposition”,
    …)
  - plot.MCMC.bs.comparison(model, …) → plot(model, type =
    “MCMC.bs.comparison”, …)
- plot.wisp’s print.plots, verbose, and species arguments, previously
  only available by name on the underlying per-type functions (and so
  only reachable through …), are now named formals on plot.wisp itself:
  plot.wisp(x, type, print.plots = TRUE, verbose = TRUE, species =
  “all”, decompose = FALSE, …). Two of these are **breaking changes**:
  - print.plots now defaults to TRUE at the plot.wisp level (every
    underlying function’s own old default was FALSE). A call that relied
    on the old default to stay silent (e.g. in a script or
    non-interactive context) must now pass print.plots = FALSE
    explicitly to keep the old behavior.
  - species now defaults to the string “all” (translated internally to
    c()) rather than c() directly; this is purely a more readable public
    spelling of the same default and changes nothing about which species
    are plotted. species is only forwarded to type %in% c(“ratecount”,
    “timeseries”, “parameters”, “species.summary”) (the types whose
    underlying function actually has a species argument); passing a
    non-default species for any other type now raises a warning and is
    otherwise ignored, rather than erroring with “unused argument” or
    being silently (and confusingly) absorbed by ….
  - verbose defaults to TRUE, matching every underlying function’s own
    existing default, so this is not a behavior change.
- **Breaking API change**: type = “decomposition” has been removed as a
  standalone plot.wisp type and folded into type = “ratecount” as a
  display mode, selected with a new decompose argument. plot(model, type
  = “decomposition”, species = “X”, …) becomes plot(model, type =
  “ratecount”, species = “X”, decompose = TRUE, …). As before, decompose
  = TRUE requires species to resolve to exactly one species level (not
  “all”, not a multi-element vector), or a clear error is raised;
  setting decompose = TRUE for any type other than “ratecount” is also
  an error. Unlike the old standalone type, the decomposed rate-count
  plots now honor the dispatcher’s shared print.plots and verbose
  arguments, just like every other type.
- Added five new, real, dispatchable S3 methods for the “wisp” class,
  for the same reason plot.wisp was made a real S3 method in the
  breaking change above (plain exportPattern-based registration is not
  reliably found by dispatch once colorspace, a Depends: of this
  package, turns a base generic into an S4 generic for the session – so
  each of these also carries an explicit @method roxygen tag and an
  explicit S3method() NAMESPACE entry):
  - print.wisp(x, …): terse console digest of a fitted model (parameter
    count, grouping levels, resampling method/counts, and a stats
    summary).
  - summary.wisp(object, …) / print.summary.wisp(x, …): summary()
    returns an object of class “summary.wisp” wrapping
    objectstatsparameters plus model-size/fit-setting context; print()
    formats it for console display.
  - residuals.wisp(object, type = c(“raw”, “log”), …): returns the
    residual summary statistics (by group) computed during the fit, for
    either raw residuals or log-residuals.
  - coef.wisp(object, …): returns the fitted parameter estimates as a
    plain named numeric vector, per the usual R coef() convention.
  - confint.wisp(object, parm = NULL, level = 0.95, …): returns the
    already-computed confidence intervals (stored at whatever alpha was
    used when sample_stats last ran, default 95%) directly when level =
    0.95; for any other level, transparently re-runs sample_stats
    (cheaply – it only re-quantiles the already-stored resampled
    parameters, without refitting) to recompute at the requested level.
  - All five give a clear, specific error (or message) rather than a
    cryptic one when called on a model fit with fit_only = TRUE (i.e. no
    statistical analysis was ever run).
- Renamed the two exported warping/sigmoid math utilities wisp.warp and
  wisp.sigmoid to wisp_warp and wisp_sigmoid (dot to underscore). This
  is a pure rename with no change to arguments or behavior; it removes
  the last two function names in the package that looked like (but never
  were) S3 methods for a wisp generic (wisp() is the model-fitting
  function itself, not a generic, and never will be one).

## Version 2.5.1 (October 3, 2026)

- Updated .Rbuildignore.

## Version 2.5 (October 3, 2026)

- Improved time series parsing.
- Updated time series plots to allow for tertiary effects.

## Version 2.4 (July 3, 2026)

- Part 2 of a major code refactor/optimization.
- Finished refactoring wspc object class and its initialization method
  into a more organized structure.
- Fixed issue in extrapolate_none.
- Added option to level preprocessing coordinate transform (the
  cortical_coordinate_transform function) based on different layers.
- Updated Roxygen2 descriptions for user-facing R functions.
- Checked all tutorials to ensure they work with current version and
  removed auto-dating, so their current working version is properly
  documented. Included footnotes on saved and pre-processed code chunks.
- Substantial code cleaning.

## Version 2.3 (June 16, 2026)

- Part 1 of a major code refactor/optimization.
- Improved algorithm in check_parameter_feasibility for finding nearby
  feasible initial parameters.
- Improved placement of check_parameter_feasibility within the various
  fit, MCMC, and bootstrap pipelines.
- Converted several large-signature stand-alone functions into wspc
  methods: compute_gamma_dispersion, make_extrapolation_pool,
  find_count_log_means, estimate_change_points, and
  estimate_initial_parameters.
- Simplified parameter vector handling.
- Wrote checks into tutorials, to ensure code changes do not change
  results.
- Replaced some R Lists with std::vector.
- Improved weight-matrix construction.
- Replaced pcg with std rng, to ease licensing issues.
- Removed stan variables from wspc objects.
- Replaced R-style masking with numeric indexing for summed count data.
- Added print.settings argument splitting_factor_colors to wisp to allow
  users to specify colors for the splitting factor in ratecount and
  timeseries plots.
- Performed extensive optimization of the code, leading to a x2 or x3
  speedup in running wisp.
- Documentation updated, but not yet rebuilt.

## Version 2.2 (June 5, 2026)

- Added model.settings argument max_bin to wisp which allows users to
  specify the largest bin number. If left at zero, the function will
  infer the maximum bin number from the count data.
- Aligned matrix nested loops to major axis (column) to improve
  computational efficiency.
- Updated plot.decomposition to allow for multiple contexts.
- Fixed plot.timeseries to plot single-block species.
- Added LRO parameter grid search option. (Very effective!)
- Updated preprocessing functions for making count data.
- Fixed batch size for forking.

## Version 2.1 (March 1, 2026)

- Added model.settings argument trtKO to wisp which enables excluding
  treatments to run reduced models.
- Replaced corrupted data file “corticallaminar_model.rds” with working
  version (hopefully for the last time!).
- Fixed memory leak related to CI computation.

## Version 2.0 (Feb 19, 2026)

- Release attached to NAR resubmission.
- Redesigned plots.
- Added attractor simulation functionality and benchmarking.
- Includes first-draft versions of package tutorials.
- Note that the package data file “corticallaminar_model.rds” was
  corrupted and should be pulled from the latest version or latest
  commit on the main branch.

## Version 1.1

- Added discrete time-series modeling functionality (the timeseries
  variable option) and plotting (the function plot.timeseries).
- Ensured code would robustly run for any data including at least count
  and bin columns, without need for context, species, ran, or
  fixed-effect variables.
- Added explicit fit_only option to wisp to avoid running any parameter
  estimation (MCMC or bootstrapping).

## Version 1.0

- Initial public release of wispack, as used in [this
  preprint](https://doi.org/10.1101/2025.06.11.659209).
- Defines the wisp function for implementing wisps.
- Introduces one-dimensional warped sigmoidal Poisson-process
  mixed-effect modeling (one-dimensional wisps) for testing for
  functional spatial effects in spatial transcriptomics data.
