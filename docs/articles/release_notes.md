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
