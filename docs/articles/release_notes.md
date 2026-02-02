# Release Notes

## Planned changes for future releases

- Autodetect fixed-effect type? E.g., fixed-effect columns of R type
  “character” or “factor” are treated as two-level categorical factors,
  and numeric columns are treated as discrete additive effects (e.g.,
  like a time series). Will need to update the RORB tutorial and any
  other tutorials which discuss time-series modeling.
- Ensure any columns flagged as a time series are also treated as
  fixed-effects.
- Finish vignettes on statistical analyses and random-effect handling.
- Add a check during the fit to see if any transition points can be
  replaced with a slope.
- Add variable check to the loading of plot.settings in wisp. Currently
  only runs check_list without doing a more substantive check. Compare
  with the more detailed checks of MCMC.settings.
- Allow MCMC.prior to be a vector of different prior distributions for
  each parameter. When completed, update the tutorial on customizing
  statistical analyses.

## Version 1.1

- Added discrete time-series modeling functionality (the **timeseries**
  variable option) and plotting (the function plot.timeseries()).
- Ensured code would robustly run for any data including at least
  **count** and **bin** columns, without need for **context**,
  **species**, **ran**, or fixed-effect variables.
- Added explicit **fit_only** option to wisp() to avoid running any
  parameter estimation (MCMC or bootstrapping).

## Version 1.0

- Initial public release of wispack, as used in [this
  preprint](https://doi.org/10.1101/2025.06.11.659209).
- Defines the wisp() function for implementing wisps.
- Introduces one-dimensional warped sigmoidal Poisson-process
  mixed-effect modeling (one-dimensional wisps) for testing for
  functional spatial effects in spatial transcriptomics data.
