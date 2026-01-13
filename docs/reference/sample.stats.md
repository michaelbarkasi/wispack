# Estimate p-values and confidence intervals from resampled parameters

Runs statistical analysis on the resampled parameters from the wisp
function. It computes p-values and confidence intervals for each
parameter, adjusting for multiple comparisons using either the
Bonferroni correction or the Holm-Bonferroni method.

## Usage

``` r
sample.stats(
  wisp.results,
  alpha = 0.05,
  Bonferroni = FALSE,
  conv.resamples.only = TRUE,
  verbose = TRUE
)
```

## Arguments

- wisp.results:

  List, output of the wisp function.

- alpha:

  Numeric value giving significance level for p-values and confidence
  intervals. Default is 0.05.

- Bonferroni:

  Logical, if TRUE, uses the Bonferroni correction for multiple
  comparisons; if FALSE, uses the Holm-Bonferroni method. Default is
  FALSE.

- conv.resamples.only:

  Logical, if TRUE, only resamples with a converged fit are used for
  statistical analysis; if FALSE, all resamples are used. Default is
  TRUE.

- verbose:

  Logical, if TRUE, prints information during the statistical analysis.

## Value

Data frame giving, for each parameter, its name, estimate, confidence
interval (CI.low, CI.high), p-value, adjusted p-value (p.value.adj),
adjusted alpha (alpha.adj), and significance level (significance).
