# Function to analyze benchmark results from attractor-based simulations

This function analyzes the results from benchmarking modeling methods on
attractor-based simulation datasets. It computes the performance
metrics: correlation coefficients for rate and random effect parameters,
false positive rates (FPR), false discovery rates (FDR), and power for
spatially variable gene (SVG) and functional spatial effect (FSE)
parameters. FPR, FDR, and power are defined as follows:

- FP: *p*-value below \\\alpha\\ for true null effect

- TN: *p*-value above \\\alpha\\ for true null effect

- TP: *p*-value below \\\alpha\\ for true non-null effect

- FN: *p*-value above \\\alpha\\ for true non-null effect

- FPR: FP / (FP + TN)

- FDR: FP / (FP + TP)

- Power: TP / (TP + FN)

The exact interpretation will depend on what's returned by the modeling
functions used. For example, `model_attractor_simulation_wisp` returns
*p*-values for the rate-effect parameter only and returns it under the
"FSE" parameter type, so FPR, FDR, and power for the "FSE" parameter
type will be computed based on the rate-effect estimates only.

## Usage

``` r
analyze_attractor_sim_benchmarks(results, sig_thresh = list(wisp = 0.05))
```

## Arguments

- results:

  The output of the `run_attractor_sim_benchmarks` function, either in
  its native form as a list, or converted into a data frame.

- sig_thresh:

  A named list specifying the significance threshold for each modeling
  method when computing FPR, FDR, and power. Default is a list with a
  single entry for the "wisp" method set to 0.05. If a value is not
  specified for a method found in results, will use 0.05.

## Value

A data frame summarizing the performance metrics for each modeling
method.

## Examples

``` r
# A toy results/times list in the same shape run_attractor_sim_benchmarks() returns.
toy_results <- data.frame(
  est    = c(0.1, 0.05, 0.02, 0.06),
  true   = c(FALSE, FALSE, TRUE, TRUE),
  param  = "FSE",
  id     = c("geneA", "geneB", "geneA", "geneB"),
  method = "wisp",
  sim    = c(1, 1, 2, 2)
)
toy_times <- data.frame(data_simulation = c(0.1, 0.1), wisp = c(0.5, 0.5))
analyze_attractor_sim_benchmarks(list(results = toy_results, times = toy_times))
#>                   wisp
#> cor_rate            NA
#> cor_ran             NA
#> FPR_svg             NA
#> FDR_svg             NA
#> power_svg           NA
#> FPR_fse            0.0
#> FDR_fse            0.0
#> power_fse          0.5
#> mean_sec_to_model  0.5
```
