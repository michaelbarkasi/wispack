# Function to analyze benchmark results from attractor-based simulations

This function analyzes the results from benchmarking modeling methods on
attractor-based simulation datasets. It computes the performance
metrics: correlation coefficients for rate and random effect parameters,
false positive rates (FPR), false discovery rates (FDR), and power for
spatially variable gene (SVG) and functional spatial effect (FSE)
parameters.

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
