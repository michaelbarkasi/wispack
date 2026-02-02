# Make plot demonstrating how the wisp function is determined by the Rt, slope, and tpoint parameters

This function takes user-provided values for wisp function parameters
and produces a panel of plots showing how the wisp model is determined
by the Rt, slope, and tpoint parameters.

## Usage

``` r
demo.sigmoid.plots(
  r = 4,
  s = 1,
  Rt = c(6, 3, 0.2, 6) * 4.65,
  tslope = c(0.4, 0.75, 1),
  tpoint = c(15, 38, 80),
  return_plots = FALSE
)
```

## Arguments

- r:

  Numeric, upper asymptote for logistic (must be a single value).

- s:

  Numeric, slope scalar at inflection point (must be a single value).

- Rt:

  Numeric vector, rate parameters for the wisp function. Degree of the
  wisp model will be length of this vector minus 1.

- tslope:

  Numeric vector, slope scalars for the wisp function. Must be one less
  than the length of Rt.

- tpoint:

  Numeric vector, transition points for the wisp function. Must be one
  less than the length of Rt.

- return_plots:

  Logical, if true, returns two separate plots; otherwise, returns
  nothing and prints combined plot

## Value

Nothing. A ggplot object is printed to console.
