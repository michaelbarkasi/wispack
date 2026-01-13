# Make plot demonstrating how the wisp function is warped by the warping factors

This function takes user-provided values for wisp function parameters
and produces a panel of plots showing how the wisp model is warped by
the warping factors

## Usage

``` r
demo.warp.plots(
  w = 2,
  point_pos = 60,
  point_neg = 40,
  Rt = c(6, 3, 0.2, 6) * 4.65,
  tslope = c(0.4, 0.75, 1),
  tpoint = c(15, 38, 80),
  w_factors = c(0.6, -0.9, 0.5)
)
```

## Arguments

- w:

  Numeric, warping factor (must be a single value).

- point_pos:

  Numeric, x coordinate at which to place positive warp segment (must be
  a single value).

- point_neg:

  Numeric, x coordinate at which to place negative warp segment (must be
  a single value).

- Rt:

  Numeric vector, rate parameters for the wisp function. Degree of the
  wisp model will be length of this vector minus 1.

- tslope:

  Numeric vector, slope scalars for the wisp function. Must be one less
  than the length of Rt.

- tpoint:

  Numeric vector, transition points for the wisp function. Must be one
  less than the length of Rt.

- w_factors:

  Numeric vector, warping factors for the wisp function. Must be
  length 3. First element is the warping factor for Rt, second for
  tslope, and third for tpoint. Note that this is more restrictive than
  the real wisp model, which not only allows for different warping
  factors across the model components (Rt, tslope, and tpoint), but also
  across the different elements within each model component as well.

## Value

Nothing. A ggplot object is printed to console.
