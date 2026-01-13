# Wisp sigmoid function, implemented in R

This function provides an R implementation of the wisp sigmoid function.
It assumes the parameters Rt, tslope, and tpoint have already been
warped by `wisp.warp`.

## Usage

``` r
wisp.sigmoid(x, Rt, tslope, tpoint)
```

## Arguments

- x:

  Numeric, input spatial position, either a scalar, vector, or 2D
  array/matrix

- Rt:

  Numeric vector, rate parameters for the wisp function. Degree of the
  wisp model will be length of this vector minus 1.

- tslope:

  Numeric vector, slope scalars for the wisp function. Must be one less
  than the length of Rt.

- tpoint:

  Numeric vector, transition points for the wisp function. Must be one
  less than the length of Rt.

## Value

Numeric vector or matrix, wisp sigmoid values for given input. Returned
object will have the same dimensions as input x.
