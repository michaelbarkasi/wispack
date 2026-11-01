# Wisp sigmoid function, implemented in R

This function provides an R implementation of the wisp sigmoid function.
It assumes the parameters Rt, tslope, and tpoint have already been
warped by `wisp_warp`.

## Usage

``` r
wisp_sigmoid(x, Rt, tslope, tpoint)
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

## Examples

``` r
# A single rate block (degree 0) is just a constant
wisp_sigmoid(c(1, 5, 10), Rt = 3, tslope = numeric(0), tpoint = numeric(0))
#> [1] 3 3 3

# Two rate blocks separated by one transition point: far from the transition
# point, the sigmoid recovers each block's asymptotic rate
wisp_sigmoid(c(-1e3, 1e3), Rt = c(1, 9), tslope = 1, tpoint = 0)
#> [1] 1 9
```
