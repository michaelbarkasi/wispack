# Wisp warping function, implemented in R

This function provides an R implementation of the warping function used
to warp model components before input into `wisp_sigmoid`.

## Usage

``` r
wisp_warp(z, b, w)
```

## Arguments

- z:

  Numeric vector or matrix, values to warp. Scalar values fine as well.
  Matrix must be 2D.

- b:

  Numeric, warping bound (must be a single value).

- w:

  Numeric, warping factor (must be a single value).

## Value

Numeric vector or matrix, warped values. Returned object will have the
same dimensions as input z.

## Examples

``` r
# No warping (w = 0) returns the input unchanged
wisp_warp(c(1, 5, 10), b = 10, w = 0)
#> [1]  1  5 10

# Positive warping factor compresses values toward the bound b
wisp_warp(c(1, 5, 10), b = 10, w = 2)
#> [1]  3.967120  9.323324 10.000000

# Works on a matrix too, preserving its dimensions
wisp_warp(matrix(1:4, 2, 2), b = 10, w = 1)
#>          [,1]     [,2]
#> [1,] 1.856463 4.814272
#> [2,] 3.450154 5.978080
```
