# Wisp warping function, implemented in R

This function provides an R implementation of the warping function used
to warp model components before input into `wisp.sigmoid`.

## Usage

``` r
wisp.warp(z, b, w)
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
