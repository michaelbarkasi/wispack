# Print a wisp model summary

Formats a `"summary.wisp"` object (as returned by
[`summary.wisp`](https://michaelbarkasi.github.io/wispack/reference/summary.wisp.md))
for console display.

## Usage

``` r
# S3 method for class 'summary.wisp'
print(x, ...)
```

## Arguments

- x:

  An object of class `"summary.wisp"`, as returned by
  [`summary.wisp`](https://michaelbarkasi.github.io/wispack/reference/summary.wisp.md).

- ...:

  Ignored; present only for S3 method consistency with
  [`print`](https://rdrr.io/r/base/print.html).

## Value

`x`, invisibly (as is standard for `print` methods).

## Examples

``` r
set.seed(1)
n_bins <- 20
bin <- sample(seq_len(n_bins), 400, replace = TRUE)
lambda <- rep(c(2, 8), each = n_bins / 2)
countdata <- data.frame(bin = bin, count = rpois(400, lambda[bin]))
model <- wisp(
  countdata,
  fit_only       = TRUE,
  model.settings = list(LRO_cost = "none"),
  plot.settings  = list(print.plots = FALSE),
  verbose        = FALSE
)
# print() dispatches to print.summary.wisp() automatically for a "summary.wisp" object
print(summary(model))
#> Summary of fitted wisp model
#> -----------------------------
#> Grouping: 1 context level(s), 1 species level(s), 1 random-effect level(s)
#> Resampling: none (fit_only run; no parameter uncertainty estimated)
#> 
#> No parameter statistics are available: this model was likely fit with fit_only = TRUE (or sample_stats()  was otherwise never run). Call sample_stats() on the fitted model to compute parameter estimates,  confidence intervals, and p-values.
```
