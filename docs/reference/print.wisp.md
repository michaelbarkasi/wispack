# Print a fitted wisp model

Terse console summary of a fitted
[`wisp`](https://michaelbarkasi.github.io/wispack/reference/wisp.md)
model: number of fitted parameters, the grouping levels
(context/species/random-effect) the model was fit over, and which
resampling method (bootstraps and/or MCMC) was used to estimate
parameter uncertainty. This is deliberately terse – use
[`summary.wisp`](https://michaelbarkasi.github.io/wispack/reference/summary.wisp.md)
for the full parameter table, or
[`plot.wisp`](https://michaelbarkasi.github.io/wispack/reference/plot.wisp.md)
to visualize the fit.

## Usage

``` r
# S3 method for class 'wisp'
print(x, ...)
```

## Arguments

- x:

  List, output of the
  [`wisp`](https://michaelbarkasi.github.io/wispack/reference/wisp.md)
  function (an object of class `"wisp"`).

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
print(model)
#> Fitted wisp model
#> -----------------
#> Parameters:    4
#> Grouping:      1 context level(s), 1 species level(s), 1 random-effect level(s)
#> Resampling:    none (fit_only run; no parameter uncertainty estimated)
#> Stats:         not computed (see sample_stats())
#> 
#> See summary() for the full parameter table, or plot() to visualize the fit.
```
