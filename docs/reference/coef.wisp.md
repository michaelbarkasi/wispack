# Fitted parameter estimates from a wisp model

Returns the fitted wisp model parameter estimates as a plain named
numeric vector, following the usual R convention for
[`coef()`](https://rdrr.io/r/stats/coef.html) methods.

## Usage

``` r
# S3 method for class 'wisp'
coef(object, ...)
```

## Arguments

- object:

  List, output of the
  [`wisp`](https://michaelbarkasi.github.io/wispack/reference/wisp.md)
  function (an object of class `"wisp"`).

- ...:

  Ignored; present only for S3 method consistency with
  [`coef`](https://rdrr.io/r/stats/coef.html).

## Value

A named numeric vector of fitted parameter estimates (names are the
parameter names, matching `object$param.names`).

## Note

`object$fitted.parameters` is used directly, rather than
`object$stats$parameters$estimate`, as the source: it's already a named
numeric vector (no extra reshaping needed to meet the
[`coef()`](https://rdrr.io/r/stats/coef.html) convention), and – unlike
`object$stats$parameters$estimate` – it's always present, even for a
model fit with `fit_only = TRUE` where no statistical analysis was ever
run. The two are otherwise identical.

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
coef(model)
#>     baseline_context_Rt_species_Tns/Blk1 
#>                                 3.890853 
#>     baseline_context_Rt_species_Tns/Blk2 
#>                                 5.061927 
#> baseline_context_tslope_species_Tns/Blk1 
#>                                 7.494465 
#> baseline_context_tpoint_species_Tns/Blk1 
#>                                10.699851 
```
