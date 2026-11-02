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

The returned vector mixes several different kinds of quantities,
distinguishable by name (matching `object$param.names`; see
`model.component.list` in the
[`wisp`](https://michaelbarkasi.github.io/wispack/reference/wisp.md)
return value for the three base component names). `"Rt"` entries are the
fitted rate on the log(sum + 1) scale of the summed count within a bin,
not the original per-instance rate – to recover the latter, exponentiate
the value, subtract one, and divide by the mean number of data rows per
bin (see the
[`wisp`](https://michaelbarkasi.github.io/wispack/reference/wisp.md)
examples, or the "Poisson processes and sigmoids" tutorial, for the full
back-transformation). `"tslope"` and `"tpoint"` entries are already on
their natural, directly interpretable scales and need no such
transformation. Whenever the random-effect grouping variable (`ran`) has
more than one level, the vector additionally contains `"wfactor_rate"`,
`"wfactor_slope"`, and `"wfactor_point"` entries (one per random-effect
level) – these are warping factors feeding into the warp function
described in the "Random effects and warping" tutorial, a different kind
of quantity altogether from the like-named base parameters; neither the
`"Rt"` back-transformation above nor "use as-is" applies to them.

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
