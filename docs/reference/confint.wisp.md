# Confidence intervals for wisp model parameters

Returns confidence intervals for the fitted wisp model parameters. These
were already computed once, at whatever significance level (`alpha`) was
in effect when
[`sample_stats`](https://michaelbarkasi.github.io/wispack/reference/sample_stats.md)
last ran on this model (the default, both for
[`sample_stats()`](https://michaelbarkasi.github.io/wispack/reference/sample_stats.md)
and for the statistical analysis
[`wisp`](https://michaelbarkasi.github.io/wispack/reference/wisp.md)
runs internally, is `alpha = 0.05`, i.e. a 95\\ columns of
`object$stats$parameters`.

## Usage

``` r
# S3 method for class 'wisp'
confint(object, parm = NULL, level = 0.95, ...)
```

## Arguments

- object:

  List, output of the
  [`wisp`](https://michaelbarkasi.github.io/wispack/reference/wisp.md)
  function (an object of class `"wisp"`).

- parm:

  Character vector, or NULL (default). Names of parameters to include
  (matching `object$param.names`); if NULL, confidence intervals for all
  parameters are returned. Follows the same `parm` semantics as
  [`confint`](https://rdrr.io/r/stats/confint.html).

- level:

  Numeric, the confidence level required. Defaults to 0.95, matching the
  confidence level
  [`wisp`](https://michaelbarkasi.github.io/wispack/reference/wisp.md)/[`sample_stats`](https://michaelbarkasi.github.io/wispack/reference/sample_stats.md)
  compute by default; see Details above for what happens when a
  different level is requested.

- ...:

  Additional arguments passed on to
  [`sample_stats`](https://michaelbarkasi.github.io/wispack/reference/sample_stats.md)
  when `level != 0.95` (e.g. `Bonferroni`, `conv.resamples.only`,
  `verbose`).

## Value

A matrix with one row per (selected) parameter and two columns, `CI.low`
and `CI.high`, row-named by parameter.

## Details

If `level = 0.95` (the default, matching
[`sample_stats()`](https://michaelbarkasi.github.io/wispack/reference/sample_stats.md)'s
own default `alpha = 0.05`), the already-stored `CI.low`/`CI.high`
columns are returned directly – exactly as they were computed during the
original fit (respecting whatever `Bonferroni`/`conv.resamples.only`
choices were made then, even if those differ from
[`sample_stats()`](https://michaelbarkasi.github.io/wispack/reference/sample_stats.md)'s
own defaults). If a different `level` is requested,
[`sample_stats`](https://michaelbarkasi.github.io/wispack/reference/sample_stats.md)
is transparently re-run (with `alpha = 1 - level`, and otherwise its own
defaults for `Bonferroni` and `conv.resamples.only`, since this model
does not record which choices were used for its originally-stored
statistics) to recompute the intervals at the requested level; this is
cheap, since it only re-quantiles the resampled parameters already
stored on the model; it does not refit anything.

## Examples

``` r
set.seed(1)
grid <- expand.grid(bin = 1:20, age = c(12, 18), ran = c("r1", "r2"), rep = 1:6)
base_rate <- ifelse(grid$bin <= 10, 3, 9)
age_effect <- ifelse(grid$age == 18, 1.5, 1)
grid$count <- rpois(nrow(grid), base_rate * age_effect)
countdata <- data.frame(
  bin = grid$bin, count = grid$count, age = grid$age, ran = grid$ran,
  species = "demoGene", context = "demoContext"
)
model <- wisp(
  countdata,
  variables      = list(timeseries = "age"),
  model.settings = list(LRO_cost = "none"),
  plot.settings  = list(print.plots = FALSE),
  MCMC.settings  = list(MCMC.steps = 100, MCMC.burnin = 20),
  bootstraps.num = 10,
  max.fork       = 1,
  verbose        = FALSE
)
confint(model)
#>                                                      CI.low     CI.high
#> baseline_demoContext_Rt_demoGene_Tns/Blk1       3.016476952  3.12833343
#> beta_Rt_demoContext_demoGene_18_X_Tns/Blk1      0.213825918  0.41373890
#> baseline_demoContext_Rt_demoGene_Tns/Blk2       4.075587853  4.16161032
#> beta_Rt_demoContext_demoGene_18_X_Tns/Blk2      0.314764126  0.39523704
#> baseline_demoContext_tslope_demoGene_Tns/Blk1   5.519569538 14.52583539
#> beta_tslope_demoContext_demoGene_18_X_Tns/Blk1  1.767844600  4.83561112
#> baseline_demoContext_tpoint_demoGene_Tns/Blk1  10.171114279 10.76825512
#> beta_tpoint_demoContext_demoGene_18_X_Tns/Blk1 -0.125822886  0.52492211
#> wfactor_point_r1_X_demoGene                    -0.042736430  0.30219300
#> wfactor_rate_r1_X_demoGene                     -0.020552776  0.11270932
#> wfactor_slope_r1_X_demoGene                     0.006024845  0.26510370
#> wfactor_point_r2_X_demoGene                    -0.233002271  0.02874054
#> wfactor_rate_r2_X_demoGene                     -0.114990011  0.01830464
#> wfactor_slope_r2_X_demoGene                     0.057493081  2.14208315
confint(model, level = 0.90)
#>                                                      CI.low     CI.high
#> baseline_demoContext_Rt_demoGene_Tns/Blk1       3.024515476  3.12277975
#> beta_Rt_demoContext_demoGene_18_X_Tns/Blk1      0.215707296  0.40962098
#> baseline_demoContext_Rt_demoGene_Tns/Blk2       4.078552713  4.15504201
#> beta_Rt_demoContext_demoGene_18_X_Tns/Blk2      0.316450735  0.39504406
#> baseline_demoContext_tslope_demoGene_Tns/Blk1   5.528391381 14.39380198
#> beta_tslope_demoContext_demoGene_18_X_Tns/Blk1  1.779599603  4.75791994
#> baseline_demoContext_tpoint_demoGene_Tns/Blk1  10.210434661 10.75089289
#> beta_tpoint_demoContext_demoGene_18_X_Tns/Blk1 -0.121094838  0.50167546
#> wfactor_point_r1_X_demoGene                    -0.031560061  0.28558263
#> wfactor_rate_r1_X_demoGene                     -0.005847111  0.11047191
#> wfactor_slope_r1_X_demoGene                     0.011613842  0.26299121
#> wfactor_point_r2_X_demoGene                    -0.226639689  0.01865384
#> wfactor_rate_r2_X_demoGene                     -0.109331077  0.01157228
#> wfactor_slope_r2_X_demoGene                     0.069949611  1.88225845
```
