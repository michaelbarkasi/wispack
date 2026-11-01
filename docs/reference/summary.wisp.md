# Summarize a fitted wisp model

Builds an object of class `"summary.wisp"`, primarily wrapping the
parameter table computed by
[`sample_stats`](https://michaelbarkasi.github.io/wispack/reference/sample_stats.md)
(`x$stats$parameters`), plus the same model-size/fit-setting context
shown by
[`print.wisp`](https://michaelbarkasi.github.io/wispack/reference/print.wisp.md).
Use
[`print.summary.wisp`](https://michaelbarkasi.github.io/wispack/reference/print.summary.wisp.md)
(dispatched automatically by
[`print()`](https://rdrr.io/r/base/print.html)) to display the result.

## Usage

``` r
# S3 method for class 'wisp'
summary(object, ...)
```

## Arguments

- object:

  List, output of the
  [`wisp`](https://michaelbarkasi.github.io/wispack/reference/wisp.md)
  function (an object of class `"wisp"`).

- ...:

  Ignored; present only for S3 method consistency with
  [`summary`](https://rdrr.io/r/base/summary.html).

## Value

An object of class `"summary.wisp"`: a list with elements `parameters`
(the parameter table from `object$stats$parameters`, or `NULL` if
statistical analysis was never run – e.g. the model was fit with
`fit_only = TRUE`), `n.context`, `n.species`, `n.ran`,
`sample.params.type`, `n.bs`, and `n.mcmc`.

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
summary(model)
#> Summary of fitted wisp model
#> -----------------------------
#> Grouping: 1 context level(s), 1 species level(s), 2 random-effect level(s)
#> Resampling: 11 bootstrap resample(s), 101 MCMC step(s) retained (sample.params.type = "bs")
#> 
#> 
#> Parameter estimates:
#> ------------------------------
#>                                         parameter    estimate       CI.low
#> 1       baseline_demoContext_Rt_demoGene_Tns/Blk1  3.00843843  3.016476952
#> 2      beta_Rt_demoContext_demoGene_18_X_Tns/Blk1  0.41785682  0.213825918
#> 3       baseline_demoContext_Rt_demoGene_Tns/Blk2  4.09250468  4.075587853
#> 4      beta_Rt_demoContext_demoGene_18_X_Tns/Blk2  0.39311418  0.314764126
#> 5   baseline_demoContext_tslope_demoGene_Tns/Blk1  5.51074770  5.519569538
#> 6  beta_tslope_demoContext_demoGene_18_X_Tns/Blk1  2.10216909  1.767844600
#> 7   baseline_demoContext_tpoint_demoGene_Tns/Blk1 10.56853427 10.171114279
#> 8  beta_tpoint_demoContext_demoGene_18_X_Tns/Blk1 -0.04871231 -0.125822886
#> 9                     wfactor_point_r1_X_demoGene  0.15730675 -0.042736430
#> 10                     wfactor_rate_r1_X_demoGene  0.08240972 -0.020552776
#> 11                    wfactor_slope_r1_X_demoGene  0.03297546  0.006024845
#> 12                    wfactor_point_r2_X_demoGene -0.18073973 -0.233002271
#> 13                     wfactor_rate_r2_X_demoGene -0.08321253 -0.114990011
#> 14                    wfactor_slope_r2_X_demoGene  0.58886926  0.057493081
#>        CI.high p.value p.value.adj  alpha.adj significance
#> 1   3.12833343      NA          NA         NA             
#> 2   0.41373890       0           0 0.01250000          ***
#> 3   4.16161032      NA          NA         NA             
#> 4   0.39523704       0           0 0.01666667          ***
#> 5  14.52583539      NA          NA         NA             
#> 6   4.83561112       0           0 0.02500000          ***
#> 7  10.76825512      NA          NA         NA             
#> 8   0.52492211       1           1 0.05000000           ns
#> 9   0.30219300      NA          NA         NA             
#> 10  0.11270932      NA          NA         NA             
#> 11  0.26510370      NA          NA         NA             
#> 12  0.02874054      NA          NA         NA             
#> 13  0.01830464      NA          NA         NA             
#> 14  2.14208315      NA          NA         NA             
#> ----
```
