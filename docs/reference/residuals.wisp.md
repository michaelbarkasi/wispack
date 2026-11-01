# Residuals from a fitted wisp model

Returns the residual summary statistics (mean, sd, variance, by group)
computed by
[`analyze_residuals`](https://michaelbarkasi.github.io/wispack/reference/analyze_residuals.md)
and stored on the fitted model, either for raw residuals or
log-residuals.

## Usage

``` r
# S3 method for class 'wisp'
residuals(object, type = c("raw", "log"), ...)
```

## Arguments

- object:

  List, output of the
  [`wisp`](https://michaelbarkasi.github.io/wispack/reference/wisp.md)
  function (an object of class `"wisp"`).

- type:

  Character string, one of `"raw"` (default) or `"log"`, specifying
  whether to return the summary statistics for raw residuals
  (`object$stats$residuals`) or log-residuals
  (`object$stats$residuals.log`).

- ...:

  Ignored; present only for S3 method consistency with
  [`residuals`](https://rdrr.io/r/stats/residuals.html).

## Value

A data frame of residual summary statistics by group (columns `group`,
`mean`, `sd`, `variance`), or an error if statistical analysis was never
run on this model (e.g. it was fit with `fit_only = TRUE`).

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
residuals(model)
#>                  group     mean       sd variance
#> 1                  all 3.725845 7.988364 63.81396
#> 2         ran_lvl_none 4.033232 6.461292 41.74830
#> 3           ran_lvl_r1 4.637785 9.219901 85.00657
#> 4           ran_lvl_r2 2.506519 8.092267 65.48478
#> 5        treatment_ref 3.223797 6.423804 41.26526
#> 6         treatment_18 4.227894 9.323712 86.93161
#> 7 demoContext_demoGene 3.725845 7.988364 63.81396
residuals(model, type = "log")
#>                  group       mean        sd   variance
#> 1                  all 0.10081184 0.2049539 0.04200609
#> 2         ran_lvl_none 0.10238083 0.1616210 0.02612136
#> 3           ran_lvl_r1 0.10026499 0.2129685 0.04535560
#> 4           ran_lvl_r2 0.09978970 0.2380999 0.05669157
#> 5        treatment_ref 0.10734713 0.2110623 0.04454727
#> 6         treatment_18 0.09427655 0.2002249 0.04009001
#> 7 demoContext_demoGene 0.10081184 0.2049539 0.04200609
```
