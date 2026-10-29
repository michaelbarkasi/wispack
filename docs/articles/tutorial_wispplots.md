# Plotting wisp models

## Introduction

Wisp models can have thousands of parameters, each with their own
intuitive physical interpretation. Schematically, the picture is more
simple. At bottom, wisp models are parameterized by a handful of
[fixed-effect spatial parameter
types](https://michaelbarkasi.github.io/wispack/articles/tutorial_Poisson.md)
– block rates, transition points, and transition slopes – plus
[random-effect
scalars](https://michaelbarkasi.github.io/wispack/articles/tutorial_warping.md)
for each of these fixed-effect spatial parameter types. However, the
fixed-effect spatial parameters multiply with the number of species, the
number of
[contexts](https://michaelbarkasi.github.io/wispack/articles/tutorial_multicontext.md),
and the number of treatment conditions. The random-effect scalars
multiply by replicate count and species. The upshot is that a wisp model
includes a dizzying array of parameter estimates and hypothesis test
results. To help with interpretation, wispack includes a number of
intuitive plotting methods for visualizing wisp model results.

## Rate-count plots

The backbone of wisp visualization is the *rate-count plot*, as in
“predicted rate, observed count”. Rate-count plots come in two forms,
one giving confidence-intervals (CI) and the other giving a scatter plot
of observations. Both forms give predicted rates across the spatial axis
for all treatment conditions.

### Confidence intervals

The CI rate-count plot is a simplified line plot showing predicted rate
as a function of spatial position, along with a 95% confidence interval
for the rate at each spatial location. As a demonstration, let’s fit a
model of the laminar axis of somatosensory cortex from the RORB
[tutorial](https://michaelbarkasi.github.io/wispack/articles/tutorial_corticallaminar.md),
under the same settings used there, except that the wisp function is
set, as in the following code chunk, to also run a thousand bootstraps
in addition to the default MCMC walk. We want bootstraps because we will
be plotting the model parameters below, and bootstraps generate more
reliable [parameter
estimates](https://michaelbarkasi.github.io/wispack/articles/tutorial_stats.md)
than the MCMC walk. As will be visible in the first rate-count plot
below (at least, when compared to the same plot in the RORB tutorial),
this has a big effect on the CI and parameter estimates. This is the
same fit used in the [stats
tutorial](https://michaelbarkasi.github.io/wispack/articles/tutorial_stats.md);
we’ll also grab the layer-boundary points from the
[CCFv3](https://doi.org/10.1016/j.cell.2020.04.007) registration here,
for use below in the rate-count plots (fitting takes a few minutes,
owing to the bootstraps):

``` r

# Set random seed for reproducibility
set.seed(123)
# Load wispack
library(wispack, quietly = TRUE)
# Grab data
layer.boundary.bins <- read.csv(
  system.file(
      "extdata",
      "S1_layer_boundary_bins.csv",
      package = "wispack"
    )
  )
countdata <- read.csv(
  system.file(
      "extdata",
      "S1_laminar_countdata_demo.csv.gz",
      package = "wispack"
    )
  )
# Fit model
model <- wisp(
    count.data     = countdata,
    variables      = list(species = "gene", ran = "mouse", timeseries = "age"),
    model.settings = list(LRO_cost = "none"),
    plot.settings  = list(
      print.plots  = FALSE,
      dim.bounds   = colMeans(layer.boundary.bins)),
    bootstraps.num = 1000,
    max.fork       = 20,
    verbose        = TRUE
  )
```

``` scroll-output
## 
## Parsing data and settings for wisp model
## ----------------------------------------
## 
## Model settings:
##  buffer_factor: 0.05
##  ctol: 1e-06
##  max_penalty_at_distance_factor: 0.01
##  LRO_cost: none
##  LROcutoff: 2
##  LROwindow_factor: 1.25
##  rise_threshold_factor: 0.8
##  max_evals: 1000
##  rng_seed: 42
##  warp_precision: 1e-07
##  round_none: TRUE
##  trtKO: none
##  max_bin: 0
##  inf_warp: 450359962.73705
## 
## Plot settings:
##  print.plots: FALSE
##  dim.bounds: 78.5, 65.25, 44, 0
##  log.scale: FALSE
##  splitting_factor: 
##  CI_style: TRUE
##  splitting_factor_colors: 120, 240
##  label_size: 5.5
##  title_size: 20
##  axis_size: 12
##  legend_size: 10
##  count_size: 1.5
##  count_jitter: 0.5
##  count.alpha.ran: 0.25
##  count.alpha.none: 0.25
##  pred.alpha.ran: 0.9
##  pred.alpha.none: 1
## 
## MCMC settings:
##  MCMC.burnin: 100
##  MCMC.steps: 1000
##  MCMC.step.size: 0.5
##  MCMC.prior: 1
##  MCMC.neighbor.filter: 1
## 
## Variable dictionary:
##  count: count
##  bin: bin
##  context: context
##  species: gene
##  ran: mouse
##  timeseries: age
##  fixedeffects: hemisphere, age
## 
## Parsed data (head only):
## ------------------------------
##   count bin context species ran hemisphere timeseries
## 1     0 100  cortex  Bcl11b   1       left         12
## 2     0  99  cortex  Bcl11b   1       left         12
## 3     0  93  cortex  Bcl11b   1       left         12
## 4     0  98  cortex  Bcl11b   1       left         12
## 5     0  94  cortex  Bcl11b   1       left         12
## 6     0  94  cortex  Bcl11b   1       left         12
## ----
## 
## Initializing Cpp (wspc) model
## ----------------------------------------
## Context grouping levels: "cortex"
## Species grouping levels: "Bcl11b" "Cux2" "Fezf2" "Nxph3" "Rorb" "Satb2"
## Random-effect grouping levels: "none" "1" "2" "3" "4"
## 
## Infinity handling:
## machine epsilon: (eps_): 2.22045e-16
## pseudo-infinity (inf_): 1e+100
## warp_precision: 1e-07
## implied pseudo-infinity for unbounded warp (inf_warp): 4.5036e+08
## 
## Extracting variables and data information:
## Max bin: 100
## Fixed effects: "hemisphere" "timeseries"
## Ref levels: "left" "12"
## Time series detected: "12" "18"
## Created treatment levels: "ref" "right" "18" "right18"
## Total rows in summed count data table: 12000
## Number of rows with unique model components: 120
## 
## Running LRO change-point detection and setting initial parameters
## ----------------------------------------
## 
## Estimated change points
## Estimated initial parameters
## Number of parameters: 300
## Size of boundary vector: 1260
## 
## Estimating model parameters
## ----------------------------------------
## 
## Running MCMC stimulations
## step: 1/1100
## step: 110/1100
## step: 220/1100
## step: 330/1100
## step: 440/1100
## step: 550/1100
## step: 660/1100
## step: 770/1100
## step: 880/1100
## step: 990/1100
## All complete!
## Acceptance rate (aim for 0.2-0.3): 0.246085
## Penalized nll: 11976.4
## 
## MCMC simulation complete... 
## MCMC run time (total), minutes: 3.779
## MCMC run time (per retained step), seconds: 0.206
## MCMC run time (per step), seconds: 0.206
## 
## Running bootstrap fits
## Performing initial fit of full data
## Penalized nll: 11976.4
## Batch: 1/50, 0.701648 sec/bs
## Batch: 5/50, 0.648595 sec/bs
## Batch: 10/50, 0.628511 sec/bs
## Batch: 15/50, 0.650552 sec/bs
## Batch: 20/50, 0.730031 sec/bs
## Batch: 25/50, 0.663684 sec/bs
## Batch: 30/50, 0.642421 sec/bs
## Batch: 35/50, 0.622834 sec/bs
## Batch: 40/50, 0.60915 sec/bs
## Batch: 45/50, 0.690686 sec/bs
## All complete!
## 
## Bootstrap simulation complete... 
## Bootstrap run time (total), minutes: 11.289
## Bootstrap run time (per sample), seconds: 0.677
## Bootstrap run time (per sample, per thread), seconds: 13.546
## 
## Computing p-values and CI for model parameters
## ----------------------------------------
## 
## Grabbing sample results, only resamples with converged fit... 
## Computing 95% confidence intervals... 
## Estimating p-values from resampled parameters... 
## 
## Recommended resample size for alpha = 0.05, 171 tests
## with bootstrapping/MCMC: 3420
## Actual resample size: 1001
## 
## Stat summary (head only):
## ------------------------------
##                                  parameter      estimate      CI.low    CI.high   p.value p.value.adj    alpha.adj significance
## 1       baseline_cortex_Rt_Bcl11b_Tns/Blk1  3.3880114451  3.30500897 3.40212815        NA          NA           NA             
## 2   beta_Rt_cortex_Bcl11b_right_X_Tns/Blk1 -0.0003773747 -0.07507784 0.07339858 0.9890110    2.967033 0.0166666667           ns
## 3      beta_Rt_cortex_Bcl11b_18_X_Tns/Blk1  0.6072410701  0.51496657 0.72377372 0.0000000    0.000000 0.0002923977          ***
## 4 beta_Rt_cortex_Bcl11b_right18_X_Tns/Blk1  0.0461743004 -0.08550940 0.14890972 0.2657343   13.818182 0.0009615385           ns
## 5       baseline_cortex_Rt_Bcl11b_Tns/Blk2  2.2781929967  1.99837574 2.55129060        NA          NA           NA             
## 6   beta_Rt_cortex_Bcl11b_right_X_Tns/Blk2  0.0822301676 -0.45766257 0.51713430 0.6013986   16.237762 0.0018518519           ns
## ----
## 
## Computing 95% CI for predicated values by bin... 
## 
## Analyzing residuals
## ----------------------------------------
## 
## Computing residuals... 
## Making masks... 
## Making plots and saving stats... 
## 
## Log-residual summary by grouping variables (head only):
## ------------------------------
##          group      mean        sd  variance
## 1          all 0.2389661 0.3918609 0.1535550
## 2 ran_lvl_none 0.2321238 0.3174750 0.1007903
## 3    ran_lvl_1 0.2384415 0.3968682 0.1575043
## 4    ran_lvl_2 0.2132537 0.3938506 0.1551183
## 5    ran_lvl_3 0.2751780 0.4605209 0.2120795
## 6    ran_lvl_4 0.2426757 0.4399106 0.1935213
## ----
## 
## Making plots
## ----------------------------------------
## 
## Making MCMC walks plots... 
## Making MCMC vs bootstrap parameter distribution comparison plots... 
## Making effect parameter distribution plots... 
## Making rate-count plots... 
## Making time series plots... 
## Making parameter plots...
```

Let’s use the plot.ratecount function to remake the RORB plot that was
made and shown in the RORB
[tutorial](https://michaelbarkasi.github.io/wispack/articles/tutorial_corticallaminar.md),
this time using the bootstrap resamples for the CI instead of MCMC
resamples:

``` r

plot.ratecount(
    model,
    dim.boundaries = colMeans(layer.boundary.bins),
    species        = "Rorb",
    verbose        = FALSE
  )
```

![](tutorial_wispplots_files/figure-html/RORB_CI-1.png)

CI rate-count plots are straightforward to interpret: The horizontal (x)
axis shows spatial position along the modeled spatial dimension, the
vertical (y) axis shows predicated rate. The solid line represents
predicated rate, after random effects have been filtered out. The shaded
ribbons show the 95% confidence interval around the predicated rate. As
explained in the RORB
[tutorial](https://michaelbarkasi.github.io/wispack/articles/tutorial_corticallaminar.md),
this confidence interval is computed without any multiple comparison
correction (unlike p-values, which are corrected for multiple
comparisons). As noted above, the different ribbon sizes in this plot
compared to the plot from the RORB
[tutorial](https://michaelbarkasi.github.io/wispack/articles/tutorial_corticallaminar.md)
are due to using bootstraps instead of MCMC samples.

Both kinds of rate-count plots separate treatment conditions. In this
data, we have two covariates, hemisphere and age, each with two levels
(left vs right for hemisphere, P12 vs P18 for age). Thus, the plot shows
the predications for the four possible combinations of covariates (i.e.,
the reference condition of left-12, plus three treatment conditions).

By default, both kinds of rate-count plots show rate as a raw value.
However, it’s possible to rescale the vertical axis into log form:

``` r

plot.ratecount(
    model,
    dim.boundaries = colMeans(layer.boundary.bins),
    species        = "Rorb",
    log.scale      = TRUE,
    verbose        = FALSE
  )
```

![](tutorial_wispplots_files/figure-html/RORB_CI_log-1.png)

### Scatter plots

The second kind of rate-count plot, the *scatter rate-count plot*, is a
more complex plot showing predicted rate without random effects as a
line, model fit to individual replicates (i.e., predicted rates with
random effects) as dashed lines, and observed counts (plus
interpolations for random-effect-free counts) as scatter points. It’s
accessed via the CI_style argument in plot.ratecount:

``` r

plot.ratecount(
    model,
    dim.boundaries = colMeans(layer.boundary.bins),
    CI_style       = FALSE,
    species        = "Rorb",
    verbose        = FALSE
  )
```

![](tutorial_wispplots_files/figure-html/RORB_ran-1.png)

As with CI rate-count plots, the horizontal axis represents spatial
position along the modeled dimension, and the vertical axis represents
predicted rate.

The main advantages of scatter rate-count plots are that they show how
well the model fits the observed data, and also show the variability in
individual replicates. Still, they can be visually crowded, so wispack
provides a function to decompose them into their components.

``` r

plots.decomp <- plot.decomposition(
    model,
    dim.boundaries = colMeans(layer.boundary.bins),
    species        = "Rorb"
  )
```

First, let’s look at separate plots of counts (observed vs extrapolated
without random effets) and predicted rates (fit of observed counts vs
predicted rates without random effects):

``` r

grid.arrange(grobs = plots.decomp[2:5], ncol = 2)
```

![](tutorial_wispplots_files/figure-html/RORB_ran_decomposed_print1-1.png)

The upper left plot (observed counts) shows one point for the total
count of RORB transcripts observed in a given replicate at a given
spatial bin. The upper right plot (extrapolated counts) shows
extrapolated counts got by taking the mean of the actual replicates
within a bin. These means stand in for observations from a hypothetical
replicate with no noise or other random effects. The lower left plot
(observed predictions) shows the model fit to the observed counts. The
lower right plot (extrapolated predictions) shows the fit to the
extrapolated counts, which stand in for the predicted rate without
random effects.

Next, let’s look at observed counts (points) and their fits (dashed
lines) on the same plot, but broken out by sample:

``` r

grid.arrange(grobs = plots.decomp[7:length(plots.decomp)], ncol = 2)
```

![](tutorial_wispplots_files/figure-html/RORB_ran_decomposed_print2-1.png)

Notice how each plot has two dashed lines which fit the observed counts
well, and two which do not. This is because the model is fit to each
sample in all treatment conditions. The treatment conditions are got by
the interaction of covariates hemisphere and age. Although each sample
includes both hemispheres, it can only include one age. The model pools
information across samples to extrapolate the fit for the missing
treatment conditions. Thus, some of the dashed lines represent fits to
observed counts in treatment conditions that were not actually observed
in that sample.

## Time-series plots

For models fit with a covariate specified as a [time
series](https://michaelbarkasi.github.io/wispack/articles/tutorial_timeseries.md),
wispack provides a function to plot predicted rates as a function of
time. This is done with the plot.timeseries function. Although there are
only two ages in the cortical laminar data, we can still use this
function to visualize how predicted rates change with age. Let’s look at
the time-series plot for RORB again:

``` r

plot.timeseries(
    model, 
    species = "Rorb",
    verbose = FALSE
  )
```

![](tutorial_wispplots_files/figure-html/RORB_timeseries-1.png)

As explained in the [time series
tutorial](https://michaelbarkasi.github.io/wispack/articles/tutorial_timeseries.md),
the rate represented by the height of the solid lines is the rate
parameter itself (for each block), not the actual predicted rate. The
actual predicated rate at a given bin (or the mean predicated rate
across an entire block) can be lower than the rate parameter, depending
on the slope parameters of the transitions. In time series plots, the
points represent the mean observed rate across an entire block (with
random effects filtered out). Thus, as is the case in this example, the
distance between the observed points and the predicted rates is not
necessarily a sign of a bad fit. As can be seen in the rate-count plots
above, the fit is quite good.

## Parameter plots

Rate-count and time-series plots show observed and predicted gene
expression rates as a function of spatial position, time, and (at least
implicitly) treatment conditions. However, with the exception of
time-series plots showing rate parameters, these plots do not directly
show the underlying rate, transition point, and slope parameters – or
the corresponding random-effect scalars – that define a wisp model.

``` r

plots.param <- plot.parameters(
    model,
    species = "Rorb",
    verbose = FALSE
  )
```

The plot.parameters function produces separate plots for the baseline
(reference condition) parameters, the treatment condition parameters,
and the random-effect scalars. Let’s look at the baseline parameters
first:

``` r

plots.param$plot_baseline_cortex_Rorb
```

![](tutorial_wispplots_files/figure-html/RORB_parameters_baseline-1.png)

There are three separate plots: Rorb.rate gives the baseline (left, P12)
block rates in log space, Rorb.tpoint gives the baseline transition
point locations, and Rorb.tslope gives the baseline transition slopes
for each of the three transition points. The horizontal axis tick labels
are formatted as \<param type\>\_T/B\<num\>, where \<param type\> is Rt
(rate), tpoint (transition point location), or tslope (transition
slope), T/B stands for “transition/block”, and \<num\> gives the
transition or block number. For example, Rt_T/B2 is the rate parameter
for block 2, while tpoint_T/B1 is the location of the first transition
point, i.e., the transition between block 1 and block 2. For each
horizontal axis tick, there is a violin plot showing the distribution of
the estimates for that parameter across the bootstraps (or the MCMC
walk), with a longer horizontal line showing the best-fit estimate and
two shorter horizontal bars showing the 95% confidence interval.

There are three true (i.e., non-baseline) treatment conditions: the left
hemisphere at P18, the right hemisphere at P12, and the right hemisphere
at P18. Let’s look first at the parameters giving the effect of
switching hemispheres in the reference condition, i.e., the shift from
left to right hemisphere at P12:

``` r

plots.param$plot_treatment_right_cortex_Rorb
```

![](tutorial_wispplots_files/figure-html/RORB_parameters_treatment1-1.png)

Notice that interpreting the parameters in the baseline condition is
straightforward. For example, Rt_T/B2 is the actual rate in block 2.
However, the parameters for the treatment conditions are additive
offsets, i.e., log-fold changes. Thus, positive parameter values
represent an *increase* in the value, while negative parameter values
represent a *decrease* in the value.

Another point of difference between the baseline and (true) treatment
condition parameters is that the latter are tested for significance
against a null hypothesis of zero (no change from baseline). The results
of these tests are shown above the parameter in the typical notation:
“ns” is not significant, “\*” is p \< 0.05, “\*\*” is p \< 0.01, and
“\*\*\*” is p \< 0.001.

``` r

plots.param$plot_treatment_18_cortex_Rorb
```

![](tutorial_wispplots_files/figure-html/RORB_parameters_treatment2-1.png)

The above plot shows the parameters for the effect of age (the shift
from P12 to P18) in the reference condition (left hemisphere). The next
plot shows the interaction effect of age in the right hemisphere:

``` r

plots.param$plot_treatment_right18_cortex_Rorb
```

![](tutorial_wispplots_files/figure-html/RORB_parameters_treatment3-1.png)

Finally, we can look at the random-effect scalars for each parameter:

``` r

plots.param$plot_ranEff_cortex_Rorb
```

![](tutorial_wispplots_files/figure-html/RORB_parameters_random-1.png)

There is one scalar for each parameter type (rate, transition point,
transition slope), with each sample (i.e., replicate) having its own set
of scalars. For example, in the above plot we can see at a glance that
samples 1 and 2 have rates that are higher than expected while samples 3
and 4 have rates that are lower than expected. The exact interpretation
of the random-effect parameters is discussed in the [warping
tutorial](https://michaelbarkasi.github.io/wispack/articles/tutorial_warping.md).
