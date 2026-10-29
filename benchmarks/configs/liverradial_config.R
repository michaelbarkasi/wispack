# Regression-check config for mclapply_vs_fork_general.R: reproduces the exact dataset/settings
# used by the original, hardcoded mclapply_vs_fork_isolated.R (see benchmarks/RESULTS.md, Round 2),
# so that running the generalized script against this config should give numbers close to that
# existing table -- confirming generalizing the script didn't change the measured behavior.
# 276 fitted parameters, same liverradial tutorial dataset.

label <- "liverradial"

count.data <- read.csv(system.file("extdata", "Droin_radial_count_data_sim.csv.gz", package = "wispack"))

data.variables <- list(
  count      = "count", bin = "bin", context = "liver", species = "gene",
  ran        = "mouse", timeseries = "ZT"
)

model.settings <- list(LRO_cost = "none", LROcutoff = 1.5)

# Same grid as mclapply_vs_fork_isolated.R's Round 2 (minus a couple of the smaller bs_num_max=100
# repeats, to keep the regression check itself reasonably fast -- bs_num_max=1000,N=20 is the point
# that most clearly separates the two arms in the original result, and a couple of bs_num_max=100
# rows are enough to check the smaller-gap end of the table too).
grid <- data.frame(
  bs_num_max = c(1000, 100, 100),
  N          = c(20,   5,   20)
)

reps <- 3
