
# Set for debugging (if needed): Sys.setenv(CXXFLAGS="-fsanitize=address -g -O1")
# ... in external terminal: cd to working directory, run "scratch.R"

# Clear project
rm(list = ls())

# Set random seed for reproducibility
# ... R only. C++ seed set in its code
ran.seed <- 123
set.seed(ran.seed)

# Load wispack
library(wispack)

# Set bootstrap chunk size
sys_name <- Sys.info()["sysname"]
if (sys_name == "Darwin" || sys_name == "Linux") {
  bs_chunksize <- 10
} else {
  bs_chunksize <- 0
}

# Load demo data
# ... from horizontal (axial) slice of mouse cortex, primary somatosensory cortex (L1 removed)
data_path <- system.file("extdata", "S1_laminar_countdata_demo.csv", package = "wispack")
countdata <- read.csv(data_path)
# ... Manual load: countdata <- read.csv("S1_laminar_countdata_demo.csv")
# ... Note on spatial coordinates: 
#  y-axes: laminar, 0 is bottom of L6b, 100 is top of L2/3
#  x-axes: columnar, 0 is most posterior, 100 is most anterior

# Grab cortical layer boundaries, gotten from CCFv3 registration
# ... Not used to fit model; used only for comparison
# ... Numbers represent the lower boundary of each layer, in bins
# ... Higher numbers closer to cortical surface
boundary_path <- system.file("extdata", "layer_boundary_bins.csv", package = "wispack")
layer.boundary.bins <- read.csv(boundary_path)
# ... Manual load: layer.boundary.bins <- read.csv("layer_boundary_bins.csv")

# Define fixed effects to test
fixed.effect.names <- c("hemisphere", "age")

# Define variables in the dataframe for the model
data.variables <- list(
    count = "count",
    bin = "bin", 
    parent = "cortex", 
    child = "gene",
    ran = "mouse",
    fixedeffects = fixed.effect.names
  )

# Model settings 
# ... all settings shown here are defaults
model.settings <- list(
    # ... these are global options needed to set up model
    buffer_factor = 0.05,                                 # buffer factor for penalizing distance from structural parameter values
    ctol = 1e-6,                                          # convergence tolerance
    max_penalty_at_distance_factor = 0.01,                # maximum penalty at distance from structural parameter values
    LROcutoff = 2.0,                                      # cutoff for LROcp, a multiple of standard deviation
    LROwindow_factor = 1.25,                              # window factor for LROcp, larger means larger rolling window
    rise_threshold_factor = 0.8,                          # amount of detected rise as fraction of total required to end run in initial slope estimation
    max_evals = 1000,                                     # maximum number of evaluations for optimization
    rng_seed = 42,                                        # random seed for optimization (controls bootstrap resamples only)
    warp_precision = 1e-7                                 # decimal precision to retain when selecting really big number as pseudo infinity for unbound warping
  )

# Setting suggestions: 
# - Recommend turning settings on data for which you have strong priors. 
# - If bootstraps not fitting due to falling off boundary (very high boundary penalty and few iterations), 
#    increase max_penalty_at_distance_factor from 0.01 to 0.1 or 1.0, so the gradient descent algorithm has 
#    more information about the boundary edge. 
# - Adjusting the LROcutoff is another way of controlling how the model finds rate transitions. Higher values 
#    mean fewer detected transitions. 

# Settings for MCMC walk
# ... all settings shown here are defaults
MCMC.settings <- list(
    MCMC.burnin = 0,
    MCMC.steps = 1e3,
    MCMC.step.size = 1.0,
    MCMC.prior = 1.0, 
    MCMC.neighbor.filter = 2
  )

# Fit model
# ... all settings shown here are defaults
laminar.model <- wisp(
    # Data to model
    count.data = countdata,
    # Variable labels
    variables = data.variables,
    # Settings used on R side
    use.median = FALSE,
    MCMC.settings = MCMC.settings,
    bootstraps.num = 0,
    converged.resamples.only = TRUE,
    max.fork = bs_chunksize,
    dim.bounds = colMeans(layer.boundary.bins),
    verbose = TRUE,
    print.child.summaries = TRUE,
    # Setting to pass to C++ model
    model.settings = model.settings
  )

# Demo plots showing anatomy of a wisp
demo.sigmoid.plots()
demo.warp.plots() 