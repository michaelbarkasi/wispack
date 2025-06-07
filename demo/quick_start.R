
# Set random seed for reproducibility
# ... R only. C++ seed set in its code
ran.seed <- 123
set.seed(ran.seed)

# Load wispack
library(wispack)

# Load demo data
countdata <- read.csv(system.file("extdata", "S1_laminar_countdata_demo_default_col_names.csv", package = "wispack"))

# Fit model
laminar.model <- wisp(countdata)

# View model 
View(laminar.model)