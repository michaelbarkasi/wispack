
# Set random seed for reproducibility
# ... R only. C++ seed set in its code
ran.seed <- 123
set.seed(ran.seed)

# Load wispack
library(wispack)

# Load demo data
# ... same rows as S1_laminar_countdata_demo.csv (used in demo/full_options.R), just renamed to
# ... wisp()'s default variable names (species, ran) so no data.variables list is needed below
countdata <- read.csv(
  system.file(
      "extdata",
      "S1_laminar_countdata_demo.csv.gz",
      package = "wispack"
      )
  )
rename.map <- c(gene = "species", mouse = "ran")
names(countdata)[match(names(rename.map), names(countdata))] <- rename.map

# Fit model
laminar.model <- wisp(countdata)

# View model 
View(laminar.model)