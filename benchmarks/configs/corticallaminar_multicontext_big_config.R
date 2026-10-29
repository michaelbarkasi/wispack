# "Big model" config for mclapply_vs_fork_general.R, built from this package's own tutorial data
# (no new data collection -- real inst/extdata rows only, no fabricated values).
#
# How this was built, following the task's suggested order of approaches:
#
#  1. Finer spatial binning (tutorial_corticallaminar.Rmd's suggestion): NOT feasible with the
#     shipped extdata. ACx_laminar_countdata_demo*.csv.gz and S1_laminar_countdata_demo.csv.gz
#     only ship the already-binned `bin` column (100 levels); the original continuous per-cell
#     laminar coordinates aren't included, so there's nothing to re-bin more finely without
#     fabricating coordinates. Lowering LROcutoff (more sensitive change-point detection) on the
#     single-context (cortex-only) cortical-laminar data was tried but, at fixed real data, this
#     approach tops out in the hundreds-to-low-thousands of parameters -- much less than approach 2.
#
#  2. Multi-context scaling (tutorial_multicontext.Rmd's suggestion): this is what's used here.
#     tutorial_multicontext.Rmd itself explains that fixed-effect parameter count scales linearly
#     with the number of `context` (cell-type) levels, and collapses the real 31 MapMyCells cell
#     types in ACx_laminar_countdata_demo_with_celltype.csv.gz down to just 2 (Glut/GABA) "to make
#     results easier to visualize and speed up model fitting" -- i.e. that collapsing is a
#     demonstration simplification, not a modeling necessity. This config instead keeps the
#     cell-type labels mostly as-is (undoing just that tutorial's binary collapse) and uses the
#     12 best-populated of the real 31 cell types (>=1278 rows each; the full breakdown is in
#     benchmarks/RESULTS_large_model.md), dropping only the extremely sparse tail (as few as 12
#     total rows for some of the 31 types) where per-bootstrap-resample fits became numerically
#     unstable even for a single fit. This reaches 4,788 fitted parameters (17x the original
#     276-parameter liverradial benchmark) at 145,776 data rows, confirmed via the same fit()-based
#     parameter count probe mclapply_vs_fork_general.R uses for every config.
#
#  3. Honesty check: 4,788 parameters clears the task's "1000+ parameters" bar comfortably, so a
#     full escalation to multicontext's independently-suggested "multiple contexts x multiple
#     genes/species x ~100 bins" combination (e.g. combining ACx and S1 cortical data, or adding
#     more cell types) was not pursued further. Scaling further is possible (confirmed up to
#     18 cell types / 6,786 parameters in exploration, see RESULTS_large_model.md), but single
#     per-bootstrap fit time grows with it too (roughly 30s/fit at 12 cell types vs. 0.15-0.3s/fit
#     for the 276-parameter model), so 12 cell types was chosen as the point that is unambiguously
#     "much bigger" while keeping the actual benchmark grid below (several grid rows x reps x 2
#     arms, each potentially minutes long) tractable to actually run.
#
# Occasional "Could not find nearby feasible parameters (boundary violation or nan rates)" warnings
# during individual bootstrap resample fits are expected and harmless here -- they're the model's
# own graceful fallback (keep the previous parameter estimate) for a resample that happens to hit a
# numerically awkward spot for one of the many context/gene/treatment combinations in a model this
# large, not a crash or a hang. This benchmark is about fork_utils.h vs. mclapply mechanics and
# memory, not about the statistical validity of any individual bootstrap estimate.

label <- "corticallaminar_multicontext_big"

.countdata_full <- read.csv(system.file(
  "extdata", "ACx_laminar_countdata_demo_with_celltype.csv.gz", package = "wispack"
))
.tab <- sort(table(.countdata_full$celltype_MMC), decreasing = TRUE)
.keep_types <- names(.tab)[seq_len(12)] # top-12 best-populated of the real 31 cell types

count.data <- .countdata_full[.countdata_full$celltype_MMC %in% .keep_types, ]

data.variables <- list(
  context    = "celltype_MMC",
  species    = "gene",
  ran        = "mouse",
  timeseries = "age"
)

# Same LRO settings tutorial_multicontext.Rmd itself uses (skip the AIC/BIC/NLL grid search,
# use a widened window factor) -- this benchmark isn't about re-tuning the change-point search,
# just about reaching a large, realistic parameter count with real tutorial settings.
model.settings <- list(LRO_cost = "none", LROwindow_factor = 1.75)

# Each individual bootstrap fit costs ~30s single-threaded at this model size (vs. ~0.15-0.3s for
# the 276-parameter liverradial model) -- bs_num_max is kept much smaller than the original
# benchmark's 100/1000 so the grid finishes in a tractable amount of wall-clock time. The grid
# still spans the same qualitative comparison: fixed bs_num_max across increasing N (fork-count
# amortization), and a doubled bs_num_max at fixed N (total-fork-count effect, mirroring the
# original benchmark's 100-vs-1000-at-N=20 comparison).
grid <- data.frame(
  bs_num_max = c(20, 20, 20, 40),
  N          = c(4,  10, 20, 20)
)

reps <- 2
