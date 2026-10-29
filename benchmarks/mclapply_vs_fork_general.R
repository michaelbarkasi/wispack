# Generalized version of mclapply_vs_fork_isolated.R: same isolated-bootstrap-step methodology
# (cpp_model$bs_batch() vs. cpp_model$fit() + mclapply(cpp_model$bs_fit(...))), but driven by a
# config file instead of a hardcoded dataset/model, and with added best-effort peak-combined-RSS
# tracking (see mem_watch_utils.R) alongside wall-clock time for every timed call.
#
# mclapply_vs_fork_isolated.R itself is left untouched as the historical record of the original
# 276-parameter liverradial result (see benchmarks/RESULTS.md, Round 2). This script generalizes
# that methodology to any wispack dataset/model configuration so it can be pointed at much larger
# models.
#
# Usage:
#   Rscript benchmarks/mclapply_vs_fork_general.R <path/to/config.R>
# or, if CONFIG_PATH is already set as a variable before sourcing this file (e.g. from another
# script), that value is used instead and commandArgs() is not consulted.
#
# Config file contract (a plain R script, sourced with source(..., local = FALSE) into the global
# environment -- see benchmarks/configs/*.R for examples):
#   count.data      -- a data.frame, already loaded (e.g. via read.csv on an inst/extdata file).
#   data.variables  -- the wisp() `variables` list (column-name mapping).
#   model.settings  -- the wisp() `model.settings` list.
#   grid            -- a data.frame with integer columns bs_num_max and N (N = max.fork / mc.cores),
#                       one row per (bs_num_max, N) configuration to benchmark.
#   reps            -- integer, number of repetitions per grid row.
#   label           -- character, short identifier used for default output filenames and in
#                       printed progress messages.
# Optional:
#   plot.settings   -- passed to wisp()'s plot.settings during cpp_model construction (only
#                       affects how the model is initialized; no plots are ever produced by this
#                       benchmark). Defaults to list(print.plots = FALSE).
#   output_csv      -- path to write the raw per-rep results to. Defaults to
#                       "benchmarks/mclapply_vs_fork_general_<label>_results.csv".
#   mem_poll_interval -- seconds between RSS polls (default 0.3).

devtools::load_all(".", quiet = TRUE)
suppressMessages(library(parallel))

# Always run this script (and point CONFIG_PATH / output_csv) from the package root -- that's the
# one documented usage, so a plain relative path works for sourcing the memory-watcher helpers too.
source("benchmarks/mem_watch_utils.R")

# --- Resolve and load the config -------------------------------------------------------------

config_path <- if (exists("CONFIG_PATH", inherits = FALSE)) {
  CONFIG_PATH
} else {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) {
    stop("Usage: Rscript benchmarks/mclapply_vs_fork_general.R <path/to/config.R>")
  }
  args[[1]]
}
cat("Loading config:", config_path, "\n")
source(config_path, local = FALSE)

required_vars <- c("count.data", "data.variables", "model.settings", "grid", "reps", "label")
missing_vars <- setdiff(required_vars, ls(envir = .GlobalEnv))
if (length(missing_vars) > 0) {
  stop("Config file is missing required variable(s): ", paste(missing_vars, collapse = ", "))
}
if (!exists("plot.settings", inherits = FALSE)) plot.settings <- list(print.plots = FALSE)
if (!exists("mem_poll_interval", inherits = FALSE)) mem_poll_interval <- 0.3
if (!exists("output_csv", inherits = FALSE)) {
  output_csv <- file.path("benchmarks", paste0("mclapply_vs_fork_general_", label, "_results.csv"))
}
stopifnot(all(c("bs_num_max", "N") %in% colnames(grid)))

# --- Shared helper: capture a freshly LRO-initialized cpp_model from inside wisp() -----------
# (identical in spirit to mclapply_vs_fork_isolated.R's capture_cpp_model(): traces wisp() right
# after its call to LRO_initial_param_ests(), grabs the cpp_model local variable at that point, and
# aborts the rest of wisp() via a sentinel error so no fitting/plotting happens inside wisp() itself.)

capture_cpp_model <- function(...) {
  ns <- asNamespace("wispack")
  wisp_fun <- get("wisp", envir = ns)
  body_exprs <- as.list(body(wisp_fun))
  idx <- which(sapply(
    body_exprs,
    function(e) grepl("LRO_initial_param_ests\\(verbose, 0, 0\\)", paste(deparse(e), collapse = " "))
  ))
  fire_idx <- idx[1] + 1
  trace(
    "wisp", where = ns, print = FALSE, at = fire_idx,
    tracer = quote({
      assign("captured_cpp_model", cpp_model, envir = .GlobalEnv)
      stop("__CAPTURED__")
    })
  )
  tryCatch(
    do.call(get("wisp", envir = ns), list(...)),
    error = function(e) if (!grepl("__CAPTURED__", conditionMessage(e))) stop(e)
  )
  untrace("wisp", where = ns)
  captured <- get("captured_cpp_model", envir = .GlobalEnv)
  rm("captured_cpp_model", envir = .GlobalEnv)
  captured
}

make_cpp_model <- function() {
  capture_cpp_model(
    count.data    = count.data, variables = data.variables, model.settings = model.settings,
    plot.settings = plot.settings, verbose = FALSE
  )
}

# --- Probe model size once, up front, for reporting -------------------------------------------

cat("Probing model size (fitting once, outside the timed benchmark)...\n")
probe_model <- make_cpp_model()
probe_fit <- time_and_track_mem(function() probe_model$fit(FALSE), poll_interval = mem_poll_interval)
probe_res <- probe_model$results()
n_params <- length(probe_res$fitted.parameters)
cat(sprintf(
  "Model '%s': %d rows, %d fitted parameters, single full fit() took %.1fs (peak RSS %.0f MB, single process, no parallelism).\n",
  label, nrow(count.data), n_params, probe_fit$sec, probe_fit$peak_rss_mb
))
single_fit_peak_rss_mb <- probe_fit$peak_rss_mb
rm(probe_model, probe_res)
gc(verbose = FALSE)

# --- Benchmark loop -----------------------------------------------------------------------------

results <- data.frame()

for (gi in seq_len(nrow(grid))) {
  bs_num_max <- grid$bs_num_max[gi]
  N          <- grid$N[gi]
  for (rep in seq_len(reps)) {

    cat(sprintf("[%s] bs_num_max=%d N=%d rep=%d/%d: bs_batch()... ", label, bs_num_max, N, rep, reps))
    cpp_model_a <- make_cpp_model()
    r_fork <- time_and_track_mem(
      function() cpp_model_a$bs_batch(bs_num_max, N, FALSE, FALSE),
      poll_interval = mem_poll_interval
    )
    rm(cpp_model_a); gc(verbose = FALSE)
    cat(sprintf("%.1fs (peak RSS %.0f MB). mclapply... ", r_fork$sec, r_fork$peak_rss_mb))

    cpp_model_b <- make_cpp_model()
    r_mcl <- time_and_track_mem(
      function() {
        cpp_model_b$fit(FALSE)
        rows <- parallel::mclapply(
          seq_len(bs_num_max) - 1L,
          function(i) cpp_model_b$bs_fit(i, TRUE),
          mc.cores = N, mc.preschedule = TRUE
        )
        do.call(rbind, rows)
      },
      poll_interval = mem_poll_interval
    )
    rm(cpp_model_b); gc(verbose = FALSE)
    cat(sprintf("%.1fs (peak RSS %.0f MB).\n", r_mcl$sec, r_mcl$peak_rss_mb))

    results <- rbind(results, data.frame(
      label = label, bs_num_max = bs_num_max, N = N, rep = rep,
      bs_batch_sec = r_fork$sec, bs_batch_peak_rss_mb = r_fork$peak_rss_mb,
      mclapply_sec = r_mcl$sec, mclapply_peak_rss_mb = r_mcl$peak_rss_mb
    ))
  }
}

write.csv(results, output_csv, row.names = FALSE)
cat("\nWrote raw results to", output_csv, "\n")

cat(sprintf("\n=== %s: isolated bootstrap-step comparison (mean over %d reps) ===\n", label, reps))
agg <- aggregate(
  cbind(bs_batch_sec, bs_batch_peak_rss_mb, mclapply_sec, mclapply_peak_rss_mb) ~ bs_num_max + N,
  data = results, FUN = mean
)
agg$speedup_mclapply_over_fork    <- agg$bs_batch_sec / agg$mclapply_sec
agg$mclapply_rss_over_fork_rss    <- agg$mclapply_peak_rss_mb / agg$bs_batch_peak_rss_mb
agg <- agg[order(agg$bs_num_max, agg$N), ]
print(agg, row.names = FALSE)

cat(sprintf(
  "\nModel size: %d fitted parameters, %d data rows. Single-process, non-parallel full fit() peak RSS: %.0f MB.\n",
  n_params, nrow(count.data), single_fit_peak_rss_mb
))
