# Benchmark: parallel::mclapply() vs. the current fork_utils.h-based wspc::bs_batch() for bootstrap-resample
# fitting, on THIS codebase (not theory).
#
# Context: bs_batch() (src/wspc.cpp) forks a fresh batch of up to `max.fork` child processes for every batch
# of `max.fork` bootstraps (see src/fork_utils.h's run_forked_batches()), i.e. it forks once per bootstrap in
# total, synchronizing (pipe write/read + waitpid) at every batch boundary. parallel::mclapply()'s default
# mc.preschedule = TRUE instead forks exactly `mc.cores` long-lived workers up front, each of which processes
# its share of the task list in a sequential loop before reporting back once. A prior stress test (see this
# project's side-session notes) already confirmed that repeated bs_fit(idx, clear_stan = TRUE) calls within one
# persistent process are memory-safe over 20,000 sequential bootstraps -- the precondition for trying this.
#
# This script is a benchmarking exercise only. It does NOT replace bs_batch()'s real implementation, and the
# .method("bs_fit", ...) exposure added to the RCPP_MODULE in src/wspc.cpp exists purely to make this possible
# from R. Nothing here is wired into wisp()'s actual fitting path. See the commit that added this file for
# more detail; this is intentionally kept out of R CMD build/check (see .Rbuildignore).
#
# Usage: Rscript benchmarks/mclapply_vs_fork_bs_batch.R

devtools::load_all(".", quiet = TRUE)
suppressMessages(library(parallel))

# --- Capture a fully LRO-initialized cpp_model exactly as bs_batch() would receive it ----------------------
#
# wisp() builds and fits its C++ model (`cpp_model`) as a local variable; it's never returned to R. To drive
# wspc::bs_fit() directly (for the mclapply prototype) using the *exact* same preprocessing wisp() itself does
# (column relabeling, fixed-effect parsing, LRO change-point search, initial parameter estimates -- all of
# which are nontrivial to reimplement faithfully by hand), this installs a one-shot trace() on wisp() that
# grabs `cpp_model` right after cpp_model$LRO_initial_param_ests(...) -- the same point bs_batch()'s own
# fit(false) call starts from -- and aborts the wisp() call immediately via a sentinel condition, so no time is
# wasted on the fit/MCMC/bootstrap work wisp() would otherwise go on to do.
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

# --- The mclapply-based prototype -----------------------------------------------------------------------
#
# Mirrors bs_batch()'s own shape: fit the full data once (the same fit() call bs_batch() makes for its last
# result row), then run all bootstrap fits via mclapply instead of run_forked_batches(). Returns a matrix of
# just the bootstrap rows (params + bnll, nll, success.code, num.evals per row), the same column layout
# bs_fit() returns per call; the full-data fit row is reported separately since there's no R-exposed way to
# pull the raw post-fit() params/diagnostics vector the way bs_batch() does internally in C++ (fitted.parameters
# via results() only gives the params, not bnll/nll/success/evals) -- immaterial to the timing comparison.
run_mclapply_bs <- function(cpp_model, bs_num_max, mc.cores, verbose = FALSE) {
  cpp_model$fit(verbose)
  rows <- parallel::mclapply(
    seq_len(bs_num_max) - 1L,
    function(i) cpp_model$bs_fit(i, TRUE),
    mc.cores = mc.cores,
    mc.preschedule = TRUE
  )
  bad <- vapply(rows, function(r) inherits(r, "try-error") || is.null(r), logical(1))
  if (any(bad)) stop(sprintf("%d/%d mclapply bootstrap fits failed", sum(bad), bs_num_max))
  do.call(rbind, rows)
}

# --- Benchmark harness -----------------------------------------------------------------------------------

time_call <- function(expr_fun) {
  t0 <- Sys.time()
  expr_fun()
  as.numeric(difftime(Sys.time(), t0, units = "secs"))
}

count_data <- read.csv(system.file("extdata", "Droin_radial_count_data_sim.csv.gz", package = "wispack"))
data.variables <- list(count = "count", bin = "bin", context = "liver", species = "gene", ran = "mouse", timeseries = "ZT")
model.settings <- list(LRO_cost = "none", LROcutoff = 1.5)

cat("Capturing a reference cpp_model (LRO + initial param ests only, no fit yet)...\n")
cpp_model_template <- capture_cpp_model(
  count.data     = count_data,
  variables      = data.variables,
  model.settings = model.settings,
  plot.settings  = list(print.plots = FALSE),
  verbose        = FALSE
)
cat("Done. n_params =", length(cpp_model_template$results()$fitted.parameters), "\n\n")

results <- data.frame()

run_phase <- function(bs_num_max, fork_values, reps, label) {
  for (N in fork_values) {
    for (rep in seq_len(reps)) {
      cat(sprintf("[%s] bs_num_max=%d N=%d rep=%d/%d: real wisp()... ", label, bs_num_max, N, rep, reps))
      t_real <- time_call(function() {
        wisp(
          count.data     = count_data,
          variables      = data.variables,
          bootstraps.num = bs_num_max,
          max.fork       = N,
          model.settings = model.settings,
          MCMC.settings  = list(MCMC.steps = 0),
          plot.settings  = list(print.plots = FALSE),
          verbose        = FALSE
        )
      })
      cat(sprintf("%.1fs. mclapply prototype... ", t_real))

      # Fresh cpp_model per timed call: bs_fit() mutates internal fit state, so reuse would let one method's
      # run warm-start off the other's results.
      cpp_model <- capture_cpp_model(
        count.data = count_data, variables = data.variables, model.settings = model.settings,
        plot.settings = list(print.plots = FALSE), verbose = FALSE
      )
      t_mcl <- time_call(function() run_mclapply_bs(cpp_model, bs_num_max, mc.cores = N))
      cat(sprintf("%.1fs.\n", t_mcl))

      results <<- rbind(results, data.frame(
        label = label, bs_num_max = bs_num_max, N = N, rep = rep,
        real_fork_sec = t_real, mclapply_sec = t_mcl
      ))
    }
  }
}

# Phase A: apples-to-apples at the real liverradial-test scale (1000 bootstraps, max.fork = 20), matching
# tests/testthat/test-slow-liverradial.R exactly.
run_phase(bs_num_max = 1000, fork_values = 20, reps = 3, label = "A_full_scale")

# Phase B: does the gap (if any) scale with N, or stay constant? Smaller bs_num_max so this finishes in a
# reasonable time; N = 1 included once at an even smaller scale since at N = 1 neither method forks at all
# (mc.cores = 1 and max.fork <= 1 both just loop in-process), so it is a weak check on bare R-loop overhead,
# not a fork comparison -- kept for completeness, not as the main signal.
run_phase(bs_num_max = 100, fork_values = c(2, 5, 10, 20), reps = 3, label = "B_scaling")
run_phase(bs_num_max = 20,  fork_values = 1,               reps = 3, label = "B_scaling_N1")

write.csv(results, "benchmarks/mclapply_vs_fork_bs_batch_results.csv", row.names = FALSE)

cat("\n=== Summary (mean over reps) ===\n")
agg <- aggregate(cbind(real_fork_sec, mclapply_sec) ~ label + bs_num_max + N, data = results, FUN = mean)
agg$speedup_mclapply_over_fork <- agg$real_fork_sec / agg$mclapply_sec
agg <- agg[order(agg$label, agg$N), ]
print(agg, row.names = FALSE)
