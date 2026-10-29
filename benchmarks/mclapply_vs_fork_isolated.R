# Follow-up to mclapply_vs_fork_bs_batch.R: that script timed the *whole* wisp() call against the mclapply
# prototype, which confounds the fork-mechanism comparison with wisp()'s own non-bootstrap overhead (LRO grid
# search, predict_rates_batch, sample.stats, residual analysis, plot generation -- none of which the mclapply
# prototype does). Confirmed this confound is real and non-trivial: even at N = 1 (neither method forks) the
# two arms differed by ~20%, which can only be that extra overhead, not a forking-mechanism effect.
#
# This script isolates just the bootstrap-fitting step itself: cpp_model$bs_batch() (the real fork_utils.h
# path) vs. mclapply() over cpp_model$bs_fit() calls, both starting from the same freshly LRO-initialized
# cpp_model and both including their own one-time full-data fit() call (bs_batch() does this internally;
# the mclapply arm does it explicitly), with nothing else around either call.

devtools::load_all(".", quiet = TRUE)
suppressMessages(library(parallel))

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

time_call <- function(expr_fun) {
  t0 <- Sys.time()
  expr_fun()
  as.numeric(difftime(Sys.time(), t0, units = "secs"))
}

count_data <- read.csv(system.file("extdata", "Droin_radial_count_data_sim.csv.gz", package = "wispack"))
data.variables <- list(count = "count", bin = "bin", context = "liver", species = "gene", ran = "mouse", timeseries = "ZT")
model.settings <- list(LRO_cost = "none", LROcutoff = 1.5)

results <- data.frame()
reps <- 3

for (bs_num_max in c(1000, 100)) {
  fork_values <- if (bs_num_max == 1000) 20 else c(2, 5, 10, 20)
  for (N in fork_values) {
    for (rep in seq_len(reps)) {
      cat(sprintf("bs_num_max=%d N=%d rep=%d/%d: bs_batch()... ", bs_num_max, N, rep, reps))
      cpp_model_a <- capture_cpp_model(
        count.data = count_data, variables = data.variables, model.settings = model.settings,
        plot.settings = list(print.plots = FALSE), verbose = FALSE
      )
      t_fork <- time_call(function() cpp_model_a$bs_batch(bs_num_max, N, FALSE, FALSE))
      cat(sprintf("%.1fs. mclapply... ", t_fork))

      cpp_model_b <- capture_cpp_model(
        count.data = count_data, variables = data.variables, model.settings = model.settings,
        plot.settings = list(print.plots = FALSE), verbose = FALSE
      )
      t_mcl <- time_call(function() {
        cpp_model_b$fit(FALSE)
        rows <- parallel::mclapply(
          seq_len(bs_num_max) - 1L,
          function(i) cpp_model_b$bs_fit(i, TRUE),
          mc.cores = N, mc.preschedule = TRUE
        )
        do.call(rbind, rows)
      })
      cat(sprintf("%.1fs.\n", t_mcl))

      results <- rbind(results, data.frame(
        bs_num_max = bs_num_max, N = N, rep = rep, bs_batch_sec = t_fork, mclapply_sec = t_mcl
      ))
    }
  }
}

write.csv(results, "benchmarks/mclapply_vs_fork_isolated_results.csv", row.names = FALSE)

cat("\n=== Isolated bootstrap-step-only comparison (mean over reps) ===\n")
agg <- aggregate(cbind(bs_batch_sec, mclapply_sec) ~ bs_num_max + N, data = results, FUN = mean)
agg$speedup_mclapply_over_fork <- agg$bs_batch_sec / agg$mclapply_sec
agg <- agg[order(agg$bs_num_max, agg$N), ]
print(agg, row.names = FALSE)
