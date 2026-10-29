# Best-effort peak-RSS tracking for the generalized mclapply-vs-fork_utils.h benchmark.
#
# Idea: while a timed call is running, a lightweight background shell process polls `ps` every
# `interval` seconds, sums resident-set-size (RSS, in KB) across the main R process and any
# process whose *parent* pid is the main R process (i.e. its direct children -- which covers both
# fork_utils.h's forked batch workers and mclapply's forked mc.cores workers, since both fork
# directly off the calling R process), and keeps a running max. We read that running max after the
# timed call ends, then tear the watcher down.
#
# This is explicitly best-effort: if `ps`/background-process launching doesn't work on this
# platform for some reason, every function here degrades to returning NA rather than erroring, so
# a memory-tracking failure never aborts the actual timing benchmark.
#
# Known limitation (stated here once, referenced from RESULTS write-ups rather than repeated):
# `ps` RSS double-counts copy-on-write pages shared between a forked child and its parent -- a
# freshly forked child's RSS includes pages it shares with (and hasn't diverged from) the parent,
# which is not really "extra" memory. Summing RSS across parent + children therefore overstates
# true combined memory use, more so right after a fork and less so the more a child's working set
# diverges from the parent's. There's no cheap, portable alternative to RSS from a plain `ps` call
# (proper shared-memory accounting needs PSS, which `ps` does not report on macOS), so this
# benchmark uses summed RSS as a best-effort proxy and flags the caveat rather than pretending it
# isn't there.

.mem_watch_script <- "
#!/bin/bash
main_pid=\"$1\"
outfile=\"$2\"
stopfile=\"$3\"
interval=\"$4\"
max_kb=0
while [ ! -f \"$stopfile\" ]; do
  sum=$(ps -eo pid,ppid,rss= 2>/dev/null | awk -v m=\"$main_pid\" '($1==m || $2==m){s+=$3} END{print s+0}')
  if [ -n \"$sum\" ]; then
    if [ \"$sum\" -gt \"$max_kb\" ] 2>/dev/null; then
      max_kb=\"$sum\"
    fi
  fi
  echo \"$max_kb\" > \"$outfile\"
  sleep \"$interval\"
done
"

#' Start a background peak-RSS watcher for the current R process + its direct children.
#' @param interval poll interval in seconds.
#' @return a list "handle" to pass to stop_mem_watcher(); handle$ok is FALSE if the watcher could
#'   not be started (e.g. `ps`/bash unavailable), in which case stop_mem_watcher() will return NA.
start_mem_watcher <- function(interval = 0.3) {
  handle <- list(ok = FALSE, outfile = NULL, stopfile = NULL, script = NULL)
  tryCatch({
    script_path <- tempfile("memwatch_", fileext = ".sh")
    writeLines(.mem_watch_script, script_path)
    Sys.chmod(script_path, "0755")
    outfile  <- tempfile("memwatch_out_")
    stopfile <- tempfile("memwatch_stop_")
    main_pid <- Sys.getpid()
    system2(
      "bash", args = c(script_path, main_pid, outfile, stopfile, interval),
      wait = FALSE, stdout = FALSE, stderr = FALSE
    )
    handle <- list(ok = TRUE, outfile = outfile, stopfile = stopfile, script = script_path)
  }, error = function(e) {
    message("mem watcher: failed to start (", conditionMessage(e), "); memory tracking disabled for this call.")
  })
  handle
}

#' Stop a watcher started with start_mem_watcher() and return the peak combined RSS observed, in MB.
#' Returns NA_real_ if the watcher never started or no reading could be recovered.
stop_mem_watcher <- function(handle) {
  peak_mb <- NA_real_
  if (isTRUE(handle$ok)) {
    tryCatch({
      Sys.sleep(0.5) # let the watcher take one more sample before we tell it to stop
      if (file.exists(handle$outfile)) {
        val_kb <- suppressWarnings(as.numeric(readLines(handle$outfile, warn = FALSE)[1]))
        if (length(val_kb) == 1 && !is.na(val_kb)) peak_mb <- val_kb / 1024
      }
    }, error = function(e) {
      message("mem watcher: failed to read peak RSS (", conditionMessage(e), ").")
    })
    tryCatch(file.create(handle$stopfile), error = function(e) NULL)
    tryCatch(unlink(c(handle$outfile, handle$stopfile, handle$script)), error = function(e) NULL)
  }
  peak_mb
}

#' Run `expr_fun()`, timing it and (best-effort) tracking peak combined RSS of the main process and
#' its direct children throughout the call.
#' @return list(sec = wall-clock seconds, peak_rss_mb = peak combined RSS in MB, or NA)
time_and_track_mem <- function(expr_fun, poll_interval = 0.3) {
  handle <- start_mem_watcher(poll_interval)
  t0 <- Sys.time()
  expr_fun()
  sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  peak_rss_mb <- stop_mem_watcher(handle)
  list(sec = sec, peak_rss_mb = peak_rss_mb)
}
