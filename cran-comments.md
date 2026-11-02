# CRAN comments

## Test environment

* local macOS (aarch64-apple-darwin20), R 4.5.2, `R CMD check --as-cran`

## R CMD check results

0 errors | 2 warnings | 3 notes

This is a new release (first submission to CRAN).

### Warnings

1. **Installation warning: unknown warning group `-Wfixed-enum-extension`.**
   This comes from R's own header `R_ext/Boolean.h` being compiled with a newer Apple clang
   that doesn't recognize a warning flag R passes it. It is a local toolchain artifact, not
   something in this package's code, and is not expected to reproduce on CRAN's build machines.

2. **`checking compiled code` flags `_exit`/`_Exit` in `wspc.o`.**
   This is intentional, not a defect. `src/fork_utils.h` forks child processes to parallelize
   bootstrap resampling (`wisp()`'s `bootstraps.num` / `max.fork` options); each forked child
   calls `_exit()`, never `exit()`, to terminate without re-running atexit/cleanup handlers
   against the child's copy of memory -- running those handlers in the child is what's actually
   dangerous (double-flushed buffers, double-closed file descriptors, destructors firing against
   inconsistent post-fork state). `parallel::mclapply`'s own C implementation uses the identical
   pattern for the same reason. This path is `#ifdef`'d out entirely on Windows (no `fork()`
   there), so it only applies on Unix-alikes where it's correct.

   This also anticipates the natural follow-up question -- why hand-roll forking at all instead
   of calling `parallel::mclapply()`, an already-reviewed base-R primitive that does the same
   fork/collect pattern? This was benchmarked directly, not assumed. At small model sizes (276
   parameters, the scale of this package's own test fixtures), `mclapply()` is in fact 5-20%
   faster. But at the model sizes this package is built for in practice (1,000+ parameters,
   tests run up to 4,788), the result reverses: `mclapply()`'s long-lived worker processes
   accumulate heap-allocator-retained memory across repeated sequential Stan/NLopt fits (not a
   Stan autodiff leak -- `clear_stan_mem()` is called after every fit in both arms) and become
   6-16% slower while using up to 1.55x the peak combined memory.

   These are not cosmetic differences at the scale this package runs at in practice. The 16%
   figure is from a benchmark deliberately scaled down to 20 bootstraps (104.5s vs. 124.7s) to
   keep it tractable, but the package's own published validation (Barkasi et al., *Nucleic
   Acids Research*, 2026) ran 10,000 bootstrap resamples, taking roughly 205 minutes -- the same
   percentage gap at that bootstrap count is tens of extra minutes, not tens of seconds.
   Likewise the 1.55x memory figure is 57.7GB vs. 89.4GB in that same benchmark configuration;
   the package's own author already runs production models that combine 50+GB across concurrent
   bootstrap fits even with the current, more memory-efficient implementation, so a 1.55x
   multiplier at that scale is the difference between fitting on one large workstation and
   needing to shard across machines.

   `fork_utils.h` avoids this reversal by forking a fresh, short-lived batch of children per
   batch of resamples, so each child's entire address space -- including retained allocator
   pages -- is reclaimed by the OS on `_exit()`. Full reproducible methodology, benchmark
   scripts, and raw results are in the package's
   GitHub repository (see `URL` in DESCRIPTION) at `vignettes/articles/forking.Rmd`,
   `benchmarks/mclapply_vs_fork_general.R`, `benchmarks/RESULTS.md`, and
   `benchmarks/RESULTS_large_model.md` -- these are excluded from the CRAN tarball itself
   (via `.Rbuildignore`) to keep the submission small, so they will not be found inside it.

### Notes

1. `checking CRAN incoming feasibility ... NOTE` / "New submission" -- expected for a first
   submission.
2. `checking top-level files ... NOTE` -- "Files 'README.md' or 'NEWS.md' cannot be checked
   without 'pandoc' being installed." Local-machine-only; this machine lacks `pandoc`.
3. `checking HTML version of manual ... NOTE` -- HTML tidy on this machine is not recent enough
   to run validation, and the `V8` package isn't installed, so math-rendering checks were
   skipped locally. Also a local-tooling gap, not a known package defect.

## Windows NLopt linkage

An earlier version of this package linked against system NLopt via `pkg-config` on Unix
(`src/Makevars`) and a hardcoded `C:/local/nlopt` path on Windows (`src/Makevars.win`), with no
`configure.win` to fetch/stage that library -- a real risk on CRAN's Windows check machine. This
has since been replaced: NLopt is now obtained through the CRAN package `nloptr` (`LinkingTo`),
which bundles/builds its own copy of NLopt and exposes its C API via registered
`R_GetCCallable` routines (`src/wisp_nlopt.h`), so neither `Makevars` file links against
`-lnlopt` or any hardcoded path on either platform anymore -- see the comments in both files.
Checked on win-builder (R-release and R-devel) on 2026-11-02: both platforms came back with only
the standard "New submission" NOTE -- no errors, warnings, or Windows-specific notes -- confirming
this links and builds cleanly on Windows.
