# mclapply vs. fork_utils.h: results

Machine: 24 logical cores (Darwin/arm64). Dataset: `Droin_radial_count_data_sim.csv.gz` (same data/settings
as `tests/testthat/test-slow-liverradial.R`, 276 fitted parameters, ~0.15-0.2s per bootstrap fit). 3 repetitions
per configuration throughout; variance across reps was small (see the raw CSVs).

## Round 1: whole `wisp()` call vs. mclapply prototype (`mclapply_vs_fork_bs_batch.R`)

| bs_num_max | N (max.fork / mc.cores) | real wisp() (s) | mclapply prototype (s) | mclapply speedup |
|---|---|---|---|---|
| 1000 | 20 | 201.2 | 157.5 | 1.28x |
| 100  | 2  | 144.4 | 130.0 | 1.11x |
| 100  | 5  | 68.7  | 55.5  | 1.24x |
| 100  | 10 | 40.1  | 31.4  | 1.28x |
| 100  | 20 | 27.3  | 18.9  | 1.45x |
| 20   | 1  | 60.3  | 50.1  | 1.20x |

**This round is confounded.** The mclapply prototype only does `fit()` + the bootstrap loop; the real `wisp()`
call also does `predict_rates_batch()`, `sample.stats()`, residual analysis, and plot generation -- none of
which involve forking at all. At N = 1 (bs_num_max = 20), *neither* arm forks, yet the "real" arm is still 20%
slower -- that 20% is pure non-bootstrap wisp() overhead, not a forking-mechanism effect. So round 1 overstates
mclapply's advantage by conflating "less work done" with "faster forking."

## Round 2: bootstrap step only, isolated (`mclapply_vs_fork_isolated.R`)

Same cpp_model (freshly LRO-initialized, not yet fit) fed to both arms. Fork arm: `cpp_model$bs_batch(bs_num_max,
N, FALSE, FALSE)` -- the real production code path, including its own internal full-data `fit()`. mclapply
arm: `cpp_model$fit(FALSE)` then `mclapply(0:(bs_num_max-1), function(i) cpp_model$bs_fit(i, TRUE), mc.cores = N,
mc.preschedule = TRUE)`. Nothing else runs in either arm.

| bs_num_max | N  | bs_batch() (s) | mclapply (s) | mclapply speedup |
|---|---|---|---|---|
| 100  | 2  | 137.1 | 130.2 | 1.05x |
| 100  | 5  | 61.3  | 55.8  | 1.10x |
| 100  | 10 | 33.7  | 32.1  | 1.05x |
| 100  | 20 | 21.7  | 19.9  | 1.09x |
| 1000 | 20 | 189.0 | 157.1 | **1.20x** |

**This is the real, apples-to-apples number: mclapply is consistently faster, by 5-20%, never slower, on this
codebase.** That's the opposite of what was seen on a different package, and a real, repeatable result on this
one -- not close to noise (reps agreed to within ~1-2s on multi-minute runs).

## Why

`fork_utils.h`'s `run_forked_batches()` forks a fresh batch of up to `max.fork` children *every batch*, so the
total number of `fork()` calls for a `bs_batch()` run is always `bs_num_max`, independent of `max.fork` -- 20
bootstraps per batch at `max.fork = 20` still means 1000 total forks for 1000 bootstraps (50 batches of 20),
and only 2 forks per batch at `max.fork = 2` still means 100 total forks for 100 bootstraps (50 batches of 2).
`mclapply(mc.preschedule = TRUE)` instead forks exactly `mc.cores` long-lived workers once, each of which
processes its share of the task list in a sequential in-process loop and reports back once.

This predicts the one place the two rounds' numbers clearly diverge: at the same `N = 20`, going from
`bs_num_max = 100` (5 batches, 100 total forks) to `bs_num_max = 1000` (50 batches, 1000 total forks) widens
the gap from 1.09x to 1.20x. `fork_utils.h`'s per-fork cost (the `fork()` syscall itself, pipe setup/teardown,
`waitpid()`) scales with total bootstrap count; `mclapply`'s fork cost is fixed at `mc.cores` and amortizes
away as the per-worker task list grows.

The gap is real but genuinely modest (5-20%, not multiples) -- and that modesty is itself explained by the
workload: each bootstrap fit here costs ~0.15-0.3s of real CPU work (NLopt + Stan-math autodiff), which is
large relative to a `fork()` call's actual cost (typically sub-millisecond to low-single-digit milliseconds
on macOS, even before counting pipe I/O). Forking overhead is a small fraction of each unit of work here, so
there isn't much headroom for *any* parallelism strategy to win big by forking less. On a different package
with much finer-grained per-task work (sub-millisecond to low-millisecond tasks), the same fork-count
difference would be a much larger fraction of total time -- and `mclapply`'s own fixed costs (task-list
serialization to send the work description to workers, `parallel:::selectChildren`-style result polling) could
plausibly dominate and flip the result, which is consistent with the opposite result seen there. In other
words: this isn't "mclapply is just better" -- it's "the fork-count difference matters more or less depending
on how expensive each unit of forked work is," and on this codebase's actual per-bootstrap cost, it modestly
favors mclapply.

## What this does and doesn't mean for `bs_batch()`

This is empirical evidence that a real `mclapply`-based rewrite of `bs_batch()` could plausibly be 5-20% faster
on workloads like this package's real tutorials/tests (hundreds to ~1000 bootstraps, hundreds of parameters).
It is NOT a recommendation to make that change: `bs_batch()`'s current fork_utils.h path is well-tested,
`mclapply` carries its own platform caveats (no forking on Windows either, same as now, but different failure
modes/warnings), and a real migration would need to actually produce the same result matrix shape bs_batch()
does (this benchmark's mclapply arm does not assemble the full-data-fit row bs_batch() includes) and get its
own correctness verification, not just a speed number. Left as a benchmarking exercise per the task; `bs_fit`'s
exposure and this benchmark code are not wired into the package's real fitting path.
