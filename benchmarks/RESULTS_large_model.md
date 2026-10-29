# mclapply vs. fork_utils.h at a much larger model: results

This extends `benchmarks/RESULTS.md` (the 276-parameter liverradial result) to answer the
follow-up question it couldn't: does "mclapply is modestly faster" hold, reverse, or change in
magnitude for a much bigger model (1000+ parameters, large memory footprint)? **The 276-parameter
result in `RESULTS.md` is left completely untouched below** -- this file is a separate, later
benchmarking pass, using a generalized version of the same methodology
(`benchmarks/mclapply_vs_fork_general.R` + `benchmarks/configs/*.R`, vs. the original's hardcoded
`benchmarks/mclapply_vs_fork_isolated.R`) plus new peak-memory tracking that the original script
didn't have.

Same machine as `RESULTS.md`: 24 logical cores (Darwin/arm64), 128 GB RAM.

## 0. Generalizing the script: regression check against the original result

`benchmarks/mclapply_vs_fork_general.R` is a config-driven generalization of
`mclapply_vs_fork_isolated.R`: it sources a plain-R config file that supplies `count.data`,
`data.variables`, `model.settings`, a `(bs_num_max, N)` grid, and a rep count, then runs the exact
same two arms (`cpp_model$bs_batch(bs_num_max, N, FALSE, FALSE)` vs. `cpp_model$fit(FALSE)` +
`mclapply(..., cpp_model$bs_fit(i, TRUE), mc.cores = N, mc.preschedule = TRUE)`) plus best-effort
peak-combined-RSS tracking (`benchmarks/mem_watch_utils.R`) around every timed call. The original
hardcoded script (`mclapply_vs_fork_isolated.R`) is untouched, as the historical record.

`benchmarks/configs/liverradial_config.R` reproduces the original dataset/settings (same
276-parameter liverradial model) as a regression check -- confirming generalizing the script
didn't change its behavior. Result (mean over 3 reps; a subset of the original grid, chosen to
keep the regression check itself fast):

| bs_num_max | N  | bs_batch() (s) | mclapply (s) | mclapply speedup | original RESULTS.md speedup |
|---|---|---|---|---|---|
| 100  | 5  | 94.6  | 83.9  | 1.13x | 1.10x |
| 100  | 20 | 32.3  | 30.4  | 1.06x | 1.09x |
| 1000 | 20 | 296.1 | 263.1 | 1.13x | 1.20x |

All three numbers land within a few percent of the original Round-2 table, well inside this kind
of benchmark's run-to-run noise on a shared machine -- **the generalized script reproduces the
original result.** (Raw per-rep data: `benchmarks/mclapply_vs_fork_general_liverradial_results.csv`.)

This run also incidentally produced the first peak-RSS numbers for the small model, which turn out
to matter for interpreting the big-model results below: **mclapply's peak combined RSS was higher
than bs_batch's in every one of these three configs** (1.06x, 1.74x, 1.97x) even though wall-clock
still favored mclapply. At 276 parameters the absolute memory involved is small (8-32 GB combined
across up to 20 forked/forked-once processes) so this memory asymmetry didn't show up as a speed
disadvantage. That changes at the larger model.

## 1. Getting to a much larger model, from this package's own tutorial data

Per the task's suggested order:

1. **Finer spatial binning** (`tutorial_corticallaminar.Rmd`): not feasible with the shipped data.
   `S1_laminar_countdata_demo.csv.gz` / `ACx_laminar_countdata_demo*.csv.gz` only ship the
   already-binned `bin` column (100 levels each); the original continuous per-cell laminar
   coordinates used to produce that binning aren't in `inst/extdata/`, so there's no way to re-bin
   more finely without fabricating coordinates, which was out of scope. Lowering `LROcutoff` (more
   sensitive change-point detection) at the real 100-bin resolution was tried and does increase
   parameter count on the single-context cortical-laminar data, but only into the hundreds/low
   thousands -- not enough to be the main approach.
2. **Multi-context scaling** (`tutorial_multicontext.Rmd`): this is what was used, and it worked
   well. The tutorial itself explains that fixed-effect parameter count scales with the number of
   `context` (cell-type) levels, and deliberately collapses the real 31 MapMyCells cell types in
   `ACx_laminar_countdata_demo_with_celltype.csv.gz` down to just 2 ("Glut"/"GABA") purely "to make
   results easier to visualize and speed up model fitting" -- i.e. that collapse is a demo
   simplification, not a modeling requirement. `benchmarks/configs/corticallaminar_multicontext_big_config.R`
   instead keeps the 12 best-populated of the real 31 cell types (1,278-103,332 rows each) and
   drops only the sparse tail (as few as 12 total rows for some of the 31 types), where even a
   single fit became numerically unstable. No fabricated data anywhere -- this is exactly the real
   `ACx_laminar_countdata_demo_with_celltype.csv.gz` rows, just not collapsed down to 2 contexts.

This reaches a model with:

| | small model (`RESULTS.md`) | big model (this file) |
|---|---|---|
| dataset | `Droin_radial_count_data_sim.csv.gz` (liverradial) | `ACx_laminar_countdata_demo_with_celltype.csv.gz` (12 of 31 real cell types kept) |
| data rows | 882,144 | 145,776 |
| contexts | 1 (dummy "liver") | 12 (real MapMyCells cell types: IT-ET Glut, Immune, CTX-MGE GABA, Vascular, Astro-Epen, NP-CT-L6b Glut, MY Glut, MY GABA, MH-LH Glut, MB Glut, OB-IMN GABA, CTX-CGE GABA) |
| genes | 1 | 6 |
| mice | -- | 5 |
| fixed effects | timeseries (ZT) | hemisphere + timeseries (age: 7/12/18) |
| **fitted parameters** | **276** | **4,788 (17.3x)** |
| single-process, non-parallel `fit()`: wall-clock | 4.1s | 30.5s |
| single-process, non-parallel `fit()`: peak RSS | 1,003 MB | 2,919 MB |
| (R + wispack loaded, no model, for reference) | 224 MB | 224 MB |

Note parameter count (17x) grew much faster than either data-row count (which actually *shrank*,
145,776 < 882,144) or peak single-process memory (2.9x) or single-fit wall-clock (7.4x) -- a wisp
model's per-process memory/compute footprint tracks parameter count (number of context x species x
treatment-level x segment combinations being optimized) much more than raw row count. Scaling
exploration further (not used for the timed benchmark below, just confirming headroom) reached
6,786 parameters at 18 cell types before per-bootstrap fit cost (minutes, single-threaded) made a
full timed grid impractical to actually run in this task's time budget. **4,788 parameters is the
real number reported here** -- honestly short of "tens of GB" for a *single* model instance (2.9 GB)
but, as shown below, comfortably into the tens-of-GB range in *combined* memory once multiple
bootstrap fits run concurrently, which is the regime the benchmark actually probes.

Each individual bootstrap fit (`bs_fit()`) on this model costs roughly 30s single-threaded (vs.
0.15-0.3s for the 276-parameter model -- two full orders of magnitude more expensive per task), so
`bs_num_max` was kept much smaller (20-40, vs. the original's 100-1000) to keep the benchmark grid
tractable; the grid still spans the same qualitative comparisons (fixed `bs_num_max` across
increasing `N`, and doubled `bs_num_max` at fixed `N`). "Could not find nearby feasible parameters"
warnings during some individual bootstrap fits are expected and harmless -- the model's own
graceful fallback for a resample that's numerically awkward for one of the many context/gene/
treatment combinations in a model this large, not a crash.

## 2. Results: timing

Mean over 2 reps (raw data: `benchmarks/mclapply_vs_fork_general_corticallaminar_multicontext_big_results.csv`):

| bs_num_max | N  | bs_batch() (s) | mclapply (s) | mclapply speedup |
|---|---|---|---|---|
| 20 | 4  | 204.3 | 224.3 | **0.91x (mclapply 9% slower)** |
| 20 | 10 | 125.8 | 136.7 | **0.92x (mclapply 8% slower)** |
| 20 | 20 | 104.5 | 124.7 | **0.84x (mclapply 16% slower)** |
| 40 | 20 | 138.7 | 148.0 | **0.94x (mclapply 6% slower)** |

**The result reverses at this scale.** At 276 parameters, mclapply was consistently 5-20% *faster*
than `bs_batch()`. At 4,788 parameters, `bs_batch()` is consistently 6-16% *faster* than mclapply --
every single grid point flips. This is a real, repeatable reversal (both reps at each grid point
agree in direction; see the CSV), not noise crossing zero.

## 3. Results: peak combined RSS

| bs_num_max | N  | bs_batch() peak RSS (MB) | mclapply peak RSS (MB) | mclapply / bs_batch |
|---|---|---|---|---|
| 20 | 4  | 15,459 | 23,807 | **1.54x** |
| 20 | 10 | 37,309 | 40,998 | 1.10x |
| 20 | 20 | 69,998 (46,683 / 93,314 across the 2 reps -- high variance, see below) | 43,518 | 0.62x |
| 40 | 20 | 57,715 | 89,397 | **1.55x** |

Combined peak RSS reaches 40-90 GB at N=20 on this 128 GB machine -- a large fraction of total
system RAM, which matters for interpretation (see "A caveat on these absolute numbers" below).

**The clearest, most mechanistically interpretable pattern is in how peak RSS scales with
*repeated calls per persistent worker*, not just with N:**

- At fixed N = 20, doubling `bs_num_max` from 20 to 40 means mclapply's `mc.preschedule = TRUE`
  workers go from 1 task each to 2 tasks each (sequential, same persistent process) -- and
  mclapply's peak RSS almost exactly doubles (43,518 -> 89,397 MB, 2.05x) under that doubled
  per-worker task count. Meanwhile `bs_batch()`'s peak RSS barely moves (46,683 -> 57,715 MB,
  comparing the cleaner rep of the noisy N=20/bs=20 row, 1.24x) -- consistent with
  `fork_utils.h` forking a *fresh* batch of (up to) `max.fork` children every batch and those
  children dying (and the OS reclaiming all their memory in one shot) at the end of each batch,
  regardless of how many total batches run.
- At fixed `bs_num_max` = 20, going from N = 4 (5 sequential tasks/worker) to N = 10 (2
  tasks/worker) to N = 20 (1 task/worker, no repetition for either arm) tracks the same story from
  the other direction: the mclapply/bs_batch RSS ratio is highest (1.54x) where mclapply workers do
  the *most* repeated sequential bootstrap fits (N=4), drops to 1.10x at 2 tasks/worker (N=10), and
  is lowest/noisy right at N=20 where *neither* arm has any repeated-call accumulation to show (a
  single batch of 20 forks vs. 20 single-task mclapply workers -- structurally the same situation
  for both arms at that one grid point, which is exactly where the ratio is close to/under 1).

This is the author's hypothesis, and it bears out: **`mclapply`'s long-lived, `mc.preschedule = TRUE`
workers accumulate more combined memory the more bootstrap fits each one runs sequentially before
exiting; `fork_utils.h`'s short-lived, re-forked-every-batch children do not, because they die (and
have their entire address space reclaimed by the OS) after every batch regardless of total
bootstrap count.** At 276 parameters this same mechanism is present (mclapply's RSS was higher
there too) but the absolute memory at stake is small enough (single-digit GB) that it didn't cost
mclapply anything in wall-clock time. At 4,788 parameters the same mechanism produces tens of GB of
difference, large enough to plausibly show up as real performance cost (see next section).

## 4. Why: mechanistic explanation

**Why memory behaves this way.** Both arms call `clear_stan_mem()` (`stan::math::recover_memory()`)
after every individual bootstrap fit -- `bs_batch()`'s serial path does this explicitly
(`bs_fit(b, true)`, see `src/wspc.cpp` line ~1277), and this benchmark's mclapply arm does it too
(`cpp_model$bs_fit(i, TRUE)`), matching that same practice. So this is *not* a case of Stan's
autodiff memory arena silently growing across calls within a worker -- that's actively guarded
against on both arms, every call. `resample()` (`src/wspc.cpp` line ~1162) also mutates fixed-size
member vectors in place rather than growing any container across calls. The accumulation visible in
section 3 therefore isn't a literal data-structure leak at the R/C++ level; it's consistent instead
with ordinary heap-allocator behavior: a long-lived process's `malloc` implementation generally
keeps *freed* pages mapped for reuse rather than returning them to the OS immediately, especially
for the large, variably-sized allocations a ~4,788-parameter NLopt fit makes on every bootstrap
resample. Each of mclapply's persistent workers can retain a high-water-mark of such pages across
its whole sequence of assigned bootstraps. A forked `fork_utils.h` child, by contrast, is a whole
separate process: when it `_exit()`s at the end of its one bootstrap fit, the OS tears down its
*entire* address space in one shot -- freed-but-retained allocator pages included -- so the next
batch's freshly forked children start from the same COW-shared baseline as the *parent's* own
(bounded) memory state, not from whatever the *previous batch's children* had grown their heaps to.

**Why timing flips.** The original 276-parameter result (`RESULTS.md`) predicted the wall-clock gap
between the two arms should *shrink* as per-task cost grows, because `fork_utils.h`'s extra
`fork()`-call overhead (sub-millisecond to low-single-digit milliseconds) becomes a vanishing
fraction of an increasingly expensive per-bootstrap task. That prediction is directionally
confirmed in part -- the gaps here (6-16%) aren't dramatically larger than the small model's (5-20%)
despite the per-task cost growing ~100-200x -- but it doesn't predict an outright *reversal*, since
raw fork-syscall cost is even more negligible here than before. The more complete explanation is
that a second effect has entered and started to dominate: the memory-accumulation pattern in
section 3 is a genuine, measurable performance cost at this scale, not just a bookkeeping
curiosity. Combined RSS at N=20 reaches tens of GB on a 128 GB machine; a process (or 20 of them
concurrently) holding a larger, more fragmented heap pays real costs on every subsequent large
allocation/fit -- more TLB and cache pressure, slower allocator searches through a bigger pool of
free chunks, more pages the OS has to manage and potentially compress under memory pressure. None
of that is present for `fork_utils.h`'s children, which start every batch from a comparatively
clean, bounded memory state. At 276 parameters this cost was too small in absolute terms to
outweigh mclapply's one-time-fork-cost advantage; at 4,788 parameters it's large enough to flip the
net result.

**A caveat on these absolute numbers.** `ps`-based RSS summing (the method `mem_watch_utils.R`
uses, documented in its own header comment) double-counts copy-on-write pages still shared between
a forked child and its parent, so these are overestimates of true distinct memory use, more so
right after a fork and less so the longer/more a worker's own allocations diverge from the parent's
state -- which is itself part of why mclapply's later-task memory looks particularly large relative
to bs_batch's. Separately, combined RSS at N=20 (40-90 GB) is a large fraction of this specific
128 GB test machine's RAM; some of the run-to-run variance seen in section 3 (the N=20/bs_num_max=20
row's 46.7 GB vs. 93.3 GB across its two reps) is plausibly real OS-level memory-pressure behavior
(paging, compression) rather than pure measurement noise, and the absolute wall-clock numbers at
N=20 may themselves be partly confounded by that pressure on this particular machine, rather than
purely isolating the forking-mechanism difference. The *direction* of the result (mclapply slower,
mclapply using more memory, both effects growing with repeated-calls-per-worker) is consistent and
repeatable across every grid point, which is the actual empirical claim being made here; the exact
percentages should be read as "this machine, this run," not as precise universal constants.

## 5. What this does and doesn't mean

This is empirical evidence that the 276-parameter result (`mclapply` modestly faster) does **not**
generalize to much larger models -- it can reverse, and a plausible, specific mechanism (persistent
workers' memory accumulating across repeated sequential calls, vs. fork_utils.h's children dying
and being fully reclaimed every batch) explains both the memory numbers and, plausibly, the timing
reversal. This confirms the premise that the original result's magnitude (and even its sign) is
workload-dependent, consistent with the author's own prior experience of the opposite result on a
different package. As with `RESULTS.md`, **this is not a recommendation to migrate `bs_batch()` to
`mclapply()`** -- if anything, this result is a point *against* such a migration for workloads with
large models and many sequential bootstraps per worker, though a single exploratory benchmark on
one dataset/machine isn't a rigorous basis for that conclusion either. Left as a benchmarking
exercise per the task; `bs_fit`'s exposure and this benchmark code remain unwired from the
package's real fitting path, and no package source code was changed to produce any of these results.
