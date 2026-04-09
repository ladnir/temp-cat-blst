# Middle Product Multipairing Status Report

## Executive Summary
- Date: April 9, 2026 (updated later the same day)
- Branch: `codex/multipairing-plan`
- Base commit: `7a31e1203dd5c46bf7073a44cc417f9c9db5be3d`
- Working tree state at update: clean

This workstream explored a direct multipairing-based middle-product path as a competitor to the existing FFT-based middle product.

The main conclusions so far are:
- The original naive direct path was the wrong performance baseline and has been replaced by real multipair-style accumulation.
- The scalar and AVX `Direct` paths now batch Miller work per output coefficient and do one final exponentiation per output.
- The `PrecomputedLines` path was initially performance-wrong because it finalized each pair term separately. That has now been fixed in the direct/reference path by delaying final exponentiation until the full output diagonal has been accumulated.
- A fresh native AVX rerun now includes the delayed-final-exp `PrecomputedLines` path. Online, `PrecomputedLines` is consistently about `22..24%` faster than AVX `Direct`, but the current FFT baseline still wins for every tested size except `n = 2`.
- For one-shot fixed-`B` calls that include offline precompute, native AVX `PrecomputedLines` is still plausible through roughly `n = 12`, but not beyond.
- A later scalar-only optimization pass added a batched scalar `miller_loop_lines_n` kernel and removed some wrapper overhead.
- That pass left scalar `Direct` essentially unchanged, but sped scalar `PrecomputedLines` up by about `1.17x..1.78x`.
- After that pass, scalar `PrecomputedLines` beats scalar FFT online through `n = 24` and wins one-shot fixed-`B` totals through `n = 32`.
- The current code is good enough to benchmark and reason about. The strongest remaining opportunities are structural reuse across adjacent output diagonals rather than more tiny wrapper cleanups.

## Original Goal
The goal was to evaluate whether the middle product should have a non-FFT path based on multipairing for small sizes.

The intended shape was:
- Keep the existing optimized middle-product formulation.
- Explore whether a direct path based on multipairing can beat FFT for small `n`.
- Pay special attention to the fixed-`B` case, where precomputing `G2` Miller-line data may help.

## Important Pivots
### 1. Naive direct path was a false start
An early direct implementation was just a nested loop over pair terms. That was not a real multipairing path and predictably lost badly.

That version should only be thought of as a correctness oracle, not as a candidate performance path.

### 2. The real direct path is “one multipair per output”
The direct path now computes each output coefficient by:
- gathering the diagonal pairs `(a_i, b_{n+k-i})`
- running one batched Miller loop over those pairs
- applying one final exponentiation at the end

### 3. `PrecomputedLines` had a major performance bug
Originally, the direct `PrecomputedLines` path evaluated precomputed lines one pair at a time and finalized each pair separately. In the dense case, that meant `n^2` final exponentiations.

That was performance-wrong.

The direct/reference `PrecomputedLines` path now does the intended thing:
- evaluate precomputed line tables in the Miller domain
- multiply Miller contributions across the output diagonal
- final exponentiate once per output

## Current Code Shape
### Files touched
- [test/MiddleProduct_reference.h](C:/Users/peter/.codex/worktrees/3234/blst/test/MiddleProduct_reference.h)
- [test/MiddleProduct_direct_tests.cpp](C:/Users/peter/.codex/worktrees/3234/blst/test/MiddleProduct_direct_tests.cpp)
- [benchmark/MiddleProduct_bench.cpp](C:/Users/peter/.codex/worktrees/3234/blst/benchmark/MiddleProduct_bench.cpp)
- [benchmark/FFT_bench.cpp](C:/Users/peter/.codex/worktrees/3234/blst/benchmark/FFT_bench.cpp)
- [blst/avx/curve/GT.h](C:/Users/peter/.codex/worktrees/3234/blst/blst/avx/curve/GT.h)

### Scalar direct path
In [test/MiddleProduct_reference.h](C:/Users/peter/.codex/worktrees/3234/blst/test/MiddleProduct_reference.h):
- `precompute_middle_product_b_direct` precomputes backend `B`, affine `B`, zero mask, and optional prepared lines.
- `middle_product_direct` converts `A` once, makes affine `A` once, compacts nonzero `A`, and then computes each output coefficient.
- `direct_multipair_window` implements the `Direct` mode.
- `direct_multipair_lines_window` implements the fixed-`B` `PrecomputedLines` mode with delayed final exponentiation.

### Scalar backend primitive
The scalar backend batched Miller engine is in [GT.cpp](C:/Users/peter/.codex/worktrees/3234/blst/blst/scalar/curve/GT.cpp):
- `MILLER_LOOP_N_MAX = 16`
- `start_dbl_n`
- `add_n_dbl_n`
- `GT::miller_loop_n`
- `GT::miller_loop_lines_n`

These are the core primitives used by the scalar direct and scalar fixed-`B` lines paths.

### AVX direct path
In [GT.h](C:/Users/peter/.codex/worktrees/3234/blst/blst/avx/curve/GT.h), `miller_loop_n` was upgraded from:
- looping over single-pair `miller_loop` calls

to:
- a real chunked batched Miller schedule

This made the AVX direct numbers much more representative.

## Correctness Validation Performed
### Local validation
Focused tests passed locally after the current changes:
- `MiddleProductDirect.ScalarMatchesScalarFft`
- `MiddleProductDirect.ScalarPrecomputedLinesMatchesScalarFft`
- `MiddleProductDirect.ScalarUpdateAndCombineReconstructsFullResult`
- `MiddleProductDirect.AvxPortableMatchesDefaultFft`
- `MiddleProductDirect.AvxPortablePrecomputedLinesMatchesDefaultFft`
- `MiddleProductCApi.*`

### Remote Linux validation
On April 9, 2026, the scalar tests were rerun on the EC2 AVX machine and passed:
- `MiddleProductDirect.ScalarMatchesScalarFft`
- `MiddleProductDirect.ScalarPrecomputedLinesMatchesScalarFft`
- `MiddleProductDirect.ScalarUpdateAndCombineReconstructsFullResult`

Later the same day, the AVX-oriented checks were rerun on the same machine and also passed:
- `MiddleProductDirect.AvxPortableMatchesDefaultFft`
- `MiddleProductDirect.AvxPortablePrecomputedLinesMatchesDefaultFft`
- `MiddleProductCApi.*`

## Benchmark Environment
### Remote host
- Host: `13.58.237.236`
- CPU: Intel(R) Xeon(R) Platinum 8488C
- Native AVX flags present, including `avx2`, `avx512f`, `avx512ifma`, `avx512bw`, `avx512vl`

### Build configuration
- Source dir used: `/tmp/blst_multipair_48cd776_src`
- Build dir used: `/tmp/blst_multipair_48cd776_build`
- Compiler: `clang 21.1.7`
- Build type: `Release`
- Key flags:
  - `BLST_AVX512=ON`
  - `BLST_DEFAULT_AVX=ON`
  - `BLST_AVX512_PORTABLE=OFF`
  - `BLST_TEST=ON`
  - `BLST_BENCHMARK=ON`

## Benchmark Results
### 1. Fresh native AVX online comparison
These numbers are from the later April 9 rerun on the current branch state. This rerun includes the delayed-final-exp `PrecomputedLines` AVX path.

Compared benchmarks:
- FFT: `bench_middle_product_one_ifft_detailed`
- Direct AVX multipair: `bench_middle_product_direct_avx`
- Direct AVX delayed-final-exp lines: `bench_middle_product_direct_avx_lines`

Online totals:

| n | FFT online ms | Direct AVX ms | AVX Lines ms | Direct / FFT | Lines / FFT |
|---|---:|---:|---:|---:|---:|
| 2 | 1.025 | 1.140 | 0.909 | 1.11x | 0.89x |
| 4 | 2.329 | 3.261 | 2.568 | 1.40x | 1.10x |
| 8 | 5.357 | 9.386 | 7.322 | 1.75x | 1.37x |
| 12 | 11.942 | 21.060 | 16.433 | 1.76x | 1.38x |
| 16 | 12.263 | 36.663 | 28.628 | 2.99x | 2.33x |
| 24 | 27.124 | 78.373 | 59.974 | 2.89x | 2.21x |
| 32 | 27.771 | 139.089 | 107.207 | 5.01x | 3.86x |
| 48 | 60.865 | 307.121 | 235.975 | 5.05x | 3.88x |
| 64 | 62.141 | 540.278 | 413.425 | 8.69x | 6.65x |

Interpretation:
- `PrecomputedLines` is the better AVX direct kernel. Across most of the tested range it is about `0.76..0.80x` the cost of AVX `Direct`.
- Online, the current FFT baseline wins for every tested size except `n = 2`, where AVX `PrecomputedLines` is modestly faster.
- The direct AVX numbers are close to the earlier snapshot, but the current `bench_middle_product_one_ifft_detailed` rerun is materially faster. Practical dispatch decisions should therefore use this newer table rather than the older AVX-only comparison.

### 2. Native AVX one-shot totals including offline precompute
If `B` is fresh and the offline precompute matters, the fixed-`B` lines path still has a more interesting niche.

Representative totals (`offline precompute + online compute`):

| n | FFT total ms | Direct AVX total ms | AVX Lines total ms | Lines / FFT |
|---|---:|---:|---:|---:|
| 2 | 1.450 | 1.173 | 1.323 | 0.91x |
| 4 | 3.612 | 3.320 | 3.280 | 0.91x |
| 8 | 8.823 | 9.496 | 8.810 | 1.00x |
| 12 | 21.034 | 21.231 | 18.726 | 0.89x |
| 16 | 21.317 | 36.894 | 31.678 | 1.49x |
| 24 | 49.796 | 78.756 | 65.121 | 1.31x |
| 32 | 50.507 | 139.581 | 114.003 | 2.26x |
| 48 | 114.975 | 307.823 | 248.135 | 2.16x |
| 64 | 116.327 | 541.202 | 425.836 | 3.66x |

Interpretation:
- For one-shot fixed-`B` calls, AVX `PrecomputedLines` is still attractive at `n = 2, 4, 12` and effectively tied at `n = 8`.
- The one-shot fixed-`B` crossover for AVX lines is still only around `n = 12..16`.
- This is the strongest remaining practical case for keeping the lines path under consideration.

### 3. Later scalar optimization pass: batched scalar lines kernel
Later on April 9, a scalar-focused optimization pass made two changes:
- Added `blst::scalar::GT::miller_loop_lines_n`, a batched scalar fixed-`B` kernel that shares the Miller schedule across many pairs instead of evaluating one full lines loop per pair term.
- Removed some wrapper overhead in the experimental direct path by avoiding per-output affine repacking and by reusing `a0`/`b_affine` in `update_direct`.

The benchmark comparison below uses the same EC2 AVX host, but still exercises the scalar backend.

Representative online totals after the scalar pass:

| n | Scalar FFT ms | Scalar Direct ms | Scalar Lines ms |
|---|---:|---:|---:|
| 2 | 2.625 | 1.547 | 1.204 |
| 4 | 5.823 | 4.291 | 3.061 |
| 8 | 13.330 | 11.624 | 7.880 |
| 12 | 29.973 | 25.884 | 16.353 |
| 16 | 30.583 | 46.064 | 27.433 |
| 24 | 67.496 | 93.162 | 57.150 |
| 32 | 69.638 | 165.061 | 99.320 |
| 48 | 153.750 | 362.404 | 213.166 |
| 64 | 156.305 | 640.909 | 371.553 |

Representative `Scalar Lines` before/after totals from the scalar pass:

| n | Lines before ms | Lines after ms | Speedup |
|---|---:|---:|---:|
| 2 | 1.410 | 1.204 | 1.17x |
| 4 | 4.086 | 3.061 | 1.34x |
| 8 | 11.721 | 7.880 | 1.49x |
| 12 | 25.896 | 16.353 | 1.58x |
| 16 | 44.520 | 27.433 | 1.62x |
| 24 | 93.506 | 57.150 | 1.64x |
| 32 | 166.875 | 99.320 | 1.68x |
| 64 | 644.322 | 371.553 | 1.73x |
| 256 | 10185.641 | 5733.120 | 1.78x |

Interpretation:
- The wrapper cleanups did not materially move scalar `Direct`; the benchmark noise is around break-even.
- The batched scalar lines kernel is the real win. It turns `PrecomputedLines` from "slightly better than direct at tiny `n`" into the clear best scalar non-FFT path.
- Online, scalar `Direct` still only beats FFT through roughly `n = 12`, but scalar `PrecomputedLines` now beats FFT through `n = 24` and loses at `n = 32`.
- A later scalar-only cleanup removed the remaining preprocessed-lines chunking in the hot path and precomputed the `A`-side line-eval form once. That pass was smaller, but still worth about `3..4%` around `n = 24..128`, with the scalar online crossover remaining between `n = 24` and `n = 32`.

### 4. Scalar one-shot totals including offline precompute
Because scalar FFT has much heavier `B` precompute than the direct kernels, the one-shot fixed-`B` story is more favorable than the online-only story.

Representative totals (`offline precompute + online compute`):

| n | Scalar FFT total ms | Scalar Direct total ms | Scalar Lines total ms |
|---|---:|---:|---:|
| 2 | 3.751 | 1.571 | 1.549 |
| 4 | 9.128 | 4.334 | 3.622 |
| 8 | 22.360 | 11.712 | 9.081 |
| 12 | 52.762 | 26.040 | 18.272 |
| 16 | 53.073 | 46.264 | 29.869 |
| 24 | 123.656 | 93.471 | 61.118 |
| 32 | 124.973 | 165.483 | 104.619 |

Representative precompute costs from the same run:

| n | Scalar FFT precompute ms | Scalar Direct precompute ms | Scalar Lines precompute ms |
|---|---:|---:|---:|
| 8 | 9.029 | 0.088 | 1.201 |
| 16 | 22.490 | 0.200 | 2.436 |
| 24 | 56.160 | 0.309 | 3.968 |
| 32 | 55.335 | 0.422 | 5.299 |
| 64 | 132.525 | 0.790 | 10.386 |

Interpretation:
- The scalar fixed-`B` lines path now wins one-shot totals through `n = 32`.
- Scalar `Direct` also remains one-shot-competitive through about `n = 24`, but `PrecomputedLines` is the stronger scalar option.
- For scalar, the offline FFT cost is now a major part of why fixed-`B` direct kernels remain interesting.

## What We Now Know
### Confirmed
- The direct multipair idea is real and measurable once implemented with batched Miller loops.
- On scalar, the small-`n` crossover is still real.
- The later scalar optimization pass showed that the real scalar headroom was in batching the lines path, not in squeezing the already-batched direct path.
- The fresh native AVX rerun shows that `PrecomputedLines` is consistently better than AVX `Direct`, but the current FFT baseline still wins online for `n >= 4`.
- For one-shot fixed-`B` calls, AVX `PrecomputedLines` remains plausible only through roughly `n = 12`.
- Fixing `PrecomputedLines` to delay final exponentiation was necessary and did help, but it did not by itself create a broad AVX online win.
- On scalar, the later batched lines kernel does materially change the practical picture: online `PrecomputedLines` is now viable through `n = 24`, and one-shot totals through `n = 32`.

### Not yet achieved
- A production-quality fixed-`B` multipair kernel.
- Reuse of overlap between adjacent output diagonals.
- A dedicated workspace API for the direct path.

## Current Weaknesses / Open Performance Questions
### 1. No cross-output reuse
The current direct path computes each output diagonal independently. Adjacent outputs share almost all of their contributing pairs, but the implementation does not currently exploit that.

### 2. Reference-path placement
The code is still living in [test/MiddleProduct_reference.h](C:/Users/peter/.codex/worktrees/3234/blst/test/MiddleProduct_reference.h), which is fine for experimentation but not the final place for a tuned implementation.

### 3. AVX benchmark story changed after the rerun
The fresh AVX rerun completed the missing lines-mode measurement, but it also showed the current `bench_middle_product_one_ifft_detailed` baseline materially faster than the earlier note. The direct AVX numbers stayed in roughly the same range, so the practical online AVX dispatch threshold moved down sharply.

### 4. `update_direct` is still simple
`update_direct` still uses a straightforward per-output path. The scalar wrapper pass cleaned up some avoidable conversions here, but benchmark impact was minimal. The main performance interest is still the full middle product, not the constant-term update.

## Recommended Next Steps
### If continuing the multipair path
1. Reconcile or stabilize the AVX FFT baseline so dispatch thresholds are based on a reproducible comparison, not on one stale table and one fresh rerun.
2. Decide whether the fixed-`B` lines path is worth keeping specifically for one-shot or very small-`n` use. The online-only case no longer justifies it by itself.
3. For scalar, the fixed-`B` lines path is now strong enough that it should be treated as the leading non-FFT candidate rather than as a curiosity.
4. If the path is kept, the next big idea is not more local micro-optimization. It is structural reuse across adjacent output diagonals.
5. Move the direct path out of `test/` into a proper implementation/workspace layer only after the algorithm shape is judged worth keeping.

### If choosing a practical dispatch strategy now
- Scalar:
  `Direct` crosses over around `n = 12..16`, while `PrecomputedLines` now crosses over around `n = 24..32`
- Native AVX, online with precomputed `B` already available:
  FFT wins from `n >= 4`
- Native AVX, one-shot fixed-`B` call:
  `PrecomputedLines` is still plausible through roughly `n = 12`

## Fast Ramp-Up Checklist For Later
When returning to this work, read in this order:
1. This report
2. [middle_product_multi_pairing.tex](C:/Users/peter/.codex/worktrees/3234/blst/docs/middle_product_multi_pairing.tex)
3. [MiddleProduct_reference.h](C:/Users/peter/.codex/worktrees/3234/blst/test/MiddleProduct_reference.h)
4. [GT.cpp](C:/Users/peter/.codex/worktrees/3234/blst/blst/scalar/curve/GT.cpp)
5. [GT.h](C:/Users/peter/.codex/worktrees/3234/blst/blst/avx/curve/GT.h)
6. [MiddleProduct_bench.cpp](C:/Users/peter/.codex/worktrees/3234/blst/benchmark/MiddleProduct_bench.cpp)
7. [FFT_bench.cpp](C:/Users/peter/.codex/worktrees/3234/blst/benchmark/FFT_bench.cpp)

## Bottom Line
The direct multipairing experiment is now in a much better state than when it started:
- the direct path is a real multipair path
- the scalar `PrecomputedLines` path now has the right delayed-final-exp structure
- the missing native AVX lines rerun is now done

The updated headline is narrower than the earlier snapshot suggested:
- online on native AVX, `PrecomputedLines` is better than `Direct`, but FFT still wins for every tested size except `n = 2`
- for one-shot fixed-`B` AVX calls, `PrecomputedLines` still has a plausible niche through about `n = 12`
- on scalar, `Direct` still only wins online through about `n = 12`, but the later batched `PrecomputedLines` kernel wins online through `n = 24` and wins one-shot totals through `n = 32`

If this path is pursued further, the next gains will have to come from exploiting shared structure across outputs, not from more local cleanup alone.
