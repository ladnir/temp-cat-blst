# Middle Product Multipairing Status Report

## Executive Summary
- Date: April 9, 2026
- Branch: `codex/multipairing-plan`
- Base commit: `7a31e1203dd5c46bf7073a44cc417f9c9db5be3d`
- Working tree state: dirty

This workstream explored a direct multipairing-based middle-product path as a competitor to the existing FFT-based middle product.

The main conclusions so far are:
- The original naive direct path was the wrong performance baseline and has been replaced by real multipair-style accumulation.
- The scalar and AVX `Direct` paths now batch Miller work per output coefficient and do one final exponentiation per output.
- The `PrecomputedLines` path was initially performance-wrong because it finalized each pair term separately. That has now been fixed in the direct/reference path by delaying final exponentiation until the full output diagonal has been accumulated.
- On native AVX hardware, the direct multipair path is only competitive at small `n`. The online crossover against the FFT path is still around `n = 12..16`.
- On scalar, after fixing the delayed-final-exp lines path, the crossover is also around `n = 12..16`.
- The current code is good enough to benchmark and reason about, but it is not the final optimized implementation.

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

This is the core primitive used by the direct scalar multipair path.

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

## Benchmark Environment
### Remote host
- Host: `13.58.237.236`
- CPU: Intel(R) Xeon(R) Platinum 8488C
- Native AVX flags present, including `avx2`, `avx512f`, `avx512ifma`, `avx512bw`, `avx512vl`

### Build configuration
- Source dir used: `/tmp/blst_scalar_multipair_src`
- Build dir used: `/tmp/blst_scalar_multipair_build`
- Compiler: `clang 21.1.7`
- Build type: `Release`
- Key flags:
  - `BLST_AVX512=ON`
  - `BLST_DEFAULT_AVX=ON`
  - `BLST_AVX512_PORTABLE=OFF`
  - `BLST_TEST=ON`
  - `BLST_BENCHMARK=ON`

## Benchmark Results
### 1. Native AVX direct multipair vs FFT
These numbers are from the improved AVX direct path after adding true batched `miller_loop_n` and caller-side scratch reuse.

Compared benchmarks:
- FFT: `bench_middle_product_one_ifft_detailed`
- Direct AVX multipair: `bench_middle_product_direct_avx`

Online totals:

| n | FFT online ms | Direct AVX online ms | Direct / FFT |
|---|---:|---:|---:|
| 2 | 1.872 | 1.140 | 0.61x |
| 4 | 4.233 | 3.328 | 0.79x |
| 8 | 9.658 | 9.531 | 0.99x |
| 12 | 21.336 | 21.323 | 1.00x |
| 16 | 22.095 | 36.964 | 1.67x |
| 24 | 48.942 | 79.064 | 1.61x |
| 32 | 49.627 | 140.978 | 2.84x |
| 48 | 109.186 | 310.936 | 2.85x |
| 64 | 111.325 | 549.557 | 4.94x |

Interpretation:
- Direct AVX is competitive up to roughly `n = 12`.
- FFT clearly wins from `n >= 16`.

Important caveat:
- These AVX results were gathered before rerunning the newly fixed delayed-final-exp `PrecomputedLines` path on AVX. They describe the improved `Direct` AVX path, not the latest lines-mode state.

### 2. Native scalar FFT vs scalar direct vs scalar delayed-final-exp lines
These numbers were gathered after fixing `PrecomputedLines` to delay final exponentiation until the full output diagonal has been accumulated.

Compared benchmarks:
- Scalar FFT: `bench_middle_product_scalar_detailed`
- Scalar direct multipair: `bench_middle_product_direct_scalar`
- Scalar precomputed-lines multipair: `bench_middle_product_direct_scalar_lines`

Online totals:

| n | Scalar FFT ms | Scalar Direct ms | Scalar Lines ms |
|---|---:|---:|---:|
| 2 | 2.637 | 1.528 | 1.408 |
| 4 | 5.987 | 4.167 | 4.062 |
| 8 | 13.670 | 12.024 | 11.200 |
| 12 | 29.972 | 25.461 | 26.339 |
| 16 | 30.623 | 44.070 | 43.991 |
| 24 | 68.416 | 94.535 | 93.077 |
| 32 | 70.370 | 164.673 | 165.747 |
| 48 | 152.058 | 364.232 | 367.632 |
| 64 | 155.031 | 636.362 | 642.119 |

Interpretation:
- Scalar direct and scalar lines both beat scalar FFT for very small `n`.
- The scalar crossover is still around `n = 12..16`.
- The new delayed-final-exp lines path is modestly better than direct at `n = 2, 4, 8`, but it does not move the crossover much.

### 3. Offline precompute costs on scalar
Representative precompute costs:

| n | Scalar FFT precompute ms | Scalar Direct precompute ms | Scalar Lines precompute ms |
|---|---:|---:|---:|
| 8 | 8.922 | 0.093 | 1.324 |
| 16 | 22.822 | 0.215 | 2.582 |
| 32 | 55.281 | 0.400 | 5.149 |
| 64 | 132.424 | 0.595 | 10.318 |

Interpretation:
- Direct and lines have far cheaper offline precompute than FFT.
- The lines precompute cost is still much smaller than FFT precompute.
- This matters most for one-shot calls with a fresh `B`.

## What We Now Know
### Confirmed
- The direct multipair idea is real and measurable once implemented with batched Miller loops.
- The small-`n` crossover is real on both AVX and scalar.
- Fixing `PrecomputedLines` to delay final exponentiation was necessary and did help, but it did not create a dramatic new win.

### Not yet achieved
- A production-quality fixed-`B` multipair kernel.
- Reuse of overlap between adjacent output diagonals.
- A dedicated workspace API for the direct path.
- A fresh native AVX benchmark for the new delayed-final-exp lines path.

## Current Weaknesses / Open Performance Questions
### 1. No cross-output reuse
The current direct path computes each output diagonal independently. Adjacent outputs share almost all of their contributing pairs, but the implementation does not currently exploit that.

### 2. Reference-path placement
The code is still living in [test/MiddleProduct_reference.h](C:/Users/peter/.codex/worktrees/3234/blst/test/MiddleProduct_reference.h), which is fine for experimentation but not the final place for a tuned implementation.

### 3. AVX lines mode needs rerun
The new delayed-final-exp lines structure applies in the direct/reference path for all contexts, but AVX lines-mode benchmarks have not yet been rerun after that fix.

### 4. `update_direct` is still simple
`update_direct` still uses a straightforward per-output path. That is fine for now because the main performance interest has been the full middle product, not the constant-term update.

## Recommended Next Steps
### If continuing the multipair path
1. Rerun native AVX benchmarks for `bench_middle_product_direct_avx_lines` with the new delayed-final-exp lines path.
2. Decide whether the fixed-`B` lines path is worth pushing further.
3. If yes, the next big idea is not more local micro-optimization. It is structural reuse across adjacent output diagonals.
4. Move the direct path out of `test/` into a proper implementation/workspace layer only after the algorithm shape is judged worth keeping.

### If choosing a practical dispatch strategy now
- Small `n`:
  direct multipairing is plausible
- Medium and large `n`:
  FFT wins clearly

The current data suggests a threshold around:
- scalar: `n = 12..16`
- native AVX direct: `n = 12..16`

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
- the native benchmark data is consistent and useful

But the headline has not changed:
- direct multipairing is a small-`n` contender
- FFT remains the clear winner once `n` gets moderately large

If this path is pursued further, the next gains will have to come from exploiting shared structure across outputs, not from more local cleanup alone.
