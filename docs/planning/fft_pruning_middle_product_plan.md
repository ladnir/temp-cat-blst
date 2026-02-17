# FFT Pruning Plan for Middle Product Pipeline

## 1) Goal and Scope

This plan defines how to add NTL-style truncated/pruned FFT ideas to BLST's middle-product application path where:

- `A` has length `n` in `G1`
- `B` has length `2n` in `G2`
- We only need the middle window of the convolution result

We will keep two callable implementations:

1. Full/basic FFT-based middle extraction (baseline, stable)
2. Specialized middle-product implementation (high-effort optimized path)

No auto-selection is included in this phase.


## 2) Current Baseline (as of this plan)

Current code already exposes:

- `MiddleProductAlgorithm::{BasicFFT, SpecializedFFT}`
- Optional precompute of `B`
- Optional timing capture via `TimingStats`

Current baseline uses power-of-two FFT sizing and includes full transforms for several steps. Performance can be jagged around size boundaries due to power-of-two rounding.


## 3) NTL Concepts to Reuse

Reference concept from NTL:

- Truncated FFT API with forward `(yn, xn)`:
  - `xn` = number of significant input coefficients
  - `yn` = number of output values to compute
- Recursive pruning skips branches that cannot affect requested outputs.
- "Admissible size" rounding (`FFTRoundUp`) to avoid tiny truncations that do not help and can regress.

Adapting this to BLST:

- Reuse truncated recursion logic and admissible-size policy.
- Keep group-specific operations and pairings in BLST representation.
- Maintain explicit CPU-only implementation now, with clear future CUDA extension points.


## 4) Problem-Specific Sizing Rules

For `|A| = n`, `|B| = 2n`:

- Full convolution support size: `3n - 1`
- Full/basic FFT path requires:
  - `N_full = pow2_ceil(3n - 1)`
- Specialized middle-product path may use a smaller working FFT size depending on algorithm internals (currently often near `2n`).

Planned sizing utilities:

- `fft_pow2_size_for_full_middle(n)` -> `N_full`
- `fft_admissible_size(x, k)` -> NTL-like rounded size for truncation (`xn`, `yn`)
- `fft_should_round_to_full(x, N)` heuristic for near-full truncation case


## 5) API and Integration Changes

### 5.1 New truncated FFT entry points (CPU)

Add explicit interfaces in `blst/fft/ctFft.h` (or adjacent helper header) for:

- `fft_truncated_forward(..., k, yn, xn, roots, ...)`
- `fft_truncated_inverse(..., k, yn, roots, inv_roots, ...)`

Constraints:

- `N = 2^k`
- `1 <= yn <= N`
- `1 <= xn <= N`
- `xn` and `yn` adjusted to admissible sizes before recursion

### 5.2 Middle product API behavior

Keep external APIs stable (`middle_product(...)`, precompute variants), but wire internals as:

- Basic FFT mode:
  - Use truncated forward on `A` (`xn = n`)
  - Use truncated forward/precompute on `B` (`xn = 2n`)
  - Keep pointwise multiply across required spectrum
  - Use inverse path and middle extraction unchanged semantically
- Specialized mode:
  - Replace current recursion blocks with explicit output-pruned strategy for middle window
  - Use admissible rounding to reduce jagged transitions

### 5.3 Timing hooks

Extend timer labels for pruning-specific visibility:

- `*.trunc_setup_ns`
- `*.trunc_forward_a_ns`
- `*.trunc_forward_b_ns`
- `*.trunc_inverse_ns`
- `*.pruned_butterfly_ns`
- `*.admissible_rounding_ns`

This enables side-by-side comparison with existing baseline timing keys.


## 6) Algorithm Work Plan (Phased)

## Phase A: Truncated Forward Infrastructure (Low Risk)

Deliverables:

- Implement generic truncated forward recursion in CPU path.
- Use it in Basic FFT path for `A` and `B`.
- Keep inverse path unchanged initially if needed for safety.

Expected result:

- Immediate reduction in work from known-zero tails (`A[n..]`, `B[2n..]`) without changing result semantics.

## Phase B: Admissible-Size Policy (Stability and Smoothness)

Deliverables:

- Add NTL-style admissible rounding helper.
- Round `xn`/`yn` to SIMD/cache-friendly multiples with near-full snap-to-`N`.

Expected result:

- Reduced benchmark jaggedness around size boundaries.
- Fewer regressions where "tiny truncation" is slower than full.

## Phase C: Specialized Middle-Window Pruning (High Value)

Deliverables:

- Rework specialized middle-product recursion to compute only influence cone for middle coefficients.
- Combine input truncation and output truncation logic in one recursion framework.
- Preserve current functional contract and output order.

Expected result:

- Largest speedup for workloads where only middle window is needed.

## Phase D: Cleanup and Hardening

Deliverables:

- Consolidate shared recursion helpers.
- Remove duplicate branches in specialized/basic internals.
- Add documentation blocks and invariants in code comments.


## 7) Correctness and Testing Plan

Primary oracle (already aligned with your requirement):

- Compute expected middle window with `O(n^2)` field arithmetic.
- Project to group and compare against FFT-based outputs.

Test suites to add/extend in `test/Fft_tests.cpp`:

1. `Fft.MiddleProductPrunedEqualsNaive_SmallDeterministic`
2. `Fft.MiddleProductPrunedEqualsNaive_Randomized`
3. `Fft.MiddleProductBasicTruncatedEqualsBasicFull`
4. `Fft.MiddleProductPrecomputeB_PruneAndNonPruneMatch`
5. `Fft.MiddleProductAdmissibleRoundingBoundaryCases`
6. `Fft.MiddleProductNonPowerOfTwoN`

Coverage dimensions:

- `n` across powers of two and non-powers
- odd/even `n`
- boundary sizes near `2^k`, `2^k +/- 1`, `3n - 1` transitions
- both algorithms and both `B` modes (precomputed vs non-precomputed)


## 8) Benchmark Plan

Use Google Benchmark filter to keep runtime short and focused:

- Basic full (current baseline)
- Basic with truncated forward
- Specialized (current)
- Specialized with pruning

Suggested benchmark dimensions:

- `n` in geometric progression with extra boundary points
- separate runs for:
  - online-only (precomputed `B`)
  - full pipeline (including precompute)

Outputs to track:

- total latency
- per-stage timer breakdown
- crossover points where specialized beats basic


## 9) CUDA Forward Plan (Not Implemented Now)

CPU remains source of truth. To keep CUDA path obvious:

- Introduce API-level placeholders mirroring truncated interfaces.
- Keep TODO markers near the new truncated APIs and middle-product entrypoints.
- Design recursion helpers so GPU kernels can replace inner butterfly blocks later without changing API contracts.

No CUDA optimization is done in this phase.


## 10) Risks and Mitigations

Risk: pruning logic bugs in edge indices.
Mitigation: exhaustive small-`n` tests against naive `O(n^2)` oracle.

Risk: truncation overhead exceeds gains for some sizes.
Mitigation: admissible rounding and benchmark-driven tuning.

Risk: specialized path diverges from basic semantics.
Mitigation: direct specialized-vs-basic equivalence tests over broad size sweep.


## 11) Execution Checklist

1. Implement truncated forward recursion and wire Basic FFT path.
2. Add admissible-size rounding utility and integrate everywhere truncation is used.
3. Rework specialized middle-product recursion for true output-pruning.
4. Add/expand correctness tests and boundary cases.
5. Add benchmark variants + timer keys and collect crossover data.
6. Document CPU complete / CUDA pending status in code comments and docs.


## 12) Definition of Done

Done when all are true:

- Both callable implementations remain available and correct.
- New pruning logic passes deterministic + randomized oracle tests.
- Boundary tests for admissible rounding pass.
- Benchmark report shows clear stage-level timing and crossover behavior.
- CPU path is production-ready; CUDA follow-up points are explicit and tracked.

