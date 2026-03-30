# AVX Tower Micro-Benchmark Plan (Code-Only)

## Scope

This plan is intentionally limited to code and micro-benchmarks for tower operations.

Target operations:

1. `Fp12` multiplication.
2. `Fp12` squaring.
3. `Fp12` cyclotomic squaring.

Supporting attribution operations:

1. `Fp2Pack8x1x1` multiply/square.
2. `Fp6Pack8x1` multiply/square.

Out of scope for this document:

1. GT GLV4 scheduling work.
2. FFT integration work.
3. API integration beyond local overloads used by tower kernels.

## Optimization Target

Primary target is relative parity with base-field SIMD gains:

1. Achieve `>= 3.0x` throughput for `Fp12Pack8x1` ops versus scalar `Fp12` reference on native AVX-512.

Secondary targets:

1. Keep bit-exact correctness versus scalar reference.
2. Keep canonical `8x1` SoA hot path.
3. Avoid extra AoS <-> SoA conversions inside timed regions.

## Benchmark Rules

Environment assumptions:

1. Native AVX-512 machine (`BLST_AVX512_PORTABLE=OFF`).
2. Fixed CPU frequency/governor where possible.
3. Release build only.

Measurement rules:

1. Benchmark scalar and pack paths on equivalent batches of 8 independent inputs.
2. Report `items_per_second` and explicit speedup (`pack_ops_per_sec / scalar_ops_per_sec`).
3. Keep conversion, random generation, and setup out of timed loops.
4. Use repetitions and aggregate reporting for stability.

Suggested benchmark families:

1. `BM_Fp12Mul_Scalar8` vs `BM_Fp12Mul_Pack8x1`.
2. `BM_Fp12Sqr_Scalar8` vs `BM_Fp12Sqr_Pack8x1`.
3. `BM_Fp12CyclotomicSqr_Scalar8` vs `BM_Fp12CyclotomicSqr_Pack8x1`.
4. Attribution: `Fp2` and `Fp6` scalar8 vs pack8x1 pairs.

## Phase 0: Baseline Micro-Bench Harness

### Code

1. Add dedicated micro-bench cases for the operation families above.
2. Ensure scalar baselines call scalar implementations directly (no runtime dispatch ambiguity).
3. Ensure pack baselines are canonical `8x1` kernels.

### Tests

1. Differential tests for each benchmarked operation against scalar reference.
2. Random vectors and deterministic edge vectors.

### Success

1. Stable benchmark outputs with low variance.
2. Clear baseline table for `Fp2/Fp6/Fp12`.

## Phase 1: Fp6 Overload Refactor for Wide Fusion

### Code

Introduce overloads in `Fp6_avx.h/.cpp`:

1. `Fp6Pack8x1::mul(Fp6Pack8x1&, const Fp6Pack8x1&, const Fp6Pack8x1&)` (existing reduced).
2. `Fp6Pack8x1::mul(UFp6Pack8x1&, const Fp6Pack8x1&, const Fp6Pack8x1&)` (new widened output).
3. `Fp6Pack8x1::sqr(Fp6Pack8x1&, const Fp6Pack8x1&)` (existing reduced).
4. `Fp6Pack8x1::sqr(UFp6Pack8x1&, const Fp6Pack8x1&)` (new widened output).
5. Add `UFp6Pack8x1` helpers needed for fused arithmetic (`add/sub/mul_by_nonresidue`).

### Tests

1. Verify `redc(widened overload) == reduced overload` for mul/sqr.
2. Verify alias behavior remains valid for reduced overloads.

### Success

1. No correctness regressions.
2. Measurable reduction in `Fp6` overhead in attribution benchmarks.

## Phase 2: Fp12 Fused Overloads Using Wide Fp6

### Code

Introduce overloads in `Fp12_avx.h/.cpp`:

1. `Fp12Pack8x1::mul(Fp12Pack8x1&, const Fp12Pack8x1&, const Fp12Pack8x1&)` (existing reduced entry).
2. `Fp12Pack8x1::mul(UFp12Pack8x1&, const Fp12Pack8x1&, const Fp12Pack8x1&)` (new widened output).
3. Optional same pattern for `sqr`.

Implementation direction:

1. Build in widened domain using `UFp6Pack8x1` intermediates.
2. Defer reduction until output boundaries.
3. Keep canonical `8x1` representation end-to-end.

### Tests

1. Differential parity against scalar `Fp12`.
2. Parity against previous pack reduced path.
3. Alias and tail safety checks.

### Success

1. Significant uplift in `Fp12` mul/sqr benchmarks.
2. Move aggregate speedup closer to `3.0x` target.

## Phase 3: Cyclotomic Square Tightening

### Code

1. Keep `Fp12Pack8x1::cyclotomic_sqr` API unchanged.
2. Rework internals to maximize fused operations and minimize reduction boundaries.
3. Reuse new overloads where they lower total reduction count.

### Tests

1. Differential parity with scalar `Fp12::cyclotomic_sqr`.
2. Random and structured subgroup inputs.

### Success

1. Cyclotomic square speedup approaches the same relative target trend as mul/sqr.

## Phase 4: Benchmark-Driven Final Tuning

### Code

1. Tighten inlining policy and temporary lifetimes for hot kernels.
2. Remove obvious register-pressure pathologies.
3. Keep semantics and API signatures stable.

### Tests

1. Re-run full micro-benchmark matrix.
2. Re-run differential correctness suite.

### Success

1. `Fp12` operations achieve or exceed `3.0x` speedup versus scalar reference.
2. No correctness regressions.

## Exit Criteria

This plan is complete when:

1. All targeted `Fp12` micro-bench operations (`mul`, `sqr`, `cyclotomic_sqr`) meet `>= 3.0x` versus scalar `Fp12`.
2. Attribution benchmarks show where remaining losses are (if any).
3. Differential tests are green for all changed kernels.
