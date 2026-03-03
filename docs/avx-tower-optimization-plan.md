# AVX Tower Throughput Optimization Plan (Fp -> Fp2 -> Fp6 -> Fp12 -> GT)

## Scope

This document defines a concrete optimization plan to recover SIMD scaling as we move up the field tower.

Target workload assumptions:

1. Native AVX-512 execution (`BLST_AVX512_PORTABLE=OFF`).
2. Independent tuples (not shared-scalar MSM style).
3. Native SoA pack path (`8x1`) as the mainline path.

GT notation reminder:

1. GT group operation uses additive notation at API level.
2. Implementation maps GT addition to multiplicative operations in cyclotomic `Fp12`.

## Problem Statement

Observed behavior:

1. `FpPack8x1` (`k=1`) is about 3x faster than 8 scalar `Fp` operations.
2. Tower throughput drops (for example near 2x in `Fp12`/GT paths) versus scalar x8.

Likely root cause:

1. Scalar `Fp12` path uses widened fused arithmetic (`__FP2x2__` path).
2. Pack canonical path reduces too frequently as it climbs the tower.
3. Extra reduction boundaries and scheduling overhead erode the base `Fp` gain.

Primary objective:

1. Preserve more widened work across `Fp2/Fp6/Fp12` before reduction.
2. Keep canonical `8x1` SoA layout end-to-end in hot code.
3. Improve GT throughput without changing existing scalar public APIs.

## Design Principles

1. Optimize based on measured attribution, not assumptions.
2. Keep canonical pack interfaces stable when possible.
3. Prefer fused wide compute over introducing shuffle-heavy detours.
4. Avoid AoS <-> SoA conversion in hot loops.
5. Defer `k`-wise interleaving changes until tower fusion gaps are closed.

## Phase 0: Attribution Baseline

### APIs

1. Add internal benchmark coverage for:
   1. `Fp2Pack8x1x1::mul`, `Fp2Pack8x1x1::sqr`
   2. `Fp6Pack8x1::mul`, `Fp6Pack8x1::sqr`
   3. `Fp12Pack8x1::mul`, `Fp12Pack8x1::sqr`, `Fp12Pack8x1::cyclotomic_sqr`
2. No public API changes.

### Data Layout

1. Compare SoA-to-SoA and scalar-to-scalar only.
2. No conversion overhead inside timed regions.

### Testing

1. Differential tests versus scalar reference for each operation.
2. Random vectors plus edge vectors.
3. Alias tests (`r == a`, `r == b`) where supported.

### Success Conditions

1. Per-layer speed table that shows where scaling drops.
2. Stable variance on native AVX-512 runs.

### What to Avoid

1. Mixing this phase with algorithmic refactors.
2. Drawing conclusions from portable-emulation timings.

## Phase 1: Add Wide Fp6 Primitives (Enabler)

### APIs

Add internal wide operations in `Fp6_avx.h/.cpp`:

1. `Fp6Pack8x1::mul_wide(UFp6Pack8x1&, const Fp6Pack8x1&, const Fp6Pack8x1&)`
2. `Fp6Pack8x1::sqr_wide(UFp6Pack8x1&, const Fp6Pack8x1&)`
3. `UFp6Pack8x1::add/sub/mul_by_nonresidue`

Keep existing reduced APIs unchanged.

### Data Layout

1. Reuse canonical `Fp6Pack8x1` and `UFp6Pack8x1`.
2. No new schedule/layout type introduced in this phase.

### Testing

1. Verify `redc(mul_wide(a,b)) == mul(a,b)`.
2. Verify `redc(sqr_wide(a)) == sqr(a)`.
3. Run randomized and adversarial vectors.

### Success Conditions

1. Bit-exact parity with current reduced outputs.
2. No regression in existing `Fp6`/`Fp12` tests.

### What to Avoid

1. Behavior changes in existing reduced routines.
2. Over-expanding temporary live ranges that force spills.

## Phase 2: Rebuild Fp12 mul/sqr with Fused Wide Flow

### APIs

Internal only, in `Fp12_avx.h/.cpp`:

1. `Fp12Pack8x1::mul_wide(UFp12Pack8x1&, const Fp12Pack8x1&, const Fp12Pack8x1&)`
2. Optional: `Fp12Pack8x1::sqr_wide(UFp12Pack8x1&, const Fp12Pack8x1&)`

`Fp12Pack8x1::mul` and `Fp12Pack8x1::sqr` remain public entry points and call fused internals.

### Data Layout

1. Stay in canonical `Fp12Pack8x1`.
2. Use `UFp12Pack8x1` intermediates to defer reductions.
3. Reduce at final output boundaries only where mathematically safe.

### Testing

1. Differential tests versus scalar `Fp12`.
2. Parity tests versus previous pack implementation.
3. Alias and tail-safety tests.

### Success Conditions

1. Throughput uplift in `Fp12` multiply/square versus current pack baseline.
2. No correctness regressions in GT pack operations that depend on `Fp12`.

### What to Avoid

1. Introducing AoS wrappers in hot pack path.
2. Prematurely adding `k`-wise complexity during this phase.

## Phase 3: Cyclotomic Square Specialization for GT Hot Path

### APIs

1. Keep `Fp12Pack8x1::cyclotomic_sqr` signature.
2. Rework internals to maximize wide/fused arithmetic reuse.

### Data Layout

1. Canonical SoA path only.
2. Any schedule kernel use must stay pack-native and avoid scalar-layout conversions.

### Testing

1. Differential tests against scalar cyclotomic square.
2. GT subgroup-consistency tests across randomized vectors.

### Success Conditions

1. Clear uplift in GT doubling throughput.
2. Correctness maintained under aliasing scenarios.

### What to Avoid

1. Route-through scalar `*_fp12_vector` conversion wrappers in pack hot loops.

## Phase 4: GT GLV4 Throughput Cleanup

### APIs

1. Keep current `GTPack8x1::mult` APIs.
2. Optionally add internal cached structures for repeated-base execution.

### Data Layout

1. Keep `GTPack8x1` as canonical SoA GT.
2. Keep basis/table values in `Fp12Pack8x1` form.

### Testing

1. One-shot and precompiled GLV4 parity tests.
2. Mixed random scalar widths and edge scalars.

### Success Conditions

1. Lower non-field overhead in GT scalar multiplication.
2. Stable gains in repeated-call scenarios.

### What to Avoid

1. Any assumption that scalars are shared across tuples.

## Phase 5: FFT-in-Exponent Path Integration

### APIs

1. Add pack-oriented helpers that accept spans of `GTPack8x1`.
2. Keep scalar API untouched.

### Data Layout

1. Prepack twiddles and operands into SoA once per stage batch.
2. No per-butterfly repacking.

### Testing

1. End-to-end FFT correctness against scalar reference:
   1. Forward FFT in exponent for G1/G2 pipeline dependencies.
   2. Inverse FFT in exponent for GT pipeline.
2. Round-trip properties and randomized vector checks.

### Success Conditions

1. End-to-end throughput tracks field-level improvements.
2. Conversion overhead is no longer dominant.

### What to Avoid

1. Hidden AoS <-> SoA transitions in butterfly loops.

## Phase 6: Aggressive Follow-ups (Only if Needed)

### APIs

1. Keep external APIs stable.
2. Add specialized internal kernels only where profiling shows bottlenecks.

### Data Layout

1. Reuse existing schedule types when justified.
2. Introduce broader accumulator-wide types only if range analysis requires it.

### Testing

1. Full differential suite and stress tests.
2. Performance counter review on native hardware.

### Success Conditions

1. Additional measurable gains after Phase 1-5 optimizations plateau.
2. No architecture regressions or correctness risk increase.

### What to Avoid

1. Complexity growth without measured wins.
2. Speculative refactors that reduce readability and maintainability.

## Implementation Order

1. Phase 0 (measurement baseline).
2. Phase 1 (wide `Fp6` primitives).
3. Phase 2 (`Fp12` fused wide flow).
4. Phase 3 (cyclotomic square optimization).
5. Phase 4 (GT GLV4 throughput cleanup).
6. Phase 5 (FFT-in-exponent integration).
7. Phase 6 only if profile-guided data justifies it.

## Exit Criteria

The program is complete when:

1. Tower-level speed loss source is fully attributed.
2. `Fp12`/GT pack throughput materially improves versus current baseline.
3. Existing correctness and API compatibility constraints are preserved.
4. FFT-in-exponent workload can run on native SoA pack kernels without conversion-heavy bottlenecks.
