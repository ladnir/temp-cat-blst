# AVX Batch Optimization (Contextualized Plan)

## Scope and current state

This plan is now aligned to the current repository state.

There are effectively two SIMD implementation styles in-tree:

1. **Shuffle-oriented packed schedules** (`4x2`, `2x2x2`, etc.).
2. **`8x1` lane-wise style** (the target for this plan).

For this workstream we only target the second style (`FpPack8x1` and tower code built on it).

Already implemented and validated:

1. `GTPack8x1` SoA API (`add`, `doubling`, `mult`, span schedulers).
2. GT packed scalar multiplication with GLV4 dispatch and masked basic fallback.
3. `Fp12Pack8x1::cyclotomic_sqr` specialized implementation.
4. `Fp12Pack8x1::frobenius_map` used by GT packed GLV basis generation.

What is **not** implemented yet:

1. `K`-wise interleaving/dispatch for `FpPack8x1` multiply+reduce.
2. Tower-level plumbing that explicitly consumes new `K`-wise `Fp8x1` kernels.

## Goal

Increase throughput of the `8x1` pipeline by running `K` independent `FpPack8x1` operations in flight, then measure scaling on real AVX-512 hardware.

Initial `K` set:

1. `K=1` baseline
2. `K=2`
3. `K=4`

(`K=8` may be evaluated later.)

## Phase 1 (next execution phase): `FpPack8x1` k-wise kernel

### APIs and file layout

Use template declarations in headers and explicit instantiations in source:

1. Add template declarations for `K`-wise kernel entry points in `blst/fields/Fp_avx.h`.
2. Implement templated bodies in `blst/fields/Fp_avx.cpp` (or a private included impl file if needed for readability).
3. Add explicit instantiations only for `K=1,2,4` in `Fp_avx.cpp`.
4. Provide thin non-template wrappers if needed for benchmark ergonomics.

Planned kernel set:

1. `template<int K> void fp8x1_mul_v2_k(UFpPack8x1* out, const FpPack8x1* a, const FpPack8x1* b);`
2. `template<int K> void fp8x1_redc_k(FpPack8x1* out, const UFpPack8x1* in);`
3. `template<int K> void fp8x1_mul_redc_k(FpPack8x1* out, const FpPack8x1* a, const FpPack8x1* b);`
4. `mul_add_redc` is intentionally deferred to a later phase.

Composition rule:

1. `mul_redc_k` is defined as `mul_v2_k` followed by `redc_k` using temporary widened values.
2. Keep `mul_v2_k` and `redc_k` independently callable for benchmarking and diagnosis.

Semantics:

1. `fp8x1_mul_v2_k`: `out[i] = mul_v2(a[i], b[i])` for `i in [0, K)`.
2. `fp8x1_redc_k`: `out[i] = redc(in[i])` for `i in [0, K)`.
3. `fp8x1_mul_redc_k`: `out[i] = redc(mul_v2(a[i], b[i]))` for `i in [0, K)`.
4. Alias policy should match existing style (safe aliasing where practical, otherwise assert/document no-overlap).

### Data layout and execution model

Keep current `8x1` representation:

1. One SIMD register per limb, 8 independent lanes.
2. No AoS/SoA conversion inside hot kernels.
3. Boundary conversion stays at API edges only.

Interleaving model:

1. Process K independent operations in a single invocation.
2. Start with k-wise wrappers around existing scalar-pack primitives (`mul_v2` and `redc`) to establish a low-risk baseline.
3. Use explicit template instantiation (`K=1,2,4`) and specialize bodies so the compiler can schedule independent chains aggressively.
4. Add runtime dispatch over larger arrays using greedy block order `4 -> 2 -> 1`.
5. If needed after measurement, move to deeper round-level interleaving inside `mul_v2`/`redc`.

Compiler optimization control policy (Phase 1):

1. Apply `noinline` on the `*_k` entry points only (not on core `mul_v2`/`redc` primitives).
2. This is an explicit control knob to limit inlining explosions, reduce accidental live-range growth, and make K-scaling behavior attributable.
3. Revisit this policy after native AVX-512 measurements.

### Testing plan

Unit tests:

1. Add canonical correctness tests for `K=1,2,4` for each primitive: `mul_v2_k`, `redc_k`, and `mul_redc_k`.
2. Compare against existing non-k-wise primitives over random and edge vectors.
3. Add aliasing tests where aliasing is supported.
4. Verify deterministic parity with current single-op path.

Failure policy:

1. If an aliasing mode is unsupported, enforce with `BLST_ASSERT` and test that contract indirectly.

### Benchmark plan

Add a focused benchmark for concept validation:

1. Measure `fp8x1_mul_v2_k` throughput for `K=1,2,4`.
2. Measure `fp8x1_redc_k` throughput for `K=1,2,4`.
3. Measure composed `fp8x1_mul_redc_k` throughput for `K=1,2,4`.
4. Report ops/s and relative speedup vs `K=1` for each kernel class.
5. Keep benchmark implementation deterministic and separate from GT benchmarks.

Important environment note:

1. Local dev machines without native AVX-512 are for correctness only.
2. Performance interpretation is only valid on native AVX-512 hardware.

### Success criteria

Phase 1 is successful if:

1. All new tests pass.
2. Existing AVX test suites remain green.
3. Benchmark builds and runs.
4. Native AVX-512 measurements show non-trivial speedup for at least one of `K=2` or `K=4` over `K=1` in at least one kernel class.
5. Results let us attribute gains/losses separately to multiply vs reduction.

## Follow-on phases (after Phase 1)

1. Thread `K`-wise `Fp8x1` kernels into `UFp2Pack8x1x1` and selected `Fp6/Fp12` operations.
2. Benchmark tower-level wins (`Fp2`, `Fp12 mul`, `Fp12 cyclotomic_sqr`).
3. Re-evaluate register-pressure cliff and whether `K=8` is viable.

## What to avoid

1. Do not add new shuffle-heavy logic to the `8x1` path.
2. Do not mix implementation styles in one hot path unless benchmark data justifies it.
3. Do not infer performance from emulation results.
4. Do not broaden scope (pairing/G1 pack integration) before `Fp8x1` `K`-wise results are in.
