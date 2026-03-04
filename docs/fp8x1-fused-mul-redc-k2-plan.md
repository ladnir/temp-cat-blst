# Fp8x1 Fused `mul+redc` K=2 Plan

## Goal

Implement a true fused `FpPack8x1` `mul+redc` path for `k=2` streams in native AVX-512, with internal-only APIs and no public API changes.

Primary objective:
- increase `mul+redc` throughput by avoiding widened intermediate materialization and reducing dependency/scheduling stalls.

## Scope

In scope:
- `FpPack8x1` / `UFpPack8x1` only.
- native AVX-512 (`BLST_AVX512`, non-portable), Linux x86_64 asm path.
- `k=2` first, then reuse for `k=4` as two `k=2` calls.

Out of scope (for this doc):
- public API/ABI additions.
- AVX portable optimization work.
- direct GT/G1/G2 refactors.

## Current Baseline (for comparison)

Representative native AVX benchmark ranges observed:
- `BM_Fp8x1MulRedc_SeqChains<2>`: ~`11.3M..11.5M` ops/s.
- `BM_Fp8x1MulRedc_HandInterleavedChains<2>`: near-flat vs seq (small gain).
- `redc_core` improved substantially with asm `k=2`, but that gain does not transfer fully to full `mul+redc` today.

Interpretation:
- current bottleneck is at composed kernel level (`mul` + reduction boundary + scheduling), not just isolated `redc_core`.

## Internal API Plan (No Public API)

### 1) Core fused primitive (internal)

```cpp
BLST_AVX_NOINLINE void fp8x1_mul_redc_core_k2(
    LFpPack8x1& out0,
    LFpPack8x1& out1,
    const FpPack8x1& a0,
    const FpPack8x1& b0,
    const FpPack8x1& a1,
    const FpPack8x1& b1);
```

Behavior:
- computes `(a0*b0)*R^-1 mod p` and `(a1*b1)*R^-1 mod p` in lazy form (`LFpPack8x1`).
- no final canonical subtraction.

### 2) Canonical fused wrapper (internal)

```cpp
BLST_FORCEINLINE void fp8x1_mul_redc_pair(
    FpPack8x1& out0,
    FpPack8x1& out1,
    const FpPack8x1& a0,
    const FpPack8x1& b0,
    const FpPack8x1& a1,
    const FpPack8x1& b1);
```

Behavior:
- calls `fp8x1_mul_redc_core_k2`, then `normalize` for both outputs.

### 3) Dispatch integration

`fp8x1_mul_redc_block<K>`:
- `K==1`: unchanged baseline C++ path.
- `K==2`: use `fp8x1_mul_redc_pair`.
- `K==4`: two `fp8x1_mul_redc_pair` calls.

## Data Layout and Scheduling Targets

### Data layout
- keep existing SoA pack layout (`FpPack8x1`, 8 limbs of `__m512i`).
- two independent streams (`stream0`, `stream1`) processed together.

### Scheduling targets
- avoid storing full widened `UFp` intermediates to memory in the hot path.
- keep enough independent work across streams to hide IFMA latency.
- control live set to avoid spill/reload cycles.
- preserve radix-48 arithmetic model used by current `redc_core`.

### Register budget strategy
- start with `k=2` only.
- explicitly bound active temporaries per round.
- prefer deterministic temp ownership per stream/step (avoid accidental long chains).

## Phased Implementation

## Phase 1: Scaffolding and Baseline Hooks

APIs:
- add internal declarations for `fp8x1_mul_redc_core_k2` + `fp8x1_mul_redc_pair`.
- add compile guards mirroring `blst_fp8x1_redc_core_k2`.

Data/layout:
- no arithmetic changes yet.
- fallback implementation uses existing `mul_v2` + current reduction path.

Testing:
- existing `fp8x1_mul_redc_k` tests must pass unchanged.
- add fused-path differential tests (random vectors + edge vectors).

Success:
- zero functional regression.
- fused hook selectable in `K==2` path under guard.

Avoid:
- no public headers/API exposure.
- no benchmark-only logic in library code.

## Phase 2: C++ Fused Baseline (Control)

APIs:
- wire `K==2` through `fp8x1_mul_redc_pair` in C++ first.

Data/layout:
- use two local widened values (`wide0`, `wide1`) as control baseline.

Testing:
- compare with old path on randomized and adversarial inputs.
- run existing Fp/Fp2/Fp12 tests that transitively depend on Fp mul.

Success:
- correctness parity.
- benchmark baseline for fused-entry path established (even if perf-neutral).

Avoid:
- no premature asm schedule complexity before control baseline is locked.

## Phase 3: Native ASM `fp8x1_mul_redc_core_k2` Prototype

APIs:
- add symbol (internal): `blst_fp8x1_mul_redc_core_k2`.
- C++ wrapper uses asm path when native AVX-512 is enabled.

Data/layout:
- first prototype may materialize some intermediates if needed for correctness.
- prioritize functional equivalence and clear round structure.

Testing:
- strict differential test vs C++ reference for millions of random lanes.
- edge patterns: max limbs, carry-heavy patterns, alternating high/low limbs.

Success:
- bit-exact canonical outputs after normalize.
- stable benchmark execution on target AVX machine.

Avoid:
- no assumption that prototype perf is acceptable.

## Phase 4: Schedule/Live-Set Optimization (Primary Perf Phase)

APIs:
- no API change.

Data/layout:
- convert prototype into a no-spill or minimal-spill schedule.
- hand interleave stream operations where compiler scheduling is insufficient.
- tune round order to overlap IFMA with carry/shift/add work.

Testing:
- perf + correctness after each schedule change.
- asm inspection and perf counters (spills, uops mix, backend stalls).

Success:
- target throughput: at least `+20%` over current `k=2` `mul+redc` baseline.
- stretch target: `+30%` or better if spill-free schedule holds.
- no correctness regressions.

Avoid:
- avoid widening scope to `k=4` asm before `k=2` is stable.
- avoid macro-generated pseudo-interleaving that obscures true schedule.

## Phase 5: Integration and Cleanup

APIs:
- keep only one production fused path (`k=2`) and one fallback path.
- preserve existing `K=1` behavior.

Data/layout:
- `K=4` composition via two `k=2` calls.

Testing:
- full AVX portable-local regression build/tests.
- native AVX benchmark suite subset:
  - `BM_Fp8x1MulRedcK<1|2|4>`
  - `BM_Fp12Mul_(Scalar8|Pack8x1)` (impact check)
  - GT scalar mult benchmarks (indirect impact check).

Success:
- no regressions.
- reproducible benchmark improvement with low variance.

Avoid:
- no legacy benchmark variant sprawl after winner is selected.

## Benchmark Plan (Minimal and Stable)

Keep a small benchmark set:
- `BM_Fp8x1MulRedcK<1|2|4>` (primary)
- `BM_Fp8x1RedcK<1|2|4>` (secondary sanity)
- optional microbench:
  - fused core only
  - table/init-free loop-only decomposition if attribution is needed.

Recommended flags:
- `--benchmark_min_time=1s`
- `--benchmark_repetitions=10`
- `--benchmark_report_aggregates_only=true`

## Correctness Matrix

Required:
- random differential tests vs baseline path.
- adversarial limb patterns.
- cross-check through tower operations (`Fp2/Fp6/Fp12`) that consume Fp mul.
- invariant checks for lazy outputs before normalize (debug only).

## Promotion Criteria

Promote fused asm path as default `K==2` implementation only when:
- all correctness gates pass.
- benchmark gain is consistent across multiple runs on native AVX hardware.
- no unstable spill behavior across compiler versions/flags used by project CI.

## What To Avoid

- exposing this as a public API before k-wise support is stable.
- mixing algorithm changes and schedule changes in one patch (hard to attribute).
- retaining obsolete benchmark variants once a winner is established.
- optimizing AVX portable path as part of this phase.
