# LUFp Benchmark Plan

## Goal

Build a coherent benchmark suite for `UFp`, `LUFp`, and `renorm` that:

- keeps benchmark code easy to read and maintain;
- separates full end-to-end timings from constituent-part timings;
- provides stage timings whose deltas add up to the full function timing;
- avoids misleading results caused by mismatched sinks, storage, or benchmark structure.

Accuracy is the top priority.

## Speed-of-Light Methodology

The benchmark suite should be interpreted with a speed-of-light mindset.

For each important kernel family:

- start from first principles;
- estimate a lower bound or roofline;
- compare the measured result against that bound;
- treat the gap as an investigation target.

Rules:

- state assumptions explicitly;
- distinguish facts from inference;
- distinguish inference from speculation;
- do not explain a gap without evidence;
- be open-minded about why the measured result is above the lower bound.

This methodology complements the benchmark suite:

- accounting benchmarks tell us where measured time goes;
- speed-of-light reasoning tells us how far we are from the best plausible outcome.

## Benchmark Types

We will use two distinct benchmark styles.

### 1. Accounting benchmarks

These are the primary benchmarks for attribution.

Properties:

- one kernel family per operation;
- one compile-time stage enum per family;
- one templated implementation with `if constexpr` gating each stage;
- stages are cumulative;
- stage deltas must add up to the full benchmark by construction.

Example:

- `RawAdd`
- `RawAddPlusBias`
- `RawAddPlusBiasPlusCarry`
- `Full`

Then:

- `bias cost = stage(Bias) - stage(RawAdd)`
- `carry cost = stage(Carry) - stage(Bias)`
- `repair cost = stage(Full) - stage(Carry)`

This is the only benchmark style that should be used for "what part of the full function costs what?"

### 2. Pure kernel microbenchmarks

These are secondary.

Properties:

- isolate one kernel such as `add`, `mul`, or `renorm`;
- useful for ceiling/throughput studies;
- useful for checking whether a specific primitive changed;
- do not need to sum exactly to the full function.

These are useful, but they are not the source of truth for constituent accounting.

They are, however, useful supporting inputs for speed-of-light estimates.

## File Layout

Planned files:

- `benchmark/Fp_avx_ufp_bench.cpp`
- `benchmark/Fp_avx_lufp_bench.cpp`
- `benchmark/Fp_avx_renorm_bench.cpp`
- `benchmark/internal/Fp_avx_bench_kernels.h`

Rules:

- `Fp_avx_ufp_bench.cpp` contains only `UFp` benchmark families.
- `Fp_avx_lufp_bench.cpp` contains only `LUFp` benchmark families.
- `Fp_avx_renorm_bench.cpp` contains only `renorm` benchmark families.
- shared setup, fixtures, sinks, and stage-gated helper kernels live in `benchmark/internal/Fp_avx_bench_kernels.h`.

Existing ad hoc benchmark experiments should be removed or moved once the new suite replaces them.

## Coding Style For Accounting Benchmarks

Each accounting family should use:

- one stage enum;
- one templated kernel implementation;
- one templated benchmark entry point.

Pattern:

```cpp
enum class UFpAddStage
{
    RawAdd,
    RawAddBias,
    RawAddBiasCarry,
    Full,
};

template <UFpAddStage Stage>
BLST_FORCEINLINE void ufp_add_accounting(
    UFpPack8x1& r,
    const UFpPack8x1& a,
    const UFpPack8x1& b)
{
    // Always do the raw add.

    if constexpr (Stage >= UFpAddStage::RawAddBias) {
        // Bias stage.
    }
    if constexpr (Stage >= UFpAddStage::RawAddBiasCarry) {
        // Carry stage.
    }
    if constexpr (Stage >= UFpAddStage::Full) {
        // Final repair stage.
    }
}
```

This keeps the implementation coherent and prevents drift between variants.

## Benchmark Families

### UFp

In `benchmark/Fp_avx_ufp_bench.cpp`:

#### A. UFp add

- `BM_UFpAdd_Full`
- `BM_UFpAdd_Accounting<Stage>`

Planned stages:

- `RawAdd`
- `RawAddBias`
- `RawAddBiasCarry`
- `Full`

#### B. UFp add microbenchmarks

These are secondary and should correspond to real sub-kernels:

- raw limb add;
- bias-only;
- carry-only;
- repair-only if it exists as a meaningful isolated piece.

The primary attribution still comes from the accounting family above.

The `UFp add` family should support a speed-of-light summary that reports:

- lower-bound estimate;
- measured full timing;
- ratio of measured to lower bound;
- what is known vs what is inferred.

#### C. Fp6 x UFp

- `BM_Fp6Mul_UFp_Full`
- `BM_Fp6Mul_UFp_Accounting<Stage>`

Planned stages:

- `PreAdd`
- `PreAddMul`
- `PreAddMulWideCombine`
- `Full`

From these and the known op counts we derive:

- average time per `Fp2` multiply;
- average time per wide combine add/sub;
- average time for the remaining stage.

These stage deltas must sum to the full `Fp6` timing.

The `Fp6 x UFp` family should also support a speed-of-light summary that reports:

- lower-bound estimate from topology and kernel structure;
- measured full timing;
- measured stage timings;
- the remaining unexplained gap, if any.

### LUFp

In `benchmark/Fp_avx_lufp_bench.cpp`:

#### A. LUFp add

- `BM_LUFpAdd_Full`
- `BM_LUFpAdd_Accounting<Stage>`

Planned stages:

- `RawAdd`
- `RawAddRenorm`

If the `LUFp` add family later needs more substructure, extend the enum rather than adding separate near-duplicate functions.

#### B. LUFp add microbenchmarks

Secondary:

- raw limb add;
- raw limb sub;
- optional exact renorm base case;
- optional generalized renorm classes.

These microbenchmarks are also useful for establishing lower-bound expectations for `LUFp`.

#### C. Fp6 x LUFp

- `BM_Fp6Mul_LUFp_Full`
- `BM_Fp6Mul_LUFp_Accounting<Stage>`

Planned stages:

- `PreAdd`
- `PreAddMul`
- `PreAddMulRawLazyCombine`
- `PreAddMulRawLazyCombineRenorm`
- `Full`

This isolates:

- multiply work;
- raw lazy combine work;
- renorm work;
- any remaining full-pipeline work.

Again, these stages must add up to the full benchmark.

The `Fp6 x LUFp` family should make it possible to report:

- expected advantage from cheaper raw lazy add/sub;
- expected renorm cost;
- measured full timing;
- the gap between expectation and measurement;
- a clear split between proven causes and hypotheses.

### Renorm

In `benchmark/Fp_avx_renorm_bench.cpp`:

- `BM_Renorm_<Class>_Full`
- `BM_Renorm_<Class>_Accounting<Stage>`

Benchmarks should be defined by exact classes, not vague names.

Examples:

- `ExactAdd`
- `Fp6MulCommon`

Planned stages:

- `Bias`
- `BiasCarry`
- `BiasCarryRepair1`
- `...`
- `Full`

This gives exact cost of:

- bias;
- carry;
- each repair stage.

It should also support a speed-of-light analysis for each renorm class:

- lower-bound estimate from bias + carry + repair structure;
- measured timing;
- unexplained gap, if any.

## Accuracy Rules

These rules apply to every benchmark family.

1. Accounting variants must have the same sink shape.

- same output storage style;
- same temporaries when possible;
- same benchmark loop shape.

2. Do not compare variants with mismatched storage behavior.

- Example: do not compare one benchmark that writes multiple large lazy vectors per iteration to one benchmark that writes one strict output and then claim the difference is arithmetic cost.

3. Keep compared hot kernels on equal footing for inlining.

- If one compared kernel is visible in the translation unit, the other should be too.

4. Benchmark only on an idle server.

5. Do not build and benchmark at the same time.

6. For close results, collect both:

- benchmark timing;
- `perf stat` style counters such as cycles and instructions, if needed later.

7. Every important family should support a speed-of-light summary.

At minimum, it should be possible to report:

- lower-bound estimate;
- measured timing;
- gap;
- known causes of the gap;
- unknown causes of the gap.

8. If a gap is not explained, say so plainly.

Do not convert an unexplained gap into an implicit explanation.

## What To Remove Or Move

The current ad hoc benchmark experiments should not remain spread across unrelated files once the new suite is in place.

Planned cleanup:

- move reusable helpers into `benchmark/internal/Fp_avx_bench_kernels.h`;
- move `UFp`-specific decomposition into `benchmark/Fp_avx_ufp_bench.cpp`;
- move `LUFp`-specific decomposition into `benchmark/Fp_avx_lufp_bench.cpp`;
- move `renorm` decomposition into `benchmark/Fp_avx_renorm_bench.cpp`;
- remove obsolete experimental benchmark paths after replacement.

## Execution Order

1. Implement `UFp add` accounting family.
2. Implement `LUFp add` accounting family.
3. Implement `renorm` accounting family.
4. Implement `Fp6 x UFp` accounting family.
5. Implement `Fp6 x LUFp` accounting family.
6. Remove superseded ad hoc benchmark code.

## Decision Rule

When there is a disagreement between:

- a pure kernel microbenchmark;
- and an accounting benchmark;

the accounting benchmark wins for constituent attribution.

When a speed-of-light estimate and measured timing disagree:

- the measured timing is the fact;
- the gap becomes an investigation item;
- any explanation of that gap must be marked as inference or speculation unless directly supported by evidence.
