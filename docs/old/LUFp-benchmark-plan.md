# LUFp Benchmark Plan

## Goal

Build a coherent benchmark suite around the canonical `LUFp v2` path.

The suite must:

- measure the real accepted `v2` implementation, not legacy `LUFp` code;
- keep constituent accounting trustworthy;
- preserve the performance of the canonical `RawAdd`, `RawAddBias`, `RawAddBiasCarry`, `Full`, and `FullSign4` stages;
- support migration to canonical `v2` `Fp6` integration;
- make focused canonical runs valid and useful.

Accuracy comes first.

## Canonical Benchmark Principle

The benchmark suite is not required to preserve legacy `LUFp` paths.

For this plan:

- canonical `LUFp` means the `v2` implementation path rooted in `benchmark/internal/Fp_avx_lufp_v2_impl.h`;
- legacy `LUFpPack8x1`, `normalize<...>`, and `renorm<...>` benchmark paths are not part of the target suite;
- benchmarks that still depend on old slow `LUFp` code should be migrated or deleted.

## Performance Guardrail

The canonical add-family is the primary performance guardrail for the migration.

Accepted baseline for `4096 = 256 * 16` add operations:

- `RawAdd`: `12216 ns`
- `RawAddBias`: `19346 ns`
- `RawAddBiasCarry`: `49990 ns`
- `Full`: `81305 ns`
- `FullSign4`: `181537 ns`

Accepted stage deltas:

- bias: `1.74 ns`
- carry: `7.48 ns`
- 1-step repair: `7.65 ns`
- 4-step repair over `RawAddBiasCarry`: `32.12 ns`

Acceptance rule:

- the canonical add-family must remain within `3%` of this baseline;
- the 1-step and 4-step repair deltas are canaries and must also remain within `3%`.

## Hot-Path Rule

The canonical `v2` implementation must stay lean and compiler-friendly.

Do not add:

- abstraction layers on the canonical hot path;
- helper indirection that risks worse inlining or constant propagation;
- cleanups that make the measured code structurally heavier;
- sharing between canonical and experimental code if it makes canonical codegen worse.

The implementation should prefer:

- direct stage code;
- compile-time stage gating with `if constexpr`;
- visible hot helpers;
- minimal temporary materialization;
- minimal non-arithmetic overhead in measured paths.

## Benchmark Methodology

We use two benchmark styles.

### 1. Accounting benchmarks

These are the source of truth for constituent attribution.

Properties:

- one kernel family per operation;
- one compile-time stage enum per family;
- one templated implementation with `if constexpr` stage gating;
- stages are cumulative within a chosen path;
- deltas are derived from measured cumulative stages.

For families with alternative terminal stages, report those terminal deltas relative to their shared predecessor.

Example:

- canonical `LUFp v2` add uses:
  - `RawAdd`
  - `RawAddBias`
  - `RawAddBiasCarry`
  - `Full`
  - `FullSign4`

Interpretation:

- bias cost = `RawAddBias - RawAdd`
- carry cost = `RawAddBiasCarry - RawAddBias`
- 1-step repair cost = `Full - RawAddBiasCarry`
- 4-step repair cost = `FullSign4 - RawAddBiasCarry`

### 2. Pure kernel microbenchmarks

These are secondary.

They are useful for:

- isolated primitive checks;
- throughput or speed-of-light studies;
- validating whether a specific primitive changed.

They are not the source of truth for constituent attribution.

## Run Profiles

Reports and validation runs must be interpreted relative to an explicit profile.

### Canonical add-only

Required families:

- `UFp add accounting`
- `LUFp v2 add accounting`

Use this profile for:

- refactors that should not perturb the canonical add path;
- benchmark canary runs;
- quick remote validation.

### Canonical core

Required families:

- `UFp add accounting`
- `LUFp v2 add accounting`
- `LUFp v2 renorm accounting`

Use this profile for:

- validating the canonical exact-add and renorm story together.

### Canonical full

Required families:

- `UFp add accounting`
- `LUFp v2 add accounting`
- `LUFp v2 renorm accounting`
- `UFp Fp6 mul accounting`
- canonical `LUFp v2 Fp6 mul accounting`

This profile becomes the full target once canonical `v2` `Fp6` integration exists.

### Exploratory

Use this only when intentionally benchmarking:

- experimental code;
- partial migrations;
- one-off focused kernels.

Exploratory runs must be labeled as exploratory in the report.

## File Layout

Canonical and migration files:

- `benchmark/Fp_avx_ufp_bench.cpp`
- `benchmark/Fp_avx_lufp_v2_bench.cpp`
- `benchmark/Fp_avx_renorm_bench.cpp`
- `benchmark/internal/Fp_avx_lufp_v2_impl.h`

Rules:

- `Fp_avx_ufp_bench.cpp` contains the `UFp` benchmark families;
- `Fp_avx_lufp_v2_bench.cpp` contains the canonical `LUFp v2` benchmark families;
- `Fp_avx_renorm_bench.cpp` contains canonical `LUFp v2` renorm benchmarks;
- `Fp_avx_lufp_v2_impl.h` contains accepted canonical `v2` code;
- canonical `Fp6` `LUFp v2` support code should live in `Fp_avx_lufp_v2_impl.h`, not in a separate experimental benchmark header.

Legacy `LUFp` benchmark helpers should not continue to grow.

## Current Canonical Families

### UFp add

In `benchmark/Fp_avx_ufp_bench.cpp`:

- `BM_UFpAddAccounting<UFpAddStage::RawAdd>`
- `BM_UFpAddAccounting<UFpAddStage::RawAddBias>`
- `BM_UFpAddAccounting<UFpAddStage::RawAddBiasCarry>`
- `BM_UFpAddAccounting<UFpAddStage::Full>`

### LUFp v2 add

In `benchmark/Fp_avx_lufp_v2_bench.cpp`:

- `BM_LUFpAddAccounting<LUFpAddStage::RawAdd>`
- `BM_LUFpAddAccounting<LUFpAddStage::RawAddBias>`
- `BM_LUFpAddAccounting<LUFpAddStage::RawAddBiasCarry>`
- `BM_LUFpAddAccounting<LUFpAddStage::Full>`
- `BM_LUFpAddAccounting<LUFpAddStage::FullSign4>`

These are the canonical `LUFp v2` add benchmarks.

### LUFp v2 renorm

In `benchmark/Fp_avx_renorm_bench.cpp`:

- `BM_LUFpV2ExactAddRenormAccounting<LUFpV2ExactAddRenormStage::Bias>`
- `BM_LUFpV2ExactAddRenormAccounting<LUFpV2ExactAddRenormStage::BiasCarry>`
- `BM_LUFpV2ExactAddRenormAccounting<LUFpV2ExactAddRenormStage::Full>`
- `BM_LUFpV2ExactAddRenormAccounting<LUFpV2ExactAddRenormStage::FullSign4>`

These benchmarks measure the same `v2` exact-add renorm primitives used by the canonical add-family path.

### UFp Fp6 mul

In `benchmark/Fp_avx_ufp_bench.cpp`:

- `BM_Fp6MulUfpAccounting<UFpFp6MulStage::PreAdd>`
- `BM_Fp6MulUfpAccounting<UFpFp6MulStage::PreAddMul>`
- `BM_Fp6MulUfpAccounting<UFpFp6MulStage::Full>`

This is the current comparison baseline for `Fp6`.

### LUFp v2 Fp6 mul

In `benchmark/Fp_avx_lufp_v2_bench.cpp`:

- `BM_Fp6MulLufpAccounting<LUFpFp6MulStage::PreAdd>`
- `BM_Fp6MulLufpAccounting<LUFpFp6MulStage::PreAddMul>`
- `BM_Fp6MulLufpAccounting<LUFpFp6MulStage::PreAddMulRawCombine>`
- `BM_Fp6MulLufpAccounting<LUFpFp6MulStage::Full>`

This family is intended to be the canonical `v2` `Fp6` benchmark path end to end.

## Remaining Migration Gaps

The suite is now aligned around canonical `v2` for add and renorm.

The remaining gaps are:

- `benchmark/Fp12_tower_avx_bench.cpp`
  - still benchmarks legacy `LUFpPack8x1` tower paths and should be migrated or deleted
- `benchmark/internal/Fp_avx_bench_kernels.h`
  - still contains legacy `LUFp` benchmark helpers that should shrink as the old path is removed

## Desired End State

The desired suite is:

- canonical `UFp` add accounting;
- canonical `LUFp v2` add accounting;
- canonical `LUFp v2` renorm accounting;
- canonical `UFp` `Fp6` accounting;
- canonical `LUFp v2` `Fp6` accounting.

Legacy `LUFp` tower or normalize benchmarks should either:

- migrate to canonical `v2`;
- or be deleted.

## Accuracy Rules

These rules apply to every benchmark family.

1. Accounting variants must keep the same sink shape.

2. Compared variants must keep comparable storage behavior.

3. Compared hot kernels must stay on equal footing for inlining and visibility.

4. Canonical-path refactors must be benchmarked, not assumed safe.

5. Benchmark only on an idle server.

6. Do not build and benchmark at the same time.

7. If a gap is unexplained, say so plainly.

## Removal Rule

Delete benchmark code when all of the following are true:

- it measures legacy `LUFp` behavior;
- it is not needed for the canonical `v2` target;
- it does not support the current benchmark methodology.

There is no requirement to keep obsolete slow benchmark paths once the canonical `v2` replacement exists.

## Execution Order

1. Keep the canonical `LUFp v2` add-family stable.
2. Update docs and reporting around the real canonical suite.
3. Replace legacy renorm benchmarks with canonical `v2` renorm benchmarks.
4. Keep `LUFp v2` `Fp6` integration in canonical `impl` code and benchmark it directly.
5. Remove legacy `LUFp` benchmark code that no longer serves the target suite.

## Decision Rule

When microbenchmarks and accounting benchmarks disagree:

- accounting benchmarks win for constituent attribution.

When a speed-of-light estimate and measured timing disagree:

- measured timing is the fact;
- the gap is an investigation item;
- unsupported causal claims are not allowed.
