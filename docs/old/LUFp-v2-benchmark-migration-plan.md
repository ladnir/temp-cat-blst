# LUFp v2 Benchmark Migration Plan

## Goal

Move the benchmark suite to a single `canonical v2` `LUFp` implementation path and remove the old slow `LUFp` benchmark code.

The desired end state is:

- `LUFp` benchmarks use the `v2` implementation path all the way through;
- `UFp` benchmarks remain the comparison baseline;
- any benchmark that still depends on the old `LUFpPack8x1` / legacy normalize path is either migrated or deleted;
- benchmark and report docs describe the actual canonical suite, not the superseded plan.

This is a migration plan, not a request to preserve old benchmark paths indefinitely.

## Non-Goals

This plan does not require:

- keeping legacy `LUFp` benchmark paths for historical comparison;
- preserving benchmarks that do not support the desired accounting methodology;
- preserving old helper code once the canonical `v2` path replaces it;
- changing the `UFp` benchmark methodology unless needed for parity with `LUFp v2`.

This plan also does not permit:

- adding abstractions, wrappers, or helper layers on the canonical hot path if they risk worse code generation;
- refactoring canonical hot code for readability at the expense of compiler optimization friendliness;
- "cleanup" changes in canonical hot code unless they are shown not to perturb the canonical benchmark path.

## Canonical Definition

For this migration, `canonical v2` means:

- the implementation rooted in `benchmark/internal/Fp_avx_lufp_v2_impl.h`;
- the accepted add-family accounting path in `benchmark/Fp_avx_lufp_v2_bench.cpp`;
- future higher-level `LUFp` integration benchmarks that are built from the same `v2` primitives and not from legacy `LUFpPack8x1` normalize/renorm code.

Anything built on:

- `LUFpPack8x1`;
- `LUFpPack8x1::normalize<...>`;
- `LUFpPack8x1::renorm<...>`;
- legacy `LUFpNormalizeBound` / `LUFpRenormClass` paths;

is non-canonical and should be considered migration or deletion work.

## Canonical Performance Guardrail

The canonical `LUFp v2` add-family path is performance-critical.

For this migration, the accepted baseline is the canonical add-family baseline from `docs/LUFp-v2-refactor-plan.md`.

Totals for `4096 = 256 * 16` add operations:

- `RawAdd`: `12216 ns`
- `RawAddBias`: `19346 ns`
- `RawAddBiasCarry`: `49990 ns`
- `Full`: `81305 ns`
- `FullSign4`: `181537 ns`

Per-add values:

- `RawAdd`: `2.98 ns`
- `RawAddBias`: `4.72 ns`
- `RawAddBiasCarry`: `12.20 ns`
- `Full`: `19.85 ns`
- `FullSign4`: `44.32 ns`

Canonical stage deltas:

- bias: `1.74 ns`
- carry: `7.48 ns`
- 1-step repair: `7.65 ns`
- 4-step repair over `RawAddBiasCarry`: `32.12 ns`

## Performance Target

The migration succeeds only if the canonical add-family remains within `3%` of the accepted baseline above.

This applies to:

- `BM_LUFpAddAccounting<LUFpAddStage::RawAdd>`
- `BM_LUFpAddAccounting<LUFpAddStage::RawAddBias>`
- `BM_LUFpAddAccounting<LUFpAddStage::RawAddBiasCarry>`
- `BM_LUFpAddAccounting<LUFpAddStage::Full>`
- `BM_LUFpAddAccounting<LUFpAddStage::FullSign4>`

Critical canaries:

- 1-step repair delta: `7.65 ns`
- 4-step repair delta: `32.12 ns`

If either critical delta moves by more than `3%`, the change fails even if the absolute full-stage timing still looks acceptable.

## Canonical Hot-Path Rules

The canonical add-family path must stay lean.

For the canonical hot path:

- `RawAdd`, `RawAddBias`, `RawAddBiasCarry`, `Full`, and `FullSign4` should only pay for their intended arithmetic stages plus small unavoidable benchmark overhead;
- no extra abstraction layers should be introduced between benchmark entry and hot arithmetic code unless they are proven codegen-neutral;
- no helper extraction should be done purely for cleanup if it risks blocking inlining, constant propagation, dead-code elimination, or register allocation quality;
- no wrapper types should be introduced on the canonical path unless they clearly compile away;
- no "shared" refactor should be accepted if it makes the canonical path structurally more indirect just to unify canonical and non-canonical code.

Compiler optimization friendliness is a requirement, not a nice-to-have.

The implementation should prefer:

- direct stage code;
- compile-time stage gating with `if constexpr`;
- visibility of hot helpers to the compiler;
- minimal control-flow and storage overhead;
- minimal temporary materialization beyond what the arithmetic actually needs.

The implementation should avoid:

- virtual or runtime-dispatched hot-path choices;
- unnecessary wrapper calls in the canonical path;
- generalization that hides constants from the compiler;
- benchmark-harness abstractions that become part of the measured hot path;
- "clean" structural symmetry between canonical and experimental code when that symmetry makes the canonical code worse.

## Code Review Rule

When changing canonical `v2` code, performance preservation wins over stylistic cleanup.

If a proposed change touches:

- `benchmark/internal/Fp_avx_lufp_v2_impl.h`
- the canonical add-family implementation in `benchmark/Fp_avx_lufp_v2_bench.cpp`
- any future canonical `v2` `Fp6` integration code

the default review posture must be:

- assume codegen can change;
- assume benchmark deltas matter;
- reject unnecessary abstraction unless there is a concrete reason to keep it;
- require benchmark confirmation for non-trivial structural edits.

## Current State

### Already canonical

- `benchmark/Fp_avx_lufp_v2_bench.cpp`
  - `BM_LUFpAddAccounting<LUFpAddStage::RawAdd>`
  - `BM_LUFpAddAccounting<LUFpAddStage::RawAddBias>`
  - `BM_LUFpAddAccounting<LUFpAddStage::RawAddBiasCarry>`
  - `BM_LUFpAddAccounting<LUFpAddStage::Full>`
  - `BM_LUFpAddAccounting<LUFpAddStage::FullSign4>`
- `benchmark/internal/Fp_avx_lufp_v2_impl.h`
  - canonical `v2` add/bias/carry/one-step repair primitives
  - canonical exact-add renorm specialization
- `benchmark/Fp_avx_renorm_bench.cpp`
  - `BM_LUFpV2ExactAddRenormAccounting<LUFpV2ExactAddRenormStage::Bias>`
  - `BM_LUFpV2ExactAddRenormAccounting<LUFpV2ExactAddRenormStage::BiasCarry>`
  - `BM_LUFpV2ExactAddRenormAccounting<LUFpV2ExactAddRenormStage::Full>`
  - `BM_LUFpV2ExactAddRenormAccounting<LUFpV2ExactAddRenormStage::FullSign4>`
  - benchmarks the same `v2` exact-add renorm primitives used by the canonical add path

### Structurally good comparison baseline

- `benchmark/Fp_avx_ufp_bench.cpp`
  - `UFp` add accounting family is already cumulative and matches the intended methodology
  - `UFp` `Fp6` accounting family is already a coherent integration benchmark

### Still legacy or mixed

- `benchmark/internal/Fp_avx_bench_kernels.h`
  - still contains legacy `LUFpPack8x1` add/sub/normalize/integration helpers
- `benchmark/Fp12_tower_avx_bench.cpp`
  - still uses legacy `LUFpPack8x1`, `normalize<HiNeg1Pos2>`, and `LUFpRenormClass::Loose`
  - does not represent canonical `v2`
- `benchmark/Fp_avx_lufp_v2_bench.cpp`
  - `Fp6 x LUFp` accounting path exists and should now be treated as the canonical `v2` benchmark path

## Progress After This Update

Completed in this migration step:

- `docs/LUFp-benchmark-plan.md`
  - updated to describe the real canonical `v2` suite and the hot-path performance rule
- `docs/LUFp-benchmark-report-spec.md`
  - updated to use run profiles and focused canonical runs instead of requiring every family in every report
- `scripts/lufp_benchmark_report.py`
  - updated to understand canonical `v2` add and renorm accounting families
  - enforces canonical add-family guardrails on non-exploratory runs
- `benchmark/Fp_avx_renorm_bench.cpp`
  - migrated from legacy normalize/prefix-lex benchmarking to canonical `v2` exact-add renorm accounting
- `benchmark/Fp_avx_lufp_v2_bench.cpp`
  - no longer depends on `benchmark/internal/Fp_avx_lufp_v2_experimental.h` for the `Fp6` benchmark path
- `benchmark/internal/Fp_avx_lufp_v2_impl.h`
  - now contains the accepted `Fp6` `LUFp v2` support code used by the benchmark path

Remaining implementation work is now concentrated in:

- deletion or migration of legacy `LUFp` tower benchmarks;
- deletion or quarantine of legacy `LUFpPack8x1` benchmark helpers once no longer needed.

## Desired End State

The benchmark suite should converge to the following structure.

### 1. Canonical `LUFp v2` add benchmarks

Keep the current `v2` add-family accounting suite as the canonical source of truth:

- `RawAdd`
- `RawAddBias`
- `RawAddBiasCarry`
- `Full`
- `FullSign4`

Interpretation:

- `Full` is the canonical exact-add finish path
- `FullSign4` is an alternative staged repair path that remains useful as a benchmarked comparison/canary

These stages are also the canonical performance guardrail for the entire migration.

Migration work must not degrade them in the name of:

- cleanup;
- sharing code with experimental paths;
- making the code more generic;
- moving logic behind extra helper layers.

### 2. Canonical `LUFp v2` renorm benchmarks

Replace legacy renorm benchmarks with `v2` renorm benchmarks built from `bench_lufp_v2` code.

The new renorm micro/accounting benchmarks should measure canonical `v2` pieces such as:

- bias-to-biased form
- carry from biased form
- one-step repair from sign
- staged sign-step repair variants where still useful

These should benchmark the same `v2` routines used by the canonical add path, not legacy normalize logic.

### 3. Canonical `LUFp v2` `Fp6` integration benchmarks

Keep the canonical `v2` `Fp6` integration path in the main `v2` implementation file and benchmark it directly.

This means:

- raw lazy combine may stay as an intermediate accounting stage;
- the final `Fp6` output path must be implemented using canonical `v2` primitives;
- `BM_Fp6MulLufpWide` and `BM_Fp6MulLufpAccounting<...>` should ultimately benchmark the canonical `v2` path end to end.

### 4. No legacy `LUFp` tower benchmarks in the canonical suite

`benchmark/Fp12_tower_avx_bench.cpp` should not continue to benchmark old slow `LUFp` code as if it were the current `LUFp` path.

For this migration:

- migrate tower-level `LUFp` benchmarks to canonical `v2` if they are still needed;
- otherwise delete the legacy `LUFp`-specific tower benchmarks;
- do not keep legacy `LUFp` tower code just because it exists today.

## Files To Change

### Must update

- `benchmark/internal/Fp_avx_bench_kernels.h`
- `benchmark/Fp12_tower_avx_bench.cpp`

### Already updated in this step

- `docs/LUFp-benchmark-plan.md`
- `docs/LUFp-benchmark-report-spec.md`
- `scripts/lufp_benchmark_report.py`
- `benchmark/Fp_avx_renorm_bench.cpp`
- `benchmark/Fp_avx_lufp_v2_bench.cpp`
- `benchmark/internal/Fp_avx_lufp_v2_impl.h`

### May add

- additional `benchmark/internal/Fp_avx_lufp_v2_*.h` headers if needed to keep canonical vs experimental code explicit during migration
- a new canonical `v2` helper header for shared `Fp6` integration code if the existing split becomes too crowded

### Should shrink or disappear

- legacy `LUFp` helper code in `benchmark/internal/Fp_avx_bench_kernels.h`
- legacy `LUFp` benchmark bodies in `benchmark/Fp12_tower_avx_bench.cpp`
- report expectations that treat legacy `LUFp` families as required forever

## Migration Phases

### Phase 1. Update docs to match reality

Rewrite the benchmark and report docs around the current real suite.

Required changes:

- define `Fp_avx_lufp_v2_bench.cpp` as the `LUFp` benchmark home
- define the canonical add-family stages exactly as they exist today
- distinguish canonical runs from exploratory/legacy runs
- stop requiring legacy `LUFp` families in every report
- define the target state as `canonical v2 all the way through`

Deliverable:

- docs no longer describe superseded benchmark names or old `LUFp` structure

### Phase 2. Replace legacy renorm benchmarks with `v2` renorm benchmarks

`benchmark/Fp_avx_renorm_bench.cpp` should be rewritten to benchmark `bench_lufp_v2` primitives.

Target benchmark families:

- canonical exact-add renorm accounting
- optional staged repair microbenchmarks for comparison

Delete:

- legacy `normalize<HiNeg1Pos2>` benchmark entries
- legacy prefix-lex repair benchmarks if they are not part of the `v2` story

Acceptance:

- renorm benchmark names and stages align with `v2` add-path terminology

### Phase 3. Make `Fp6 x LUFp` canonical in `v2`

Keep `LUFp v2` `Fp6` integration in canonical `v2` implementation code.

Required work:

- implement the `Full` `Fp6` path using canonical `v2` code
- keep accepted `Fp6` helpers in canonical `impl` code
- avoid reintroducing a separate experimental benchmark header for accepted `Fp6` support

Constraint:

- this phase must not add overhead to the canonical add-family path in order to make `Fp6` integration cleaner or more uniform

Acceptance:

- `BM_Fp6MulLufpAccounting<...>` `Full` benchmarks canonical code only
- `BM_Fp6MulLufpWide` benchmarks the same canonical `v2` full path

### Phase 4. Remove legacy `LUFp` tower benchmarks

Audit `benchmark/Fp12_tower_avx_bench.cpp`.

For each `LUFp`-specific benchmark:

- migrate it to canonical `v2`, or
- delete it if it no longer serves the benchmark plan

Do not preserve:

- old slow `LUFp` experiments that are neither canonical nor needed for current analysis

Acceptance:

- no benchmark remains in the canonical suite solely because it was historically useful

### Phase 5. Simplify report generation

Once the benchmark migration lands:

- update the report script to expect the canonical suite
- remove legacy-family assumptions
- report canonical `LUFp v2` add, renorm, and `Fp6` families explicitly
- mark any remaining experimental families as optional, not required

## Deletion Rule

During this migration, delete benchmark code when all of the following are true:

- it benchmarks legacy `LUFp` behavior;
- it is not part of the canonical `v2` target;
- it does not provide required comparison value for the current benchmark methodology.

There is no requirement to keep old slow benchmark code once the canonical `v2` replacement exists.

There is also no requirement to preserve a refactor if it makes the canonical `v2` path slower.

## Acceptance Criteria

This migration is complete only when:

1. the benchmark docs describe the current canonical `v2` suite accurately;
2. canonical `LUFp` add benchmarks use only `v2` code;
3. canonical `LUFp` renorm benchmarks use only `v2` code;
4. canonical `LUFp` `Fp6` integration benchmarks use only `v2` code;
5. legacy `LUFpPack8x1` benchmark paths are removed from the canonical benchmark suite;
6. any remaining experimental `LUFp` benchmarks are clearly marked experimental and are not required by the report spec;
7. the canonical add-family remains within `3%` of the accepted baseline;
8. the canonical raw add, bias, carry, and repair stages remain attributable to the intended arithmetic work rather than avoidable structural overhead.

## Immediate Next Steps

1. Rewrite `docs/LUFp-benchmark-plan.md` around the real `v2` suite and this target state.
2. Rewrite `docs/LUFp-benchmark-report-spec.md` so focused canonical runs are valid.
3. Replace `benchmark/Fp_avx_renorm_bench.cpp` with `v2` renorm benchmark families.
4. Define the accepted canonical `v2` `Fp6` full path and migrate `BM_Fp6MulLufpAccounting` to it.
5. Audit and delete or port legacy `LUFp` paths in `benchmark/Fp12_tower_avx_bench.cpp`.
