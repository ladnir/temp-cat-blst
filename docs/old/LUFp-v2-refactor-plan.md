# LUFp v2 Refactor Plan

## Goal

Reorganize the `LUFp v2` benchmark code so there is:

- one canonical implementation file
- one experimental implementation file
- benchmark harness files that call into those implementations

This step is organizational only.

## Non-Goals

This step must not:

- change behavior
- refactor logic
- simplify code
- rename algorithms
- change benchmark methodology
- change benchmark sinks
- change benchmark data generation

No "while we're here" changes are allowed.

## File Roles

### Canonical implementation

File:

- `benchmark/internal/Fp_avx_lufp_v2_impl.h`

This file contains the active implementation used by the default `LUFp v2` benchmark path.

Anything used by the accepted add-family benchmark path belongs here.

### Experimental implementation

File:

- `benchmark/internal/Fp_avx_lufp_v2_experimental.h`

This file contains any non-canonical alternative implementations.

Anything that is not the active default path but is still worth preserving for comparison belongs here.

### Benchmark harness

Files:

- `benchmark/Fp_avx_lufp_v2_bench.cpp`
- any other benchmark `.cpp` files that invoke canonical or experimental routines

These files contain:

- benchmark fixtures
- data generation
- sinks
- stage enums
- benchmark registration

They must not hide algorithm choices.

## Classification Rule

For this refactor, code is classified exactly as it exists today.

- If the current default benchmark path uses it, it is canonical.
- If it is an alternative or old path, it is experimental.
- If it is only benchmark plumbing, it stays in the benchmark harness.

No code changes are allowed during classification.

## Explicit Reachability Rule

For this refactor, the canonical set is defined operationally as:

- every type, helper, and routine reachable from these benchmark entries exactly as they exist today:
  - `BM_LUFpAddAccounting<LUFpAddStage::RawAdd>`
  - `BM_LUFpAddAccounting<LUFpAddStage::RawAddBias>`
  - `BM_LUFpAddAccounting<LUFpAddStage::RawAddBiasCarry>`
  - `BM_LUFpAddAccounting<LUFpAddStage::Full>`
  - `BM_LUFpAddAccounting<LUFpAddStage::FullSign4>`

Everything not reachable from that set is non-canonical for the purpose of this refactor.

This means:

- canonical = the currently accepted/default add-family implementation path
- experimental = alternative or older implementations not on that path
- harness = fixtures, enums, data generation, sinks, benchmark wrappers, and registration

## FP6 Scope Clarification

`Fp6`-mul-related code may remain in the current implementation header without becoming canonical for this step.

In particular, any helper used only by:

- `BM_Fp6MulLufpWide`
- `BM_Fp6MulLufpAccounting<...>`

is out of scope for the add-family acceptance baseline in this document.

For this refactor step:

- add-family reachability determines canonical status
- `Fp6` paths that are not on the add-family reachability graph should be treated as experimental unless and until a separate `Fp6` baseline is defined

## Preservation Rule

The phrase “worth preserving for comparison” is made explicit as follows:

- if a routine is not on the canonical add-family path but is still referenced by an active benchmark, it moves to `experimental`
- if a routine is not on the canonical add-family path and is not referenced by any active benchmark, it should not be moved during this step unless there is an explicit cleanup task

This refactor is classification-only.

- no deletion decisions are made here
- no behavior changes are made here

## Canonical Before/After Performance Baseline

The accepted pre-refactor baseline is the stable add-family run on the idle AWS server.

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

Stage deltas:

- bias: `1.74 ns`
- carry: `7.48 ns`
- 1-step repair: `7.65 ns`
- 4-step repair over `RawAddBiasCarry`: `32.12 ns`

## Performance Acceptance Criteria

This refactor passes only if the post-refactor benchmark numbers remain within `3%` of the accepted baseline above.

This requirement applies to:

- `BM_LUFpAddAccounting<LUFpAddStage::RawAdd>`
- `BM_LUFpAddAccounting<LUFpAddStage::RawAddBias>`
- `BM_LUFpAddAccounting<LUFpAddStage::RawAddBiasCarry>`
- `BM_LUFpAddAccounting<LUFpAddStage::Full>`
- `BM_LUFpAddAccounting<LUFpAddStage::FullSign4>`

The critical canary metrics are:

- 1-step repair delta: `7.65 ns`
- 4-step repair delta: `32.12 ns`

If either of these moves by more than `3%`, the refactor fails.

## Validation Procedure

After the move-only refactor:

1. rebuild `blst_benchmarks`
2. rerun:
   - `BM_LUFpAddAccounting.*`
3. compare against the baseline values in this document
4. stop if any regression exceeds `3%`

## Sequencing

1. Move canonical code into `Fp_avx_lufp_v2_impl.h`
2. Move non-canonical code into `Fp_avx_lufp_v2_experimental.h`
3. Leave benchmark harness logic in the benchmark `.cpp`
4. Rebuild
5. Rerun the add-family benchmarks
6. Verify the acceptance criteria
