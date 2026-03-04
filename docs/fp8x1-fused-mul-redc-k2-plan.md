# Fp8x1 Fused mul+redc Kernel Plan (K=1 -> K=2)

## Goal

Optimize the fused `Fp8x1` `mul+redc` kernel for throughput.

Primary objective:
- maximize arithmetic throughput while keeping register pressure low enough to scale from `k=1` to `k=2` interleaving.

This plan is kernel-focused only.

## Scope (Strict)

In scope:
- benchmark-only kernel development and evaluation.
- native AVX-512 (`BLST_AVX512`, non-portable), Linux x86_64 target.
- `k=1` kernel first, then `k=2` interleaved kernel.

Out of scope:
- library/public API changes.
- dispatch/integration wiring.
- AVX portable optimization.
- GT/G1/G2 integration.

## Algorithm Choice

Use an integrated Montgomery **CIOS-style** kernel with in-place accumulator flow.

Why this choice:
- lower live state than full-product (`z0..z15`) style.
- better fit for `k=2` without spills.
- clear round structure for hand scheduling and asm mapping.

Not the initial target:
- full-FIOS with aggressive per-term interleaving.
- staged `mul_v2` then `redc_core`.

## Kernel Dataflow (Per Stream)

For `n=8` limbs:
- maintain rolling accumulator `t[0..8]` (9 limbs).
- each outer round `i`:
1. accumulate `a_i * b_j` terms into `t`.
2. compute Montgomery factor `m = t0 * n0' mod 2^48`.
3. fold `m * p_j` into `t`.
4. shift/advance accumulator window logically (ring index, not bulk moves).

Boundary:
- keep output lazy (`LFpPack8x1`) in core benchmark.
- measure normalize separately and combined.

## Register-Pressure Strategy

Design constraints for `k=1` implementation:
- avoid materializing full widened product arrays.
- avoid duplicate temporaries (`w*`, `r*`) after reduction.
- keep constants minimal in live registers:
  - pin only always-used constants (`mask`, `montw`, `zero`).
  - load/use `p[j]` with short live range.
- use ring-indexed accumulator to avoid shuffle/move storms.

Target for `k=2`:
- replicate same round machine for stream A/B.
- interleave at round/sub-round granularity to hide IFMA latency.
- keep stream state independent and bounded to avoid spills.

## Implementation Phases

## Phase 1: K=1 CIOS Baseline (Benchmark-Only)

Deliverables:
- benchmark-only fused `k=1` CIOS kernel.
- A/B benchmarks:
  - current `k=1` (`fp8x1_mul_redc_k<1>`)
  - fused core only (`lazy`)
  - fused core + normalize

Correctness:
- deterministic differential check inside benchmark path against current kernel output.

Success:
- correctness parity.
- measurable gain vs current `k=1`.

## Phase 2: K=1 Schedule Tuning

Focus:
- reduce live ranges and redundant register copies.
- tune instruction ordering inside rounds.
- keep constant pressure low.
- use objdump/perf evidence to drive changes (not intuition).

### Current Linux Findings (NoInline C++ vs K1 ASM)

Measured on native AVX-512 Linux with:
- `BM_Fp8x1MulRedcCoreK1_FusedBaseline` (C++ noinline core)
- `BM_Fp8x1MulRedcCoreK1_FusedAsm` (hand ASM core)

Throughput snapshot (latest A/B):
- current non-fused baseline (`BM_Fp8x1MulRedcK1_Current`): `~12.71M/s`
- C++ noinline core: `~16.24M/s`
- C++ noinline core + normalize: `~13.74M/s`
- ASM core: `~15.16M/s`
- ASM core + normalize: `~12.48M/s`

Objdump summary:
- both kernels execute the same IFMA count (`232 x vpmadd52*`).
- C++ noinline core: `0` memory-operand `vpmadd52`.
- ASM core: `72` memory-operand `vpmadd52` (reduced from earlier `128`).
- C++ uses stack scratch (`[rsp]`) to keep IFMA operands register-register.
- ASM avoids stack (`0` `%rsp` refs), but still pays higher memory-form IFMA pressure.

### Findings Summary

- C++ is faster today for both core and core+normalize.
- C++ is winning mostly because it achieves better register assignment in the hot path:
  - register-register IFMA in redc rounds (`0` memory-operand `vpmadd52`)
  - bounded stack scratch to keep IFMA operands resident
- ASM currently avoids stack spills entirely, but pays for it with too many memory-form IFMA ops (`72`), which lowers effective throughput.
- arithmetic work is equivalent (`232` IFMA each), so the gap is scheduling/operand-form, not operation count.

### What Must Change For ASM To Win

1. Reduce memory-form IFMA further (target: near zero in redc hot rounds).
2. Allow selective, bounded scratch if it enables register-register IFMA (the same tradeoff C++ is using).
3. Tighten constant handling:
   - keep `p[j]` constants in registers only for the shortest useful window
   - avoid repeated reload/rematerialization in inner reduction rounds
4. Keep carry-chain scheduling overlapped with independent IFMA groups (avoid clustering dependency-heavy ops).
5. Re-measure after each change with objdump counters + throughput and reject edits that only move pressure around.

### Implications For K=2

- K=2 can provide more ILP, but it also doubles live-state pressure.
- If K=1 still relies on many memory-form IFMA ops, naive K=2 duplication will likely amplify memory pressure and fail to scale.
- K=2 should be built from the K=1 winner pattern:
  - round-level interleaving of two streams
  - explicit live-range budgeting
  - selective scratch allowed when it removes memory-form IFMA
- Practical gate before full K=2 rollout:
  - bring K=1 ASM closer to C++ operand form (substantially below current `72` memory-form IFMA) to create headroom for two-stream scheduling.

### C++ Optimization Track (Kernel-Only)

C++ is currently the fastest K1 baseline. Further C++ work should target:
1. Keep IFMA operands register-register in redc rounds.
2. Minimize stack traffic by shortening temporary lifetimes (`y`/carry staging).
3. Preserve deterministic round structure so compiler keeps a stable schedule.
4. Hoist constant setup pressure:
   - keep broadcast constants in the smallest possible live window.
   - avoid re-materializing constant vectors in hot sub-rounds.
5. Re-check instruction mix and variance after each change; reject changes that only move work between stack and memory-form IFMA.

Validation:
- throughput benchmark with low variance (`min_time=1s`, `repetitions=10`).
- no correctness drift.

Success:
- improved fused `k=1` throughput with stable variance.

## Phase 3: K=2 Interleaved Kernel (Benchmark-Only)

Focus:
- hand-interleaved CIOS rounds for two streams.
- no wrapper-only chaining; explicit interleaved body.

Validation:
- compare `k=2` interleaved vs sequential two `k=1` calls.
- inspect asm/perf for spills and backend stalls.

Success:
- clear gain over sequential `k=2` baseline.
- spill behavior controlled.

## Phase 4: K=2 Schedule Refinement

Focus:
- refine port utilization and dependency spacing.
- tune unroll factor and constant handling.
- carry forward K1 lessons:
  - prefer register-register IFMA.
  - use bounded scratch only when it removes memory-form IFMA pressure.

Success:
- sustained gain with low variance across repeated runs.

## Benchmark Matrix (Minimal)

Keep only kernel-attribution benchmarks:
- `BM_Fp8x1MulRedcK1_Current`
- `BM_Fp8x1MulRedcCoreK1_FusedBaseline`
- `BM_Fp8x1MulRedcK1_FusedBaseline`
- `BM_Fp8x1MulRedcCoreK2_*` (once implemented)
- `BM_Fp8x1MulRedcK2_*` (once implemented)

Run with:
- `--benchmark_min_time=1s`
- `--benchmark_repetitions=10`
- `--benchmark_report_aggregates_only=true`

## Current Targets

Near-term:
- improve fused `k=1` combined path (`core + normalize`) over current `k=1`.

Medium-term:
- achieve strong `k=2` gain from interleaving with bounded register pressure.

## What To Avoid

- no integration work while kernel design is still moving.
- no dispatch/API churn.
- no benchmark variant sprawl after winners are identified.
- no assumptions about compiler auto-interleaving; verify with asm/perf.
