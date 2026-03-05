# Fp8x1 Fused mul+redc Kernel Plan (K=1 -> K=2)

## Goal

Optimize the fused `Fp8x1` `mul+redc` kernel for throughput.

Primary objective:
- maximize `k=1` arithmetic throughput with a no-spill hand ASM path that meets or beats the C++ baseline.

`k=2` interleaving work is explicitly deferred to a later phase.

This plan is kernel-focused only.

## Scope (Strict)

In scope:
- benchmark-only kernel development and evaluation.
- native AVX-512 (`BLST_AVX512`, non-portable), Linux x86_64 target.
- `k=1` kernel optimization only in this phase.

Out of scope:
- library/public API changes.
- dispatch/integration wiring.
- additional integration churn in this phase (existing wiring remains as-is).
- AVX portable optimization.
- GT/G1/G2 integration.

## Algorithm Choice

Current working baseline is fused **Karatsuba-style `mul_v2` + `redc_core`** (`lazy` output).

Why this choice:
- CIOS was implemented and benchmarked, but underperformed current Karatsuba-based fused kernels.
- current best throughput path is already Karatsuba-based in both C++ and K1 ASM explorations.
- this keeps K1/K2 comparisons anchored to the fastest known arithmetic structure.

Not the current target:
- CIOS/FIOS redesign unless new evidence shows a clear win.
- architecture-portable rework.

## Kernel Dataflow (Per Stream)

For `n=8` limbs:
- form widened product via Karatsuba-style `mul_v2` structure (`z0..z15` equivalent state).
- apply Montgomery `redc_core` to produce `LFpPack8x1`.
- normalize as a separately measured boundary step.

Boundary:
- keep output lazy (`LFpPack8x1`) in core benchmark.
- measure normalize separately and combined.

## Register-Pressure Strategy

Design constraints for `k=1` implementation:
- keep widened/intermediate state live ranges as short as possible.
- avoid duplicate temporaries (`w*`, `r*`) after reduction.
- keep constants minimal in live registers:
  - pin only always-used constants (`mask`, `montw`, `zero`).
  - load/use `p[j]` with short live range.
- avoid memory-form IFMA in hot reduction rounds whenever possible.

`k=2` strategy notes are kept for later, but implementation is deferred.

## Implementation Phases

## Phase 1: K=1 Fused Karatsuba Baseline (Benchmark-Only)

Deliverables:
- benchmark-only fused `k=1` Karatsuba+redc kernel baseline and checks.
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
- `BM_Fp8x1MulRedcCoreK1_FusedAsmObjClone` (object-clone ASM anchor)
- `BM_Fp8x1MulRedcCoreK1_FusedAsmCandidate` (experimental candidate)

Throughput snapshot (latest A/B):
- current non-fused baseline (`BM_Fp8x1MulRedcK1_Current`): `~12.71M/s`
- C++ noinline core: `~16.24M/s`
- C++ noinline core + normalize: `~13.72M/s`
- ObjClone ASM core: `~16.34M/s`
- ObjClone ASM core + normalize: `~13.63M/s`

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
2. Keep the hot path spill-free while still increasing register-register IFMA share.
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
  - keep no-spill allocation discipline from K1
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

### Completed K=1 Idea Results (A/B Log)

Bench setup for all entries:
- `--benchmark_min_time=0.2s` (or `1s` for stability checks), aggregate-only, Linux AVX-512.
- compare `ObjClone` vs `Candidate` on:
  - core only
  - core + normalize

1. Idea #1 (no-stack-spill macro candidate)
   - core: `16.3169M/s` (ObjClone) vs `15.1324M/s` (Candidate), `-7.3%`
   - core+normalize: `13.6141M/s` vs `12.4639M/s`, `-8.5%`
   - Result: rejected.
2. Idea #2 (tail/carry overlap tweak on ObjClone-derived candidate)
   - core: `16.4122M/s` vs `16.3479M/s`, `-0.39%`
   - core+normalize: `13.7330M/s` vs `13.6239M/s`, `-0.79%`
   - Result: rejected.
3. Idea #3 (move-pressure cleanup in one hot cluster)
   - core: `16.3437M/s` vs `16.2279M/s`, `-0.71%`
   - core+normalize: `13.6283M/s` vs `13.5991M/s`, `-0.21%`
   - Result: rejected.

### Introspection (Why Small Edits Lost)

- ObjClone and baseline both run at high IPC (`~2.12-2.15`) with very low branch-miss (`~0.2-0.3%`).
- This is not a branch/front-end bottleneck.
- Local instruction cleanups mostly perturb a near-critical backend schedule without reducing the real dependency depth.
- Practical implication: prioritize structural schedule changes over micro-edits.

### Next K=1 Design Queue (Current)

1. Two-round software-pipelined REDC (partial attempt tested)
   - Keep rounds `i` and `i+1` in flight, moving `u`/carry availability earlier relative to IFMA groups.
   - Expected: `+1% to +4%`.
   - Attempt A (partial round0/1 pipeline) result:
     - core: `16.3290M/s` (ObjClone) vs `15.0450M/s` (Candidate), `-7.9%`
     - core+normalize: `13.6284M/s` vs `12.3698M/s`, `-9.2%`
   - Status: rejected; keep idea open only for a full-round reschedule, not partial pair edits.
2. Liveness-driven register map (single-home assignment)
   - Reduce copy pressure by fixing stable register ownership across the full core.
   - Expected: `+1% to +3%`.
   - Attempt A (round0 high-term single-home rewrite) result:
     - core: `16.3306M/s` (ObjClone) vs `16.2498M/s` (Candidate), `-0.49%`
     - core+normalize: `13.6349M/s` vs `13.6141M/s`, `-0.15%`
   - Status: rejected for now.
   - Attempt B (late REDC constant pinning + register-carry for one long scratch value) result:
     - structural: stack touches reduced `32 -> 29` in candidate core asm.
     - stability run (`min_time=1s`, `repetitions=20`):
       - core: `16.4545M/s` (ObjClone) vs `16.4517M/s` (Candidate), `-0.02%`
       - core+normalize: `13.7289M/s` vs `13.6877M/s`, `-0.30%`
   - Status: rejected for throughput gate (no sustained speedup).
3. Hybrid scratch design
   - Keep hottest chain values in registers; allow a small planned scratch ring for cold values.
   - Expected: `+2% to +5%`.
4. Streaming fused mul->REDC schedule
   - Start REDC work as soon as required limbs are ready (reduce peak live state).
   - Expected: `+3% to +8%` (high effort/risk).
   - Attempt A (C++ streaming prototype, diagonal-on-demand schoolbook) result:
     - core: `16.3480M/s` (ObjClone) vs `10.0602M/s` (Streaming), `-38.5%`
     - core+normalize: `13.6984M/s` vs `9.2190M/s`, `-32.7%`
   - Status: rejected in current form (extra mul cost dominates despite earlier REDC start).
5. Internal dual-pack interleaving for k=1 throughput mode
   - Two independent packs in one body to add ILP.
   - Expected: `+5% to +12%` throughput (larger redesign).
6. Bounded-lazy tail specialization
   - Reduce terminal carry/mask work if tighter lazy bounds are proven safe.
   - Expected: `+1% to +4%`.

### Deferred Idea (Not In This Phase)

1. Two-stream interleaving (`k=2`) to hide IFMA latency with cross-stream ILP.
   - Expected: `+3% to +7%` throughput when register pressure is controlled.
   - Status: deferred until the current `k=1` hand ASM path is stabilized and faster than C++ baseline.

## Phase 3: K=2 Interleaved Kernel (Benchmark-Only, Deferred)

Focus:
- hand-interleaved fused rounds for two streams (built from K1 winner pattern).
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
  - keep no-spill register allocation discipline.

Success:
- sustained gain with low variance across repeated runs.

## Benchmark Matrix (Minimal)

Keep only kernel-attribution benchmarks:
- `BM_Fp8x1MulRedcCoreK1_FusedBaseline`
- `BM_Fp8x1MulRedcK1_FusedBaseline`
- `BM_Fp8x1MulRedcCoreK1_FusedAsmObjClone`
- `BM_Fp8x1MulRedcK1_FusedAsmObjClone`
- `BM_Fp8x1MulRedcCoreK1_FusedAsmCandidate`
- `BM_Fp8x1MulRedcK1_FusedAsmCandidate`

Deferred for later phase:
- `BM_Fp8x1MulRedcCoreK2_*`
- `BM_Fp8x1MulRedcK2_*`

Run with:
- `--benchmark_min_time=1s`
- `--benchmark_repetitions=10`
- `--benchmark_report_aggregates_only=true`

## Current Targets

Near-term:
- beat C++ on fused `k=1` core while keeping no spills.
- then improve fused `k=1` combined path (`core + normalize`).

Medium-term (deferred):
- achieve strong `k=2` gain from interleaving with bounded register pressure.

## What To Avoid

- no integration work while kernel design is still moving.
- no dispatch/API churn.
- no benchmark variant sprawl after winners are identified.
- no assumptions about compiler auto-interleaving; verify with asm/perf.
