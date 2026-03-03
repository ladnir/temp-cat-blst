# AVX `redc_core` K=2 Spill Analysis and Optimization Plan

## Scope
This note is only about improving `BM_Fp8x1RedcCore_HandInterleavedChains<2>` by reducing register pressure and stack spill traffic in the hand-interleaved `k=2` path.

## Baseline Data (Linux AVX512 Native)

Benchmark command:

```bash
./out/build/linux/benchmark/blst_benchmarks \
  --benchmark_filter='BM_Fp8x1RedcCore_(Seq|HandInterleaved)Chains<2>$' \
  --benchmark_min_time=1s \
  --benchmark_repetitions=10 \
  --benchmark_report_aggregates_only=true
```

Current means (latest runs):
- `SeqChains<2>`: `~18.1M ops/s`
- `HandInterleavedChains<2>`: `~22.9M ops/s`
- `AsmK2`: `~25.1M ops/s`
- `AsmK2V2`: `~26.0M ops/s`

Observed uplift:
- C++ hand `k=2`: about `+26%` over `SeqChains<2>`
- asm `k=2`: about `+39%` over `SeqChains<2>`
- asm `k=2 v2`: about `+44%` over `SeqChains<2>`

Recent change notes:
- Linux benchmark-only asm `k=2` path (`blst_redc_core_k2_asm`) is now the fastest measured variant.
- The asm path is spill-free in the kernel body.
- The latest asm revision reduced constant-load pressure by reusing `p[i]` in a temp zmm for each lo/hi pair.
- A follow-up revision removed extra `vmovdqa64` for `montw/bmask` in each round (restored memory-form ops), and hoisted final carry-mask to one loaded zmm reused by both carry chains.
- Added `AsmK2V2` A/B kernel with:
  - dual-stream round macro,
  - shared-constant scheduling across both streams per round (`p[i]` loaded once, applied to stream0+stream1),
  - 6-slot lazy queues (p0 high consumed immediately into next limb).

## Objdump Findings

### 1) Standalone scalar kernel (`UFpPack8x1::redc_core`)
Symbol: `blst::UFpPack8x1::redc_core(blst::LFpPack8x1&, blst::UFpPack8x1 const&)`

- instruction count: `361`
- IFMA count: `136` (`72 lo + 64 hi`)
- stack memory refs: `0`

Conclusion: the base single-stream kernel itself is already spill-free.

### 2) Hand `k=2` benchmark hot loop (latest)
Symbol: `BM_Fp8x1RedcCore_HandInterleavedChains<2>`
Hot range: `0x427750..0x429140`

- instruction count: `1359` (benchmark body includes harness loop/control)
- IFMA count: `272` (exactly `2 x 136`)
- stack memory refs: `292`
- stack `vmovdqa*` refs: `208`

Conclusion: the C++ hand benchmark body remains stack-heavy because it includes harness state and large inlined setup, but this is not the limiting factor for the dedicated asm kernel path.

### 3) Linux asm `k=2` kernel (`blst_redc_core_k2_asm`)

- IFMA count: `272` (`2 x 136`)
- stack memory refs in kernel body: `0`
- best measured throughput: `~23.3M redc_core/s`

Conclusion: we reached a spill-free `k=2` implementation and improved throughput over C++ hand `k=2`, but we are still backend-limited.

## Why Current No-Spill `k=2` Is Still Limited
The current asm schedule is spill-free, so the dominant limit is no longer stack traffic. The remaining issue is backend pressure from:
- IFMA dependency/carry chains, and
- constant-feed/scheduling overhead around IFMA and carry operations.

Key measured evidence (single-core `perf stat`, Linux AVX machine):
- asm before constant-reuse rewrite:
  - IPC about `1.97`
  - `topdown.memory_bound_slots` about `5.36e9`
  - `uops_dispatched.port_2_3_10` about `13.8e9`
- asm after constant-reuse rewrite:
  - IPC about `2.63`
  - `topdown.memory_bound_slots` about `1.57e9`
  - `uops_dispatched.port_2_3_10` about `9.03e9`
  - throughput: about `23.26M/s`
- asm after round/mask cleanup:
  - IPC about `2.49`
  - `topdown.memory_bound_slots` about `2.23e9`
  - throughput: about `25.06M/s` (`perf` run) / `25.15M/s` (5x aggregate benchmark run)
- asm `k=2 v2`:
  - IPC about `2.36`
  - `topdown.memory_bound_slots` about `1.44e9`
  - `topdown.backend_bound_slots` about `46.03e9`
  - throughput: about `26.11M/s` (`perf` run), `~25.96M/s` (10x aggregate run)

Interpretation:
- The prior asm variant was overusing load ports due to repeated memory-form IFMAs on constants.
- Reusing constants in-register per lo/hi pair removed much of that load-port pressure.
- After this fix, the kernel is mainly backend/compute constrained (not memory-latency constrained).

## Throughput Context
From IFMA microbenchmarks on the same host, `VPMADD52` can reach much higher rates under synthetic independent chains.

For current asm `k=2` redc-core:
- effective IFMA rate is about `23.3M * 136 = 3.17 G IFMA/s`

So there is still headroom versus pure IFMA microbench ceilings, but `redc_core` includes significant non-IFMA carry/shift/add work and strict dependencies.

Important interpretation:
- `redc_core` is not an IFMA-only kernel; non-IFMA ops (add/shift/mask/carry) and memory traffic remain large.
- So IFMA roofline alone is not an achievable bound for full `redc_core`.

## Concrete Optimization Plan (Before Implementation)

### Phase 1: Keep no-spill asm as baseline
Current baseline is `blst_redc_core_k2_asm` at about `23.3M/s`.

Design constraints:
- keep two streams interleaved (`k=2`)
- keep arithmetic identical
- keep IFMA count identical (`272` in helper)
- avoid giant all-at-once variable sets

Primary tactic:
- replace long-lifetime high accumulators (`y0..y14` / `z0..z14`) with a sliding-window representation
- consume carry/high contributions as soon as legal each round
- do not retain values past their last use

### Phase 2: Reduce constant-feed overhead further (still no spill)
Target the round schedule so one loaded `p[i]` constant can feed both streams in the same round.

Per-round intent:
1. derive `u0` and `u1` for both streams
2. load each `p[i]` once and apply to both streams before reloading next constant
3. perform carry/shift/add steps with bounded live set
4. preserve spill-free kernel body

### Phase 3: Verify with assembly + perf gates
After rewrite, collect Linux `objdump` and `perf stat` versus current asm baseline.

Required gates:
- preserve IFMA total (`272`)
- keep kernel spill-free
- improve throughput beyond current `~23.3M/s`
- no correctness regressions against `UFpPack8x1::redc_core`

### Phase 4: Re-benchmark only the critical set
Run only:
- `BM_Fp8x1RedcCore_SeqChains<2>`
- `BM_Fp8x1RedcCore_HandInterleavedChains<2>`
- `BM_Fp8x1RedcCore_AsmK2`

Acceptance target:
- first checkpoint: `>= 25M ops/s` asm `k=2`
- stretch: `>= 27M ops/s`
- `>= 30M ops/s` is still aspirational without deeper arithmetic/schedule changes

Checkpoint status:
- `>= 25M ops/s`: achieved.
- `AsmK2V2` currently improves `AsmK2` by roughly `+3.4%`.
- Latest `AsmK2V2` revision improves `AsmK2` by roughly `+2.5% .. +3.3%` depending on run mode.

## Current Status
- We now have a coherent explanation with measured counters:
  - `k=2` interleaving helps.
  - spill-free asm helps more.
  - constant-feed/load-port pressure was a major limiter and was reduced.
  - remaining limiter is backend arithmetic/dependency scheduling.
- We are not done optimizing `redc_core`, but the next target is clear: reduce round scheduling overhead while preserving spill-free behavior.

## ASM Update

We implemented and tuned a Linux benchmark-only asm kernel:
- file: `benchmark/redc_core_k2_asm.S`
- symbol: `blst_redc_core_k2_asm`
- exposed as separate benchmark `BM_Fp8x1RedcCore_AsmK2` (A/B alongside C++ baselines).

Current schedule choice:
- lazy-`y` queue model (close to scalar arithmetic),
- spill-free hot body,
- in-round constant reuse for each lo/hi IFMA pair.

Measured outcome:
- `SeqChains<2>`: ~`18.1M/s`
- `HandInterleavedChains<2>` (C++): ~`22.9M/s`
- `AsmK2`: ~`25.1M/s` (best current path)

## What Not To Do
- do not broaden to `k=4/k=8` until `k=2` is clean
- do not mix in mul+redc fusion yet
- do not infer performance from Windows asm for this decision; use Linux asm only

## Repro Commands Used for Analysis

```bash
# Benchmark
./out/build/linux/benchmark/blst_benchmarks \
  --benchmark_filter='BM_Fp8x1RedcCore_(SeqChains<2>|HandInterleavedChains<2>|AsmK2)$' \
  --benchmark_min_time=1s \
  --benchmark_repetitions=10 \
  --benchmark_report_aggregates_only=true

# Perf counters (single core)
perf stat -e \
cycles,instructions,uops_executed.core,uops_executed.stall_cycles,\
idq_uops_not_delivered.core,cycle_activity.stalls_total,\
cycle_activity.cycles_mem_any,resource_stalls.scoreboard,resource_stalls.sb,\
topdown.backend_bound_slots,topdown.memory_bound_slots,topdown.bad_spec_slots \
taskset -c 0 ./out/build/linux/benchmark/blst_benchmarks \
  --benchmark_filter='BM_Fp8x1RedcCore_AsmK2$' \
  --benchmark_min_time=2s \
  --benchmark_repetitions=1 \
  --benchmark_report_aggregates_only=true

# ASM disassembly
objdump -d -M intel --start-address=0x433c95 --stop-address=0x434300 \
  out/build/linux/benchmark/blst_benchmarks
```
