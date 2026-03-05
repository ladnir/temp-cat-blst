# Fp8x1 mul_v2 K=2 Plan (C++ First, ASM Second)

## Goal

Improve `Fp8x1` `mul_v2` throughput for `k=2` blocks.

Priority order:
1. C++ `k=2` mul optimization.
2. Hand ASM `k=2` mul, with **v1 cloned from compiled C++ object**.

## Scope

In scope:
- `mul_v2` kernel work for `k=2`.
- benchmark-first development and A/B evaluation.
- Linux AVX-512 throughput benchmarking.
- local Windows debug correctness checks (portable AVX emulation path).

Out of scope (this phase):
- integration/dispatch changes.
- `redc_core` redesign (reuse existing `k=2` path).
- fused `mul+redc k=2` redesign.

## Current Baseline

Current `k=2` mul path is C++ and effectively executes two independent `mul_v2` calls in sequence:
- [Fp_avx.cpp](/c/Users/peter/repo/blst/blst/fields/Fp_avx.cpp:2569)
- [Fp_avx.cpp](/c/Users/peter/repo/blst/blst/fields/Fp_avx.cpp:2620)

Primary benchmark entry:
- [Fp_avx_k_bench.cpp](/c/Users/peter/repo/blst/benchmark/Fp_avx_k_bench.cpp:490)

## Measurement Rules

Use strict A/B and keep only winners.

Benchmark settings (Linux AVX-512):
- `--benchmark_min_time=1s`
- `--benchmark_repetitions=10` (or `20` for close calls)
- `--benchmark_report_aggregates_only=true`

Core correctness checks:
- local Windows debug checks first.
- remote Linux correctness checks for any asm touch.

Decision rule:
- keep change only if sustained throughput gain exceeds run variance.

## Phase 1: C++ K=2 mul Optimization

Objective:
- improve `k=2` throughput in C++ before any asm authoring.

Deliverables:
1. No-inline C++ anchor for `k=2` mul candidate(s) to stabilize codegen for analysis.
2. A/B benchmarks versus current `fp8x1_mul_v2_k<2>`.
3. Differential correctness checks against current path.

Candidate directions (C++):
1. Interleave two streams inside one `k=2` body (do not just call `mul_v2` twice unchanged).
2. Shorten live ranges for cross-term temporaries.
3. Reduce memory-form IFMA and unnecessary move pressure.
4. Keep deterministic schedule shape for compiler stability.

Exit criteria for Phase 1:
- C++ `k=2` candidate beats baseline with stable variance.

## Phase 2: ASM v1 from C++ Object (ObjClone)

This phase starts only after a stable C++ winner exists.

Hard requirement for ASM v1:
- v1 must be a direct clone from compiled C++ object output of the selected C++ `k=2` function.
- same instruction stream structure, same register allocation pattern, same stack usage pattern.
- no speculative manual rewrites in v1.

Process:
1. Build Linux AVX-512 target with fixed flags.
2. Extract symbol and disassembly from object file.
3. Write `.S` body as ObjClone (`mul k=2` core symbol).
4. Add self-check and benchmark wiring.
5. Verify parity and benchmark parity vs C++ object anchor.

Source-of-truth artifact:
- disassembly of compiled C++ `k=2` symbol from object file.

ObjClone acceptance:
- correctness parity (bit-exact).
- throughput at least parity with C++ anchor (within noise).

## Phase 3: ASM Optimization (Post-ObjClone)

Only after ObjClone parity.

Rules:
1. One structural change at a time.
2. A/B each attempt against ObjClone anchor.
3. Remove losing candidates immediately.

Focus areas:
1. reduce memory-form IFMA in hot regions.
2. improve dependency spacing (carry chain overlap).
3. tighten constant live windows.
4. preserve no-spill behavior where possible.

## Benchmark Matrix (Minimal)

Keep this small and stable:
1. current `BM_Fp8x1MulV2K/2` baseline.
2. C++ `k=2` candidate anchor.
3. ASM ObjClone `k=2`.
4. ASM candidate `k=2` (single active experiment slot).

## Artifacts and Logging

For each accepted/rejected attempt, record:
1. benchmark command.
2. aggregate throughput (`items_per_second`).
3. delta vs anchor.
4. correctness status.
5. keep/reject decision.

Store logs under `docs/` plan notes (no integration churn in this phase).

## Stop Conditions

Stop and reassess if any of the following hold:
1. C++ phase cannot produce a repeatable gain.
2. ObjClone cannot reach correctness parity.
3. ASM edits repeatedly trade wins in one metric for regressions in primary throughput.

## Immediate Next Actions

1. Add explicit C++ `k=2` mul candidate anchor function and self-check.
2. Add A/B benchmark entries for baseline vs C++ candidate.
3. Run Linux AVX-512 baseline + candidate and record first result table.

## Phase 1 Status (March 4, 2026, local Windows debug)

Implemented in benchmark-only scope (no `blst/fields` integration yet):
1. Added benchmark-only `k=2` C++ wrappers:
   - `mul_v2_k2_cpp_baseline_once`
   - `mul_v2_k2_cpp_candidate_once`
2. Added one-time differential self-check:
   - `mul_v2_k2_selfcheck_once` compares candidate against baseline path.
3. Simplified A/B to two benchmarks only:
   - `BM_Fp8x1MulV2K<2>` (current best)
   - `BM_Fp8x1MulV2K2_CppCandidate` (new candidate)

Portable emulation enablement:
1. `benchmark/Fp_avx_k_bench.cpp` guard now includes `BLST_AVX512_PORTABLE` so these C++ benchmarks run on local Windows portable builds.

Validation performed:
1. Benchmark self-check ran during benchmark execution (no mismatches).
2. `ctest --test-dir out/build/x64-Debug --output-on-failure -R "Fp|Fr"` passed (133/133).

Local benchmark command:
1. `./out/build/x64-Debug/benchmark/blst_benchmarks.exe --benchmark_filter="BM_Fp8x1MulV2K<2>|BM_Fp8x1MulV2K2_CppCandidate" --benchmark_min_time=0.1s --benchmark_repetitions=10 --benchmark_report_aggregates_only=true`

Candidate iteration log:
1. Iteration A (cross-block call-order only): candidate was slower than current best (`-1.88%` in local debug run).
2. Iteration B (streamed dual-kernel body, flattened/no inner `mul_v2` calls): candidate became faster in local debug portable runs.
3. Iteration C (remote AVX-512 release, flattened body): rejected; `-16.1%` vs current best and very high spill pressure.
4. Iteration D (staged scheduler with noinline stream helpers, remote AVX-512 release): accepted as current best candidate.
5. Iteration E (`v3` per-stream fused helper): rejected on remote AVX-512 (`v3` slower than `v2`).
6. Iteration F (`v4` cross-block wavefront schedule): rejected on remote AVX-512 (`v4` slower than `v2`).
7. Iteration G (`v5` hybrid lockstep schedule): accepted as new C++ candidate leader (small but repeatable gain vs `v2`).
8. Iteration H (`v6` full schoolbook schedule): accepted as new C++ candidate leader (beats `v5` in pinned remote runs).
9. Iteration I (`v6a` dual-stream diagonal interleave on schoolbook): rejected (regresses vs `v6`).
10. Iteration J (`v6b` schoolbook pair-carry cadence): accepted as new C++ candidate leader (small but repeatable gain vs `v6`).
11. Iteration K (`v6c` phase-split schoolbook): rejected (no stable gain vs `v6b` in clean serial run).
12. Iteration L (`v6d` seam-only cross-block pipeline): rejected (near-parity vs `v6b`, below acceptance threshold).
13. Iteration M (`v6e` unroll-by-2 seam pipeline): rejected (codegen bloat and clear regression vs `v6b`).
14. Iteration N (`v6b` ASM lane, trampoline baseline): infrastructure-only step, superseded by ObjClone.
15. Iteration O (`v6b` ASM ObjClone): implemented as real `.S` body matching C++ `v6b` object bytes.

Local debug `items_per_second` (higher is better), Iteration B:
1. command:
   - `./out/build/x64-Debug/benchmark/blst_benchmarks.exe --benchmark_filter="BM_Fp8x1MulV2K<2>|BM_Fp8x1MulV2K2_CppCandidate" --benchmark_min_time=0.1s --benchmark_repetitions=10 --benchmark_report_aggregates_only=true`
2. result:
   - `BM_Fp8x1MulV2K<2>/256`: `14.3087k/s`
   - `BM_Fp8x1MulV2K2_CppCandidate/256`: `16.7221k/s` (`+16.87%`)

Interim decision:
1. Keep benchmark-only scaffolding (current-best + candidate + self-check).
2. Retire flattened candidate variants that trigger heavy spill traffic.
3. Keep staged/noinline-stream-helper candidate as active experiment.
4. Confirm stability with repeated isolated remote runs before any integration decisions.

Remote AVX-512 release A/B (March 5, 2026, isolated runs):
1. Baseline command:
   - `./out/build/linux/benchmark/blst_benchmarks --benchmark_filter='BM_Fp8x1MulV2K<2>' --benchmark_min_time=1s --benchmark_repetitions=20 --benchmark_report_aggregates_only=true`
2. Candidate command:
   - `./out/build/linux/benchmark/blst_benchmarks --benchmark_filter='BM_Fp8x1MulV2K2_CppCandidate' --benchmark_min_time=1s --benchmark_repetitions=20 --benchmark_report_aggregates_only=true`
3. Results:
   - `BM_Fp8x1MulV2K<2>/256_mean`: `42.1866M/s`
   - `BM_Fp8x1MulV2K2_CppCandidate/256_mean`: `42.9450M/s`
   - delta: `+1.80%` candidate over baseline.

Remote AVX-512 release `v2` vs `v3` (March 5, 2026, isolated runs):
1. `baseline`:
   - `BM_Fp8x1MulV2K<2>/256_mean`: `37.4671M/s`
2. `v2` (`BM_Fp8x1MulV2K2_CppCandidate`):
   - `42.7719M/s` (`+14.16%` vs baseline)
3. `v3` (`BM_Fp8x1MulV2K2_CppCandidateV3`):
   - `39.5914M/s` (`+5.67%` vs baseline, `-7.44%` vs v2)
4. Codegen note:
   - `v2` and `v3` helper bodies both compiled without `[rsp+...]` spill traffic, so regression appears schedule-related (lost cross-stream overlap), not spill-related.

Remote AVX-512 release `v2` vs `v4` (March 5, 2026, isolated pinned runs):
1. command shape:
   - `taskset -c 0 ./out/build/linux/benchmark/blst_benchmarks --benchmark_filter='^BM_Fp8x1MulV2K2_CppCandidate/256$' --benchmark_min_time=1s --benchmark_repetitions=20 --benchmark_report_aggregates_only=true`
   - `taskset -c 0 ./out/build/linux/benchmark/blst_benchmarks --benchmark_filter='^BM_Fp8x1MulV2K2_CppCandidateV4/256$' --benchmark_min_time=1s --benchmark_repetitions=20 --benchmark_report_aggregates_only=true`
2. result:
   - `v2`: `42.8583M/s`
   - `v4`: `40.5354M/s`
   - delta: `v4` is `-5.42%` vs `v2`.
3. codegen note:
   - `v4` has no `[rsp+...]` spill references in its hot symbol, but its driver prologue/addressing is materially heavier than `v2` (many callee-saved pushes and extra pointer arithmetic), and the schedule does not recover enough IFMA latency to offset that overhead.

Remote AVX-512 release `v2` vs `v5` (March 5, 2026, pinned runs):
1. `v5` schedule:
   - keep split Karatsuba helpers (`low_high`, `middle`)
   - reorder per block to stream-local lockstep:
     - `low_high(s0) -> middle(s0) -> low_high(s1) -> middle(s1)`
2. head-to-head run (single invocation, pinned, 30 reps, `1s`):
   - `v2`: `42.3314M/s`
   - `v5`: `42.8201M/s`
   - delta: `+1.15%` for `v5`.
3. paired stability sample (5 pinned paired runs, `2s`, 1 rep each):
   - `v2` items/s: `43.0991`, `43.0197`, `42.9154`, `40.5342`, `43.0013`
   - `v5` items/s: `43.3897`, `43.2499`, `43.1806`, `43.1899`, `43.2106`
   - median delta: about `+0.49%` for `v5`.
4. codegen note:
   - `v5` hot symbol has no `[rsp+...]` spill references.
   - compared to `v2`, `v5` driver prologue is lighter (no extra callee-saved pushes), which likely contributes to the small uplift.

Remote AVX-512 release `v5` vs `v6` (March 5, 2026, pinned runs):
1. `v6` schedule:
   - full schoolbook stream kernel (more IFMA than Karatsuba split),
   - two streams per block in lockstep (`stream0`, `stream1`).
2. head-to-head run (single invocation, pinned, 20 reps, `1s`):
   - `v5`: `42.1773M/s`
   - `v6`: `43.4644M/s`
   - delta: `+3.05%` for `v6`.
3. head-to-head run (single invocation, pinned, 30 reps, `1s`):
   - `v5`: `42.4666M/s`
   - `v6`: `43.5896M/s`
   - delta: `+2.64%` for `v6`.
4. paired stability sample (5 pinned paired runs, `2s`, 1 rep each):
   - `v5` items/s: `42.4892`, `42.4418`, `42.4477`, `42.0817`, `42.5477`
   - `v6` items/s: `43.1838`, `38.7860`, `43.2896`, `43.2587`, `43.6737`
   - one clear outlier on `v6` run 2; median trend still favors `v6` (`~+1.9%`).
5. codegen note:
   - `v6` hot symbol has no `[rsp+...]` spill references.
   - disassembly shows `128 x vpmadd52luq + 128 x vpmadd52huq` in the hot `v6` function body (higher IFMA count than Karatsuba split), but with a regular schedule that improved measured throughput on this host.

Remote AVX-512 release `v6` vs `v6a` (March 5, 2026, pinned runs):
1. `v6a` schedule:
   - same full-schoolbook math as `v6`,
   - dual-stream diagonal interleave inside one kernel body.
2. head-to-head run (single invocation, pinned, 20 reps, `1s`):
   - `v6`: `43.3863M/s`
   - `v6a`: `40.7092M/s`
   - delta: `v6a` is `-6.17%` vs `v6`.
3. codegen note:
   - `v6a` hot symbol has no `[rsp+...]` spill references and same IFMA count as `v6`
     (`128 x vpmadd52luq + 128 x vpmadd52huq`),
   - but `v6a` emits noticeably more move traffic (`vmovdqa64`), consistent with the regression.

Remote AVX-512 release `v6` vs `v6b` (March 5, 2026, pinned runs):
1. `v6b` schedule:
   - full-schoolbook math like `v6`,
   - uses a 2-column carry cadence ("pair-carry") instead of single-column carry handoff.
2. head-to-head run (single invocation, pinned, 30 reps, `1s`):
   - `v6`: `41.3906M/s`
   - `v6b`: `42.7214M/s`
   - delta: `+3.22%` for `v6b`.
3. validation rerun (isolated, pinned, 20 reps, `1s`):
   - `v6`: `43.4024M/s`
   - `v6b`: `43.6978M/s`
   - delta: `+0.68%` for `v6b`.
4. variance note:
   - combined invocations showed occasional large outliers on this host, so keep/reject decisions are based on isolated pinned reruns.
5. codegen note:
   - `v6b` and `v6` are both spill-free in the hot schoolbook body,
   - improvement appears to come from a slightly better carry/update cadence rather than reduced IFMA count.

Remote AVX-512 release `v6b` vs `v6c` (March 5, 2026, pinned runs):
1. `v6c` schedule:
   - full-schoolbook math, split into three phases per stream:
     `cols 0..5`, `cols 6..10`, `cols 11..14`.
2. clean comparison command (single process, serial):
   - `taskset -c 0 ./out/build/linux/benchmark/blst_benchmarks --benchmark_filter='^BM_Fp8x1MulV2K2_CppCandidateV6(B|C)/256$' --benchmark_min_time=2s --benchmark_repetitions=20 --benchmark_report_aggregates_only=true`
3. result:
   - `v6b`: `43.5058M/s`
   - `v6c`: `43.4302M/s`
   - delta: `v6c` is `-0.17%` vs `v6b`.
4. decision:
   - reject `v6c` for now (within noise/slightly slower, no clear win over current leader).

## Revised `v6d` Proposal (Post-`v6c`)

Objective:
1. keep `v6b` arithmetic shape (same IFMA count),
2. avoid `v6c`-style phase barriers,
3. test only seam-level overlap where it might hide carry latency without bloating driver overhead.

Design (`v6d`):
1. unit of work remains per-stream full schoolbook with pair-carry cadence (`v6b` inner shape unchanged).
2. add a minimal cross-block seam pipeline only:
   - compute stream0 of block `i`,
   - compute stream0 of block `i+1` (prefetch/load overlap opportunity),
   - compute stream1 of block `i`,
   - compute stream1 of block `i+1`,
   - advance by one block and repeat.
3. do not split diagonals into phases; no extra local scopes for column groups.
4. no new temp buffers; only existing `out/a/b` pointers and fixed-stride increments.

Hard constraints:
1. no benchmark process parallelism on the server.
2. A/B must run in one serial command, optionally with `sleep 1-2` between runs.
3. reject immediately if codegen shows stack-spill growth or significantly higher move traffic.

Measurement protocol for `v6d`:
1. command shape:
   - `taskset -c 0 ./out/build/linux/benchmark/blst_benchmarks --benchmark_filter='^BM_Fp8x1MulV2K2_CppCandidateV6(B|D)/256$' --benchmark_min_time=2s --benchmark_repetitions=20 --benchmark_report_aggregates_only=true`
2. validation:
   - run the same command a second time.
3. accept only if:
   - both runs favor `v6d`,
   - median/mean uplift is outside observed noise (target `>= +0.5%`).

Expected risk/reward:
1. potential upside: small (`~0.5% to 1.5%`) if seam overlap improves issue smoothness.
2. main risk: extra pointer/live-range pressure recreates `v4`-like overhead and loses performance.

Remote AVX-512 release `v6b` vs `v6d` (March 5, 2026, pinned runs):
1. `v6d` schedule:
   - keep `v6b` pair-carry stream kernel unchanged,
   - pipeline only at block seams:
     - prime `stream0(block0)`,
     - per next block: `stream0(next)` then `stream1(current)`,
     - drain final `stream1(last)`.
2. command used (single remote command, serial):
   - `taskset -c 0 ./out/build/linux/benchmark/blst_benchmarks --benchmark_filter='^BM_Fp8x1MulV2K2_CppCandidateV6(B|D)/256$' --benchmark_min_time=2s --benchmark_repetitions=20 --benchmark_report_aggregates_only=true`
   - rerun same command after `sleep 2` in the same remote shell command.
3. run 1 result:
   - `v6b`: `43.5575M/s`
   - `v6d`: `43.6848M/s`
   - delta: `+0.29%` for `v6d`.
4. run 2 result:
   - `v6b`: `43.6348M/s`
   - `v6d`: `43.6378M/s`
   - delta: `+0.01%` for `v6d`.
5. decision:
   - reject `v6d` for now; uplift is below acceptance threshold (`>= +0.5%`) and effectively within noise.

Remote AVX-512 release `v6b` vs `v6e` (March 5, 2026, pinned runs):
1. `v6e` schedule:
   - `v6d` seam pipeline structure,
   - unroll seam loop by 2 seams per iteration.
2. command used (single remote command, serial):
   - `taskset -c 0 ./out/build/linux/benchmark/blst_benchmarks --benchmark_filter='^BM_Fp8x1MulV2K2_CppCandidateV6(B|E)/256$' --benchmark_min_time=2s --benchmark_repetitions=20 --benchmark_report_aggregates_only=true`
   - rerun same command after `sleep 2` in the same remote shell command.
3. run 1 result:
   - `v6b`: `43.3694M/s`
   - `v6e`: `38.0969M/s`
   - delta: `-12.16%` for `v6e`.
4. run 2 result:
   - `v6b`: `43.3537M/s`
   - `v6e`: `42.4986M/s`
   - delta: `-1.97%` for `v6e`.
5. codegen note:
   - `v6e` hot symbol grew substantially vs `v6b` in this build:
     - IFMA count in candidate body: `v6b=256+256`, `v6e=640+640`,
     - move count: `v6b vmovdqa64=198`, `v6e vmovdqa64=506`,
     - stack refs in candidate body remained similar (`[rsp]` count `27` in both).
6. decision:
   - reject `v6e`; unrolling-by-2 increased front-end/code-size pressure and regressed throughput.

Remote AVX-512 release `v6b` C++ vs `v6b` ASM lane (March 5, 2026, pinned runs):
1. ASM lane shape (`benchmark` scope only):
   - new symbol `blst_fp8x1_mul_v2_k2_v6b_asm` in `.S`,
   - v0 implementation is a trampoline jump to `blst_fp8x1_mul_v2_k2_v6b_cpp_entry`
     (used to validate asm wiring and A/B harness).
2. command used (single remote command, serial):
   - `taskset -c 0 ./out/build/linux/benchmark/blst_benchmarks --benchmark_filter='^BM_Fp8x1MulV2K2_CppCandidateV6B(Asm)?/256$' --benchmark_min_time=2s --benchmark_repetitions=20 --benchmark_report_aggregates_only=true`
   - rerun same command after `sleep 2` in the same remote shell command.
3. run 1 result:
   - `v6b`: `38.6408M/s`
   - `v6b-asm`: `38.4063M/s`
   - delta: `-0.61%` for asm lane.
4. run 2 result:
   - `v6b`: `38.4909M/s`
   - `v6b-asm`: `38.5575M/s`
   - delta: `+0.17%` for asm lane.
5. decision:
   - keep asm lane as infrastructure only; no performance win yet.
   - next asm step (if pursued): replace trampoline with true ObjClone body from `v6b` object output.

Remote AVX-512 release `v6b` C++ vs `v6b` ASM ObjClone (March 5, 2026, pinned runs):
1. ObjClone implementation:
   - replaced trampoline with full `.S` body generated from C++ object symbol bytes.
2. ObjClone verification:
   - extracted symbol sizes in executable: `v6b cpp = 0x9a8`, `v6b asm = 0x9a8`.
   - byte hashes for both function bodies are identical:
     - `sha256 = f1ee984e3762886dcfe3f46a1783973f1f4ebece73c0682f0bbcd062d0d188ed`
3. command used (single remote command, serial):
   - `taskset -c 0 ./out/build/linux/benchmark/blst_benchmarks --benchmark_filter='^BM_Fp8x1MulV2K2_CppCandidateV6B(Asm)?/256$' --benchmark_min_time=2s --benchmark_repetitions=20 --benchmark_report_aggregates_only=true`
   - rerun same command after `sleep 2` in the same remote shell command.
4. run 1 result:
   - `v6b`: `41.1137M/s`
   - `v6b-asm`: `43.4122M/s`
   - note: run 1 had high `v6b` variance (`cv ~6%`), treat as unstable.
5. run 2 result:
   - `v6b`: `43.5760M/s`
   - `v6b-asm`: `43.4492M/s`
   - delta: `-0.29%` for asm.
6. additional single-run sample (`1s`, 30 reps):
   - `v6b`: `43.1192M/s`
   - `v6b-asm`: `43.7634M/s`
   - delta: `+1.49%` for asm.
7. decision:
   - ObjClone parity achieved at instruction-stream level; throughput remains within host noise band.
   - proceed to manual ASM edits only with strict one-change A/B.

Direct ASM rewrite + opt #1 (`unroll-by-2`) (March 5, 2026, pinned runs):
1. direct-ASM baseline step:
   - rewrote `benchmark/Fp_avx_mul_v2_k2_v6b_asm.S` from `.byte` blob into explicit mnemonic Intel-syntax instructions.
   - kept original loop shape first (`n_blocks` counted loop), no semantic changes.
   - spill check in symbol body: no `[rsp+...]` references.
2. opt #1 attempt:
   - hand-authored block-level unroll-by-2 in ASM (two full block bodies per loop iteration + odd-block tail).
   - no C++ donor/trampoline in this step.
3. command used (single remote command, serial):
   - `taskset -c 0 ./out/build/linux/benchmark/blst_benchmarks --benchmark_filter='^BM_Fp8x1MulV2K2_CppCandidateV6B(Asm)?/256$' --benchmark_min_time=0.2s --benchmark_repetitions=8 --benchmark_report_aggregates_only=true`
4. result:
   - `v6b` (C++): `44.1684M/s`
   - `v6b-asm` (opt #1): `38.2291M/s`
   - delta: `-13.45%` for opt #1 ASM.
5. analysis:
   - regression matches prior `v6e` behavior: code-size/front-end pressure dominates any loop-control savings.
6. decision:
   - reject opt #1 for this kernel.
   - reverted ASM to single-copy direct mnemonic baseline.
7. quick post-revert confirmation:
   - command: same shape, `--benchmark_min_time=0.2s --benchmark_repetitions=6`
   - `v6b` (C++): `44.3671M/s`
   - `v6b-asm` (reverted direct baseline): `44.4096M/s` (parity, slight noise-level edge).

Direct ASM micro-optimization trials (`free/repurpose regs`) (March 5, 2026, pinned runs):
1. trial A (`r8` countdown loop-control):
   - replaced `rax` increment/compare loop bookkeeping with `r8` countdown (`dec/jne`).
   - short run showed neutral-to-mixed behavior; no repeatable uplift.
2. trial B (`zmm20/zmm21` second-stream preload):
   - preloaded second stream heads early (`[rsi-0x200]`, `[rdx-0x200]`) into free regs and consumed via reg-reg move at seam.
   - result (`0.2s`, 8 reps): `v6b=44.3497M/s`, `asm=44.0966M/s` (`-0.57%`).
3. post-trial cleanup:
   - reverted preload path.
   - rerun with only `r8` tweak still trailed in a stable sample (`v6b=44.4773M/s`, `asm=44.0838M/s`).
   - reverted loop-control tweak as well to exact direct baseline.
4. final confirmation after full revert:
   - `v6b=44.2936M/s`, `asm=44.4197M/s` (baseline parity/noise-level edge restored).
5. decision:
   - reject both micro-opts; keep exact direct baseline schedule.
