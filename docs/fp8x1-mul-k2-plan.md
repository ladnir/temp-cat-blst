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
