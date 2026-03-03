# AVX-512 IFMA Throughput Plan (Revised)

This plan is ordered by expected throughput impact, with explicit arithmetic invariants first.
The key objective is to recover relative SIMD speedup at tower level (Fp12/GT), not only at Fp.

## Scope

- 8x1 lane-wise SoA hot path.
- Large-batch throughput objective (FFT-in-the-exponent / batched GT work).
- No shuffle-heavy layout transforms in hot arithmetic unless data proves a net win.

## Canonical Numeric Model (Must Be Fixed Before Implementation)

Use one explicit model and keep it consistent across code and tests.

### Stored representation

- `Fp`: radix-48 limbs, canonical reduced form.
- `LFp`: radix-48 limbs, lazy-reduced form (not canonical), bounded.
- `AccFp`: widened signed accumulator limbs for sum-of-products.

### Compute representation

- IFMA (`madd52`) may use aligned/shifted temporaries internally.
- Temporary compute bounds are separate from stored bounds and must be documented.

### Invariants table (parameters and constants)

- Core parameters: `t` (limb count), `R = 2^48`, `c_l` (lazy range multiplier), `B_acc` (accumulator limb bound), `ACC_MULADD_BUDGET`.

- `Fp`: each limb `< 2^48`; value canonical mod `p`.
- `LFp`: each limb `< 2^48`; value in `[0, c_l * p)`.
- `AccFp`: signed limbs with `|limb| < B_acc`; non-canonical, non-modular until reduced.
- `Wide`: 2t-limb widened product/accumulator input to reduction core.

## Bounds and Safety Workflow (Required)

Do not set lazy budgets heuristically without this process.

1. Fix the limb/radix model and carry schedule in writing.
2. Derive a conservative per-op growth bound `Delta` for:
   - `acc_mul_add`
   - `acc_mul_sub`
   - `acc_add/sub`
3. Derive safe budget from `B_acc / Delta` with margin.
4. Define `acc_renorm` exactly (carry propagation rules, post-condition bounds).
5. Validate with debug instrumentation and adversarial fuzz:
   - track max absolute accumulator limbs
   - assert invariants continuously
   - run high-value edge patterns and long random compositions

Initial defaults are allowed (`ACC_MULADD_BUDGET=32`), but must be justified and then tuned.

### Practical policy

- Start conservative immediately (safety first).
- Instrument from day one (debug max-limb telemetry and invariant asserts).
- Refine early, not late: tune bounds/budget after Phase 1 microbench + stress results, before broad Fp12/GT rollouts.
- Keep release behavior stable while tuning by changing only documented constants and re-running the same A/B suite.

## Concrete Invariants (Initial Conservative Constants)

- Initial implementation constants (conservative, then tuned by measurement):
- `R = 2^48`
- `c_l = 4` for `LFp` (target value range `[0, 4p)`)
- `B_acc = 2^61` for signed accumulator limbs
- `ACC_MULADD_BUDGET = 32`
- `acc_renorm` postcondition target: each accumulator limb satisfies `|limb| < 2^50`, and represented value is unchanged.
- These are starting constraints, not final optima. Keep debug instrumentation enabled until final constants are locked.
- Bound tuning checkpoint: complete at end of Phase 1, and re-check after first Phase 2 implementation.

## Unified Dataflow and API

Use one conceptual pipeline across phases.

- Internal fast path values: `LFp`
- Sum-of-products: `AccFp` / widened `Wide`
- Reduction: `reduce_wide_core`
- Canonicalization: `normalize` only at boundaries

### Core API set

- `mont_redc_core(LFp& out, const UFp& in)`
- `mont_redc_core(LFp& out, const AccFp& in)`
- `mont_redc_core_x2(LFp& out0, LFp& out1, const UFp& in0, const UFp& in1)`
- `mont_redc_core_x2(LFp& out0, LFp& out1, const AccFp& in0, const AccFp& in1)`
- `normalize(Fp& out, const LFp& in)`
- `acc_zero(AccFp&)`
- `acc_mul_add(AccFp&, const LFp&, const LFp&)`
- `acc_mul_sub(AccFp&, const LFp&, const LFp&)`
- `acc_add(AccFp&, const AccFp&)`
- `acc_sub(AccFp&, const AccFp&)`
- `acc_renorm(AccFp&)`
- `acc_renorm` required postcondition: bounded limbs (`|limb| < renorm_bound`) with identical represented integer value.

### Leaf multiply path

Keep one mainline leaf path:

- `mont_mul_core(LFp& out, const LFp& a, const LFp& b)` (fused Montgomery multiply)

Legacy `mul_v2 -> redc` can remain as reference/debug transition path, but not as primary optimization target.

## Phased Rollout (Risk-Gated)

Rollout order is intentionally narrow-to-wide:

1. Fp only
2. Fp2
3. Fp12
4. GT loop and exponentiation

Each phase must pass correctness + invariant checks before widening scope.

### Phase 0: Instrumentation Baseline

Deliverables:

- Stable microbench suite (cycles/op) for:
  - `Fp mul`, `redc`, `Fp2 mul/sqr`, `Fp12 mul/sqr`, `Fp12 cyclotomic_sqr`, GT LoopOnly
- Low-noise run mode (core pinning, fixed frequency where available)

Success:

- Repeated-run noise under ~1-2%
- Clear bottleneck split (mul vs reduction vs memory)

Tests:

- Repeat benchmark runs and enforce CV/variance threshold.
- Validate checksum sinks and DCE resistance in microbench harness.

### Phase 1: `LFp` + split reduction on Fp only

Deliverables:

- `LFp` type and explicit invariants
- `mont_redc_core` + `normalize`
- Fp differential tests vs existing canonical path

Success:

- No invariant violations under stress
- Fp microbench with measurable decomposition of reduction vs normalize cost
- Bound telemetry captured (`max_abs_limb`, renorm frequency/headroom), with first bound/budget refinement recorded before Phase 2 expansion.

Tests:

- Differential tests vs canonical `Fp` path.
- Invariant tests with debug limb bound checks on random + adversarial vectors.

### Phase 2: `AccFp` + reduce-once policy in Fp2

Deliverables:

- `AccFp` primitives + `acc_renorm`
- Fp2 mul/sqr rewritten to accumulate then reduce once per component

Success:

- Correctness parity
- Reduction count and cycles/op drop at Fp2 level
- Post-refinement invariants remain stable under long compositions after applying tuned constants.

Tests:

- Fp2 mul/sqr differential tests.
- Long random composition tests (thousands of chained ops) with invariant assertions.

### Phase 3: Fused Montgomery leaf (`mont_mul_core`)

Deliverables:

- Replace wide-product materialization in mainline Fp mul
- Validate no spills/regression in generated code

Success:

- Fp `mul_redc` cycles/op improvement
- Better register behavior in hot path
- Must-not-spill gate: no stack spills in inner loop of `mont_mul_core` (asm/perf inspection).
- Load/store retired per op must not regress materially vs baseline (target <= +5%).

Tests:

- Differential tests against legacy `mul_v2 -> redc` reference.
- asm/perf validation for spill behavior under release flags.

### Phase 4: Dual-stream reduction and mul where profitable

Rationale:

- Use x2 streams to hide latency of carry-dependent reduction steps, but only after Phase 3 has reduced register pressure enough to avoid spills.

Deliverables:

- `mont_redc_core_x2` first (high confidence)
- optional dual-stream mul only if spill-free gains remain

Success:

- x2 microbench win and net kernel-level win in Fp2/Fp12

Tests:

- Side-by-side x1 vs x2 microbench and correctness checks.
- Spill guard re-check for x2 kernels.

### Phase 5: Fp12 integration with accumulators

Deliverables:

- Rewrite Fp12 mul/sqr and cyclotomic_sqr around `AccFp`/`AccFp2`/`AccFp4`
- Reduce once per final component boundary

Success:

- Material Fp12 gain toward target relative speedup
- GT loop benefits without semantic regressions

Tests:

- Fp12 mul/sqr/cyclotomic_sqr differential tests.
- Random GT-derived inputs for subgroup-relevant correctness validation.

### Phase 6: GT loop upgrades

Deliverables:

- Integrate optimized Fp12 kernels into GT scalar-mult loop
- keep additive GT semantics explicit in docs/APIs

Success:

- End-to-end GT throughput improvement with bounded variance

Tests:

- GT scalar-mult differential tests against scalar reference.
- Bench split into table-init / predecode / loop-only to localize wins.

### Phase 7: Optional algorithmic upgrades

- Cyclotomic-specific multiplication in GT loop
- Windowed GLV4 recoding with reusable precompiled plans

Apply only if profiling shows arithmetic bottleneck remains after Phases 1-6 and table bandwidth is not already limiting.

Tests:

- Equivalence tests for windowed recoding and precompiled plans.
- Bench across reuse patterns to confirm cache/table behavior.

## Benchmark Methodology

Primary metrics:

- cycles/op for micro-kernels
- GT LoopOnly throughput
- GT full scalar-mult throughput
- optional hardware counters for spills/load-store pressure

### IFMA Roofline Baseline (native AVX-512, recorded 2026-03-03)

These numbers are the current reference for IFMA-only behavior (`VPMADD52LOQ` / `VPMADD52HIQ`) from dedicated microbenchmarks:

- Preliminary (single-step chain, more loop-overhead-sensitive):
  - `BM_IFMA52Latency_Lo_DepChain`: `2.99 ns` per `madd52`
  - `BM_IFMA52Latency_Hi_DepChain`: `2.98 ns` per `madd52`
  - `BM_IFMA52Latency_LoHi_DepChain`: `4.06 ns` per 2 `madd52` ops
- Independent-stream throughput plateau:
  - `Lo`: ~`6.6-6.8 G madd52/s` (best observed ~`6.72 G/s`)
  - `Hi`: ~`6.5-6.9 G madd52/s` (best observed ~`6.82 G/s`)
  - `Lo+Hi`: ~`6.6-6.8 G madd52/s` overall (best observed ~`6.84 G/s`)

Use `6.8e9 madd52/s` as the current practical IFMA ceiling for planning.

### Core discovery: dependency chains collapse IFMA throughput

Unrolled dependency-chain measurements (`unroll=8`) show a much lower effective IFMA rate than independent streams:

- `BM_IFMA52Latency_Lo_DepChain`: `10.8 ns` for 8 ops -> ~`0.742 G madd52/s`
- `BM_IFMA52Latency_Hi_DepChain`: `10.8 ns` for 8 ops -> ~`0.740 G madd52/s`
- `BM_IFMA52Latency_LoHi_DepChain`: `19.7 ns` for 16 ops -> ~`0.813 G madd52/s`

Compared to the independent-stream plateau (~`6.8 G madd52/s`), dependency chains are about `8x-9x` slower.

This is a primary design driver:

- Minimize long carry-dependent chains in hot loops.
- Reduce reduction boundaries and defer canonicalization.
- Favor accumulator/lazy strategies that increase independent work between renorm/reduce steps.
- Treat low-level IFMA scheduling as secondary unless a kernel is already close to the IFMA ceiling.

### K-stream finding (LoHi chain PoC): scale streams, not source-order interleaving

`BM_IFMA52LoHi_(SeqChains|HandInterleavedChains)<k>` with `unroll=8`:

- `k=2`:
  - Seq: ~`1.781 G madd52/s`
  - Hand-interleaved: ~`1.786 G madd52/s`
- `k=4`:
  - Seq: ~`3.571 G madd52/s`
  - Hand-interleaved: ~`3.560 G madd52/s`
- `k=8`:
  - Seq: ~`6.005 G madd52/s`
  - Hand-interleaved: ~`6.016 G madd52/s`

Conclusions:

- Increasing independent streams (`k`) is the dominant lever; `k=8` approaches IFMA ceiling behavior.
- Manual source-level "hand interleaving" vs sequential chain blocks is a wash at this stage (differences within noise, ~sub-1%).
- For now, prioritize creating enough independent arithmetic streams in real kernels over hand-stamped instruction ordering.

### Real-kernel benchmark priority (next)

Use this order to maximize actionable signal:

1. `mul+redc` (first target): closest to real hot path cost and captures interaction between multiply and reduction chains.
2. `redc_core`: isolate reduction/carry latency to quantify how much of `mul+redc` is reduction-bound.
3. `redc` (`redc_core + normalize`): quantify boundary canonicalization overhead separately.
4. `mul` (`mul_v2` wide-only): useful as an upper bound, but lowest priority for end-to-end optimization decisions.

Primary decision metric:

- If `mul` is near IFMA roofline but `mul+redc` is not, reduction/lazy strategy is the main optimization target.

### Comparing kernel results to IFMA ceiling

For a kernel with `N_madd52` IFMA ops per operation:

- IFMA-only lower bound:
  - `t_ifma_floor = N_madd52 / 6.8e9` seconds/op
- Measured overhead factor:
  - `overhead = t_measured / t_ifma_floor`

Interpretation:

- `overhead ~= 1`: near IFMA-bound.
- Large `overhead`: bottleneck is elsewhere (reduction/carry chain, extra non-IFMA ops, spills, memory).

Parameter sweeps:

- `ACC_MULADD_BUDGET`: 8, 16, 32, 64
- renorm placement policy
- K interleave (only after kernels are lean): 1, 2, 4, 8

A/B rule:

- Keep benchmark set small and stable while iterating.
- Remove low-value legacy variants once a direction is clearly dominated.

## What to Avoid

1. Canonicalization after every arithmetic operation.
2. Mixing multiple inconsistent limb/radix models.
3. Introducing shuffle-heavy transforms in hot SoA arithmetic.
4. Depending on compiler ILP without checking spills/asm.
5. Expanding rollout beyond phase gate without invariant evidence.
