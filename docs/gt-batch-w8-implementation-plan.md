# GT Batch W=8 Implementation Plan (v2)

## Introduction and Context

This document is a build plan for a new SIMD-oriented GT batch path in `blst`.

For newcomers: the current codebase already has scalar/asm and AVX paths, but most of the
existing AVX integration targets single-operation latency (for example, one pairing at a time).
Recent measurements showed two important things:

- Low-level AVX kernels can be very strong.
- End-to-end pairing speedups are highly sensitive to higher-layer orchestration, data movement,
  and machine/platform details.

As a result, this v2 plan focuses on a separate fixed-width batch engine (`W=8`) that is designed
for throughput, with strict apples-to-apples performance accounting.

What this plan changes:

- Adds a dedicated `W=8` tower path (`Fp -> Fp2 -> Fp6 -> Fp12 -> GT`).
- Keeps AoS API boundaries for callers while using packed internal blocks.
- Requires scalar-equivalent metrics so reported speedups are comparable and defensible.

What this plan does not change:

- It does not redesign scalar APIs.
- It does not switch the default scalar execution path.
- It does not introduce variable-width template dispatch.

Read this plan as an implementation sequence with acceptance gates. Each phase has explicit
correctness and performance exits before moving to the next layer.

## Scope

This plan defines a fixed-width (`W=8`) batch engine for independent GT tuples, optimized for throughput.

Primary operations:

- GT multiply batch (`Fp12` multiply)
- GT cyclotomic-square batch (`Fp12` cyclotomic square)
- GT exponentiation batch with independent per-lane scalars

Non-goals:

- No scalar API redesign
- No variable-width templates
- No default-path replacement during this phase

Hard constraints:

- Fixed width of 8 lanes end-to-end
- Build bottom-up from `Fp -> Fp2 -> Fp6 -> Fp12 -> GT`
- No hidden scalar fallback inside AVX APIs
- No exceptions, no heap allocation in hot loops

---

## Key Principles

1. Lane stability:
- Lane `i` always maps to tuple `i` in each 8-lane block.

2. Strong types:
- No type aliases for schedule/layout types in the hot path.
- Layout types are explicit structs with documented invariants.

3. Scalar-equivalent accounting:
- Performance claims must compare equivalent work.
- For `Fp mul`, AVX comparisons must include reduction (`mul + redc`) when scalar baseline includes reduction.

4. End-to-end first:
- Kernel-only wins are not sufficient.
- We gate with AoS API throughput (pack + kernel + unpack), not kernel-only cycles.

5. Reduction strategy for AVX:
- AVX kernels should prefer widened/unreduced intermediates (`UFp`-style) across multiply-add chains.
- Reduce at explicit boundaries, not after every multiply, when correctness bounds allow.
- Do not mechanically mirror scalar reduction cadence in AVX code.
- Any early-reduction choice must be justified by measured cycles/SEO, not by scalar-code similarity.

---

## Measurement Contract (Required)

All benchmark tables must include these columns:

- Operation name
- Path (`scalar`, `avx_kernel`, `avx_e2e`)
- Cycles total
- Scalar-equivalent operations (SEO)
- Cycles per SEO
- Speedup vs scalar (`scalar cycles/SEO` divided by `candidate cycles/SEO`)

### Scalar-equivalent definitions

- `FpMulReducedSEO`: one fully reduced `Fp` multiplication result.
  - Scalar reference: `Fp::mul`
  - AVX equivalent: `FpPack8x1::mul_*` followed by `UFpPack8x1::redc`

- `GtMulSEO`: one `GT` multiplication result (`out[i] = a[i] * b[i]`)
- `GtCyclotomicSqrSEO`: one `GT` cyclotomic square result
- `GtExpSEO`: one `GT` exponentiation result for one base/scalar pair

### Benchmark environment capture

Each run must capture:

- CPU model, flags, core count, sockets
- Hypervisor presence
- Governor / turbo state
- Compiler and flags
- OS/kernel

Recommended commands:

```bash
lscpu
uname -a
cat /etc/os-release
grep -En "model name|cpu MHz|microcode" /proc/cpuinfo | head -n 80
cat /sys/devices/system/cpu/intel_pstate/no_turbo 2>/dev/null
cat /sys/devices/system/cpu/amd_pstate/status 2>/dev/null
cat /sys/devices/system/cpu/cpufreq/boost 2>/dev/null
cat /sys/devices/system/cpu/cpufreq/policy0/scaling_governor 2>/dev/null
```

---

## AVX Reduction Guidance (Design Constraint)

AVX performance depends heavily on reduction placement. The default strategy is:

- Keep intermediate values in widened form (`UFp`, `UFp2`, `UFp6`, etc.) through dependent multiply/add/sub chains.
- Perform a single reduction at the latest safe point that preserves correctness and range bounds.
- Treat reduction as a scheduled resource, not a per-operation reflex.

Implication for implementation:

- A naive line-by-line scalar port (with immediate reduction after each multiply) is not an acceptable AVX baseline.
- AVX code may differ structurally from scalar code to preserve widened pipelines, while still matching scalar results at API boundaries.

Validation requirement:

- For each major kernel, record both:
  - chosen reduction boundaries
  - an alternative earlier-reduction variant (if implemented)
  and keep the faster one by cycles/SEO, with correctness parity.

---

## API Surface (AoS Boundary)

Target public APIs:

- `gt_mul_batch(GT* out, const GT* a, const GT* b, size_t n)`
- `gt_cyclotomic_sqr_batch(GT* out, const GT* a, size_t n)`
- `gt_exp_batch(GT* out, const GT* a, const byte* scalars, size_t scalar_stride, size_t nbits, size_t n)`

Rules:

- API inputs/outputs are AoS.
- Internally process by 8-lane packed blocks.
- Tail policy default is masked AVX lanes.
- Scalar tail fallback is allowed only as a separately measured path.

---

## Phase Plan

## Phase 0: Contracts and Baseline

Deliverables:

- Freeze metric contract (SEO + cycles/SEO).
- Baseline scalar benchmarks for all target APIs.
- Baseline AVX kernel-only and AVX end-to-end benchmarks.

Exit criteria:

- Reproducible baseline numbers with environment metadata.
- CI job can run correctness and perf smoke checks.

---

## Phase 1: Fp x8 Core

Current codebase status (from review):

- Already present in `blst/fields/Fp_avx.h/.cpp`: `FpPack8x1` and `UFpPack8x1` arithmetic core (`add/sub/as*`, `mul_v2`, widened add/sub/as*, `redc`, conversions).
- Already present in `test/Fp_avx_tests.cpp`: strong differential coverage for conversions, reduced ops, widened ops, alias-equivalence for `mul` paths.
- Missing for batch-engine construction: production AoS pack/unpack helpers and lane-mask selection helpers (currently ad hoc/test-local patterns).

Phase 1 execution scope:

- Keep existing arithmetic kernels unchanged unless a correctness bug is found.
- Add only the low-level control/data-movement primitives needed to make Phases 2-7 mechanical.
- Add SEO-correct Fp throughput measurements that include both kernel-only and AoS end-to-end paths.

Phase 1 APIs to add (style aligned with `Fp_avx.h` static methods):

```cpp
struct FpPack8x1 : public AvxPack<NWORDS>
{
    static void from_aos(FpPack8x1& r, const Fp* in8);
    static void to_aos(Fp* out8, const FpPack8x1& a);

    static void from_aos_partial(FpPack8x1& r, const Fp* in, size_t n_active, const Fp& fill);
    static void to_aos_partial(Fp* out, size_t n_active, const FpPack8x1& a);

    static void select(FpPack8x1& r, const FpPack8x1& a, const FpPack8x1& b, u8 lane_mask);
};

struct UFpPack8x1 : public AvxPack<2 * NWORDS>
{
    static void select(UFpPack8x1& r, const UFpPack8x1& a, const UFpPack8x1& b, u8 lane_mask);
};
```

Lane contract for `lane_mask`:

- Bit `i` controls lane `i` (`0 <= i < 8`).
- `0` selects from `a`; `1` selects from `b`.
- Must be implemented with SIMD blend/mask ops, not scalar per-lane loops in hot paths.

Data layout and boundary rules:

- Preserve existing `8x1w` limb layout exactly.
- `from_aos`: lane `i` maps to `in8[i]`.
- `from_aos_partial`: lanes `[0, n_active)` from input, lanes `[n_active, 8)` filled with caller-provided neutral (`Fp::zero()` for add chains, `Fp::one()` for mul chains).
- `to_aos_partial`: write only `[0, n_active)`.
- No heap allocation and no hidden scalar arithmetic fallback inside these helpers.

Phase 1 files expected to change:

- `blst/fields/Fp_avx.h`
- `blst/fields/Fp_avx.cpp`
- `test/Fp_avx_tests.cpp`
- `benchmark/avxbls_bench.cpp`

Phase 1 testing plan:

- Keep all existing `FpAvx` tests passing.
- Add round-trip tests for `from_aos`/`to_aos` using random and edge vectors.
- Add partial-pack tests for `n_active = 0..8` with both fill values (`0`, `1`).
- Add `select` tests for `FpPack8x1` and `UFpPack8x1` over all 256 masks (or full deterministic subset if runtime becomes problematic).
- Add aliasing tests for `select` (`r==a` and `r==b`).
- Add tail-shape differential test: partial pack -> `mul_v2` -> `redc` -> partial unpack vs scalar lane-by-lane reference.

Phase 1 perf gates:

- Report `FpMulReducedSEO` only (`mul_v2 + redc`), never unreduced mul alone.
- Publish both:
  - `avx_kernel`: pre-packed input, `mul_v2 + redc`
  - `avx_e2e`: `from_aos + mul_v2 + redc + to_aos`
- Include SEO-normalized cycles and speedup vs scalar `Fp::mul`.
- Prefer widened fused chains (`mul...add/sub...redc`) in all microbench variants.

Phase 1 success conditions:

- Bit-exact parity with scalar `Fp` across new and existing tests.
- Deterministic lane behavior for all masks and partial sizes.
- Reproducible cycles/SEO numbers across repeated runs.
- `avx_e2e` must demonstrate positive throughput scaling with larger `n` (block streaming), even if absolute speedup is platform-dependent.

What to avoid in Phase 1:

- No template/general-width abstractions (`W` remains hard-coded at 8).
- No redesign of scalar `Fp` API.
- No `4x2w` refactors in this phase unless needed for build breakage.
- No embedding GT/Fp2/Fp6 logic into Phase 1; keep it strictly Fp x8 foundation.
- No reduction-after-every-multiply rewrites that discard widened-pipeline advantages.

---

## Phase 2: Fp2 x8

Current codebase status (from review):

- `blst/fields/Fp2_avx.h/.cpp` already contains substantial AVX kernels across multiple schedules (`8x1x1`, `4x2x1`, `2x4x1`, `2x2x2`, `1x4x2`).
- For the canonical batch path (`8x1x1`), widened operations are present (`UFp2Pack8x1x1::mul`, `add`, `sub`, `asx4`, `mul_by_u_plus_1`), but reduced-layer convenience and AoS boundary helpers are incomplete.
- In `test/Fp2_avx_tests.cpp`, several useful helpers are still test-local (pack/unpack and local widened->reduced glue), which should be promoted to production APIs for batch layering.

Phase 2 execution scope:

- Standardize `Fp2Pack8x1x1` / `UFp2Pack8x1x1` as the canonical W=8 Fp2 layer used by later Fp6/Fp12/GT batch code.
- Add missing reduced wrappers and boundary/control primitives so upper layers can avoid schedule-specific ad hoc code.
- Do not refactor or remove existing non-canonical schedules in this phase; keep them intact for compatibility.

Phase 2 APIs to add (or complete):

```cpp
struct Fp2Pack8x1x1 : public AvxPackTuple<FpPack8x1, 2>
{
    static void from_aos(Fp2Pack8x1x1& r, const Fp2* in8);
    static void to_aos(Fp2* out8, const Fp2Pack8x1x1& a);

    static void from_aos_partial(Fp2Pack8x1x1& r, const Fp2* in, size_t n_active, const Fp2& fill);
    static void to_aos_partial(Fp2* out, size_t n_active, const Fp2Pack8x1x1& a);

    static void add(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a, const Fp2Pack8x1x1& b);
    static void sub(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a, const Fp2Pack8x1x1& b);
    static void mul(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a, const Fp2Pack8x1x1& b);
    static void sqr(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a);
    static void conjugate(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a);
    static void mul_by_u_plus_1(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a);
    static void select(Fp2Pack8x1x1& r, const Fp2Pack8x1x1& a, const Fp2Pack8x1x1& b, u8 lane_mask);
};

struct UFp2Pack8x1x1 : public AvxPackTuple<UFpPack8x1, 2>
{
    static void redc(Fp2Pack8x1x1& r, const UFp2Pack8x1x1& a);
    static void select(UFp2Pack8x1x1& r, const UFp2Pack8x1x1& a, const UFp2Pack8x1x1& b, u8 lane_mask);
};
```

Data layout and tail policy:

- Lane `i` always maps to tuple `i`: `in8[i] <-> lane i`.
- Component mapping: `[0]` is real part, `[1]` is imaginary part.
- Partial-pack fill for inactive lanes:
  - additive chains: `Fp2::zero()`
  - multiplicative chains: `Fp2::one()` (not `Fp::one()` broadcast to both components)
- No scalar fallback inside these AVX APIs.

Phase 2 files expected to change:

- `blst/fields/Fp2_avx.h`
- `blst/fields/Fp2_avx.cpp`
- `test/Fp2_avx_tests.cpp`

Phase 2 testing plan:

- Keep all existing `Fp2Avx` tests passing.
- Add AoS round-trip and partial-tail tests for `8x1x1`.
- Add differential tests vs scalar `Fp2` for `sub/mul/sqr/conjugate/mul_by_u_plus_1`.
- Add lane-mask `select` tests for reduced and widened forms (full mask coverage or deterministic exhaustive subset).
- Add aliasing tests for supported `select` in-place cases.
- Add tail-shape differential tests: partial pack -> op -> partial unpack vs scalar reference.

Benchmark policy for Phase 2:

- Benchmark runs are not required at this stage.
- No performance gate blocks Phase 2 completion.
- Performance measurements resume in Phase 3+ once Fp6/Fp12 integration paths are in place.

Exit criteria:

- Full scalar parity for canonical `8x1x1` API surface.
- Clean AoS boundary for Fp2 batch blocks (`from_aos*`/`to_aos*` complete).
- No layout conversion inside arithmetic kernels.

---

## Phase 3: Fp6 x8

Types/APIs:

- `struct Fp6Pack8x1`
- `add`, `sub`, `mul`, `sqr`, `mul_by_nonresidue`

Testing:

- Differential vs scalar `Fp6`
- Stress tests for intermediate growth and reduction correctness

Perf gates:

- Kernel cycles/SEO and AoS end-to-end cycles/SEO both required
- Reduction boundaries in `UFp6` paths must be documented and tested.

Exit criteria:

- Scalar parity under random/adversarial vectors
- Measurable throughput gain in AoS end-to-end mode

---

## Phase 4: Fp12 x8

Types/APIs:

- `struct Fp12Pack8x1`
- `mul`, `sqr`, `cyclotomic_sqr`, `conjugate`

Testing:

- Differential vs scalar `Fp12`
- Cross-check vs existing AVX kernels where equivalent

Perf gates:

- End-to-end AoS cycles/SEO is required for acceptance
- Pack/unpack cost must be reported explicitly
- Avoid scalar-style eager reduction inside Fp12 chains unless it wins on cycles/SEO.

Exit criteria:

- Correctness parity
- Positive end-to-end throughput gain vs scalar

---

## Phase 5: GT Batch Kernels (mul/sqr)

Types/APIs:

- `gt_mul8(...)`
- `gt_cyclotomic_sqr8(...)`

Testing:

- Differential vs scalar GT operations
- In-place/out-of-place cases where supported

Perf gates:

- Must report both:
  - `avx_kernel` cycles/SEO
  - `avx_e2e` cycles/SEO (AoS boundary included)

Exit criteria:

- Correctness parity
- End-to-end speedup maintained for practical `n`

---

## Phase 6: GT Exponentiation x8 (independent scalars)

Types/APIs:

- `gt_exp8(GT* out8, const GT* base8, const byte* scalars8, size_t scalar_stride, size_t nbits)`

Testing:

- Differential vs scalar `GT::mult` for:
  - random scalars
  - `0`, `1`, max, sparse, dense patterns
  - multiple `nbits`

Perf gates:

- Throughput over scalar loop reported as cycles/SEO
- Optional pre-recode path measured separately

Exit criteria:

- Scalar parity across scalar classes
- Throughput gain on realistic batch sizes

---

## Phase 7: n-Size Scheduler and Tail Handling

Types/APIs:

- Finalize `gt_*_batch(..., n)` public entry points

Policy:

- Chunk by contiguous 8-lane blocks
- Tail default: masked AVX lanes
- Scalar tail path, if retained, is explicit and separately benchmarked

Testing:

- `n=0..16`, non-multiples of 8, large `n`
- Determinism checks

Exit criteria:

- Correct tail behavior for all shapes
- Linear scaling trend without hidden cliff at tails

---

## Phase 8: Parallel Throughput Scaling

Types/APIs:

- API unchanged; parallelism internal only

Policy:

- Contiguous block partitioning
- Thread-local scratch only
- No global mutable scratch

Testing:

- Determinism across thread counts
- Scaling sweep across cores

Exit criteria:

- Strong multicore scaling until hardware limits
- No correctness drift under concurrency

---

## Phase 9: Hardening and Acceptance

Deliverables:

- Freeze naming and contracts
- Document aliasing and lane contracts in headers
- CI matrix for AVX512 native and AVX512 portable modes

Testing:

- Long randomized soak tests
- Regression checks vs scalar references

Exit criteria:

- Correctness parity accepted
- Perf targets met with cycles/SEO accounting
- No scalar API regression

---

## Style Requirements

- POD structs
- Static methods
- Explicit out parameters
- Clear layout contracts in comments
- `static_assert` for trivial default constructibility/copyability
- `BLST_ASSERT` for debug-time contract checks

---

## Risks and Mitigations

Risk 1:

- Overstating speedups via non-equivalent work comparisons

Mitigation:

- Enforce SEO accounting and reduced-equivalent comparisons only

Risk 2:

- Kernel wins lost by AoS boundary overhead

Mitigation:

- Require AoS end-to-end metrics from Phase 3 onward

Risk 3:

- Tail behavior hiding AVX regressions

Mitigation:

- Masked AVX tail as default; scalar tail measured separately
