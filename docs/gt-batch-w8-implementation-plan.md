# GT Batch W=8 Implementation Plan

## Scope

This plan defines a fixed-width (`W=8`) batch engine for independent GT tuples, optimized for throughput.

Primary operations:

- Group addition in GT (`Fp12` multiply)
- Group doubling in GT (`Fp12` cyclotomic square)
- GT scalar multiplication with independent per-lane scalars (exponentiation)

Non-goals:

- No scalar API redesign
- No variable `W` templates
- No immediate integration into default scalar execution path

Design rule:

- Lane `i` always corresponds to tuple `i` in a batch block of 8, from `Fp` through `GT`.

---

## Phase 0: Contracts and Baseline

### APIs

Define target batch API surface first (AoS boundary):

- `gt_mul_batch(GT* out, const GT* a, const GT* b, size_t n)`
- `gt_cyclotomic_sqr_batch(GT* out, const GT* a, size_t n)`
- `gt_exp_batch(GT* out, const GT* a, const byte* scalars, size_t scalar_stride, size_t nbits, size_t n)`

### Data Layout

- External: existing AoS (`GT*`).
- Internal: packed 8-lane blocks.

### Testing

- Build scalar-vs-batch differential harness for all three operations.
- Add perf harness with stable warmup and pinning policy.

### Success Conditions

- Baseline correctness harness ready.
- Baseline throughput numbers recorded.

### What to Avoid

- Changing scalar API behavior.
- Defining hot-kernel contracts late.

---

## Phase 1: Fp x8 Core

### APIs

Using `Fp_avx.h` style (POD structs, static methods, explicit out args):

- `FpPack8x1::add`
- `FpPack8x1::sub`
- `FpPack8x1::mul` + reduction path
- `select` / mask-blend helpers required by exponentiation control flow

### Data Layout

- Strict lane-stable `8x1` layout.
- No schedule-mixed operand packing in the primary batch path.

### Testing

- Random differential tests vs scalar `Fp`.
- Edge vectors (`0`, `1`, `p-1`, carry-heavy limbs).
- Aliasing tests (`r==a`, `r==b` where supported).

### Success Conditions

- Bit-exact parity with scalar `Fp`.
- No correctness regressions under high-volume random tests.

### What to Avoid

- Per-operation pack/unpack in hot loops.
- Data-dependent branches on secret values.

---

## Phase 2: Fp2 x8

### APIs

- `struct Fp2Pack8x1x1`
- `add`, `sub`, `mul`, `sqr`, `conjugate`, `mul_by_u_plus_1`, `select`

### Data Layout

- Two `FpPack8x1` components (`c0`, `c1`), lane-stable.

### Testing

- Differential tests vs scalar `Fp2`.
- Algebraic identity checks across random vectors.

### Success Conditions

- Full parity with scalar `Fp2`.
- Stable vector path behavior on edge inputs.

### What to Avoid

- Accidental fallback to schedule-oriented mixed layouts for main path.

---

## Phase 3: Fp6 x8

### APIs

- `struct Fp6Pack8x1`
- `add`, `sub`, `mul`, `sqr`, `mul_by_nonresidue`

### Data Layout

- Three `Fp2Pack8x1x1` components (`c0`, `c1`, `c2`), lane-stable.

### Testing

- Differential tests vs scalar `Fp6`.
- Stress tests for temporary growth and reduction correctness.

### Success Conditions

- Scalar parity maintained under randomized and adversarial vectors.

### What to Avoid

- Any AoS conversion inside arithmetic kernels.

---

## Phase 4: Fp12 x8 and GT Add/Doubling

### APIs

- `struct Fp12Pack8x1`
- `mul`, `sqr`, `cyclotomic_sqr`, `conjugate`
- `gt_mul8(...)`
- `gt_cyclotomic_sqr8(...)`

### Data Layout

- Two `Fp6Pack8x1` components (`c0`, `c1`), lane-stable.

### Testing

- Differential vs scalar `Fp12` and GT operations.
- Cross-check against existing `Fp12_avx` results where comparable.

### Success Conditions

- Correctness parity and measurable single-thread throughput gain.

### What to Avoid

- Using wrapper-style per-op repacking as the core implementation strategy.

---

## Phase 5: GT Exponentiation x8 (Independent Scalars)

### APIs

- `gt_exp8(GT* out8, const GT* base8, const byte* scalars8, size_t scalar_stride, size_t nbits)`

Optional helper:

- Pre-recoded per-window digit buffers for `8` lanes.

### Data Layout

- Packed GT lane blocks.
- Independent scalars per lane.

### Testing

- Differential vs scalar `GT::mult` across:
  - random scalars
  - `0`, `1`, max, sparse, dense patterns
  - multiple `nbits` values

### Success Conditions

- Scalar parity for all scalar classes.
- Throughput gain over scalar loop on realistic batch sizes.

### What to Avoid

- Assuming same scalar across lanes.
- Secret-dependent branching or table indexing.

---

## Phase 6: `n`-size Batch Scheduler

### APIs

- Finalize public batch entry points from Phase 0.

### Data Layout

- Chunk into 8-lane blocks.
- Tail handling for `n % 8` via masked lanes or scalar fallback.
- Pack once per block, unpack once per block.

### Testing

- All `n` shapes:
  - `n=0..16`
  - non-multiples of 8
  - large `n`
- In-place and out-of-place variants where supported.

### Success Conditions

- Correct outputs for all tail scenarios.
- Linear scaling trend with `n`.

### What to Avoid

- Heap allocation in inner loops.
- Hidden global mutable scratch state.

---

## Phase 7: Parallel Throughput Scaling

### APIs

- Keep batch API stable; parallelism remains internal.

### Data Layout

- Partition work by contiguous block ranges.
- Thread-local scratch only.

### Testing

- Determinism across thread counts.
- Throughput scaling sweep across cores.

### Success Conditions

- Strong multicore scaling until hardware limits.
- No correctness drift under concurrency.

### What to Avoid

- False sharing.
- Tiny task granularity that dominates with scheduler overhead.

---

## Phase 8: Hardening and Acceptance

### APIs

- Freeze kernel contracts and naming.
- Document lane and aliasing contracts in headers.

### Data Layout

- Lock lane-stable invariants in comments and tests.

### Testing

- Long randomized soak runs.
- CI matrix for AVX512 native and portable modes.
- Regression gates against scalar references.

### Success Conditions

- Correctness parity accepted.
- Throughput targets met.
- No scalar API regression.

### What to Avoid

- Silent fallback paths that hide performance regressions.

---

## Style Requirements (Match `Fp_avx.h`)

- POD structs
- Static methods
- Explicit out parameters
- Clear per-function layout contracts in comments
- `static_assert` for trivial copyability/default constructibility
- `BLST_ASSERT` for debug-time contract checks

---

## Key Risk and Mitigation

Risk:

- Reusing existing schedule-packed tower kernels as the primary batch representation can break lane-stable throughput design.

Mitigation:

- Keep a dedicated lane-stable `x8` tower path for batch kernels.
- Use existing AVX kernels as arithmetic references and building blocks, not as representation constraints.
