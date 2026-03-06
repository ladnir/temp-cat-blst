# AVX-512 `LUFp` Plan

This note proposes a new lazy widened domain for the AVX-512 IFMA path:
`LUFp` = lazy widened `Fp`.

The purpose of `LUFp` is to make widened linear-combination regions cheaper than
the current strict `UFp` path, while keeping the existing `UFp` contract intact.

## Motivation

Current `UFpPack8x1` add/sub kernels are not "just limb-wise add/sub."
They restore a strict bounded widened representation after every operation:

- carry propagation across all 16 limbs
- conditional correction by `p * 2^384`
- final limb masking back to radix-48 form

This is correct and useful, but expensive. Profiling in GT / Fp12 code shows:

- `UFpPack8x1::sub` is a major hotspot
- `UFpPack8x1::add` is also noticeable
- these costs are paid repeatedly in Karatsuba-style combine regions

The likely reason is not that the current kernels are poorly written. The issue is
that the current `UFp` contract is tight, so every combine step pays to re-enter it.

## Goal

Introduce a second widened domain with a looser contract:

- `UFp`: strict widened domain, unchanged
- `LUFp`: lazy widened domain for local signed linear combinations

`LUFp` should make add/sub-style combine regions cheaper, then convert back to
strict `UFp` only at explicit boundaries.

## Non-Goals

- Do not weaken the existing `UFp` contract in place.
- Do not make `redc_core` silently accept values outside its proven input model.
- Do not add per-object runtime metadata for drift or bounds.
- Do not start with a worst-case generic renorm for all callers.
- Do not redesign the tower around `LUFp` before the bounds and APIs are stable.

## Representation Model

`LUFpPack8x1` should use the same physical storage layout as `UFpPack8x1`:

- 16 limbs
- radix `R = 2^48`
- one widened value per AVX lane

The semantic difference is that `LUFp` is a lazy signed widened domain:

- limbs may be negative
- limbs may temporarily exceed 48 bits
- value is not required to satisfy the strict `UFp` bounded contract after every op

This is the key design choice. A cheap `LUFp::sub` is only interesting if it is
allowed to stay negative without immediately adding `p * 2^384`.

## Type Naming and Documentation

`LUFp` by itself is not descriptive enough, so the type must be documented very
explicitly in the header.

Suggested documentation summary:

- "Lazy widened companion for `UFpPack8x1` used for local signed linear
  combinations before returning to the strict widened domain."
- "Same physical 16-limb radix-48 layout as `UFpPack8x1`, but semantic range is
  relaxed and signed."
- "Not accepted by existing `redc_core` unless a dedicated lazy-widened reduction
  contract is later proven."

## Core Contract

The initial contract should be conservative and simple:

- `UFpPack8x1` remains the only widened type directly accepted by existing
  `redc_core`.
- `LUFpPack8x1` is closed under cheap local add/sub/fused-combine operations.
- `LUFpPack8x1::renorm<B>` reduces limb growth while staying in the `LUFp` domain.
- `LUFpPack8x1::normalize(UFpPack8x1&, const LUFpPack8x1&)` returns to the strict
  widened `UFp` contract.

The most important invariant is:

- `LUFp` may drift within a documented signed bound
- `normalize(UFp&, const LUFp&)` must restore the current `UFp` contract exactly

## API Shape

The API should follow the existing house style.

Suggested first version:

- `struct LUFpPack8x1 : public AvxPack<2 * NWORDS> {};`
- `static void add(LUFpPack8x1& r, const LUFpPack8x1& a, const LUFpPack8x1& b);`
- `static void sub(LUFpPack8x1& r, const LUFpPack8x1& a, const LUFpPack8x1& b);`
- `template <LUFpRenormClass B> static void renorm(LUFpPack8x1& r, const LUFpPack8x1& a);`
- `static void normalize(UFpPack8x1& r, const LUFpPack8x1& a);`
- `static void select(LUFpPack8x1& r, const LUFpPack8x1& a, const LUFpPack8x1& b, u8 lane_mask);`

Then add mixed and fused helpers:

- `add/sub` overloads for mixed `UFp` / `LUFp` inputs
- `aaasx2`
- `aassx2`
- `asx4`
- `ssaax2`
- `ssasx2`

## Mixed Overloads

Callers should not need an explicit `promote(...)` in hot code.

Even though `UFp` is a valid `LUFp` value semantically, forcing explicit promotion
risks unnecessary copies or awkward call sites.

Preferred approach:

- provide overloads taking `UFp` inputs directly
- share one internal implementation for both `UFp` and `LUFp` sources

Examples:

- `add(LUFp&, const UFp&,  const UFp&)`
- `add(LUFp&, const UFp&,  const LUFp&)`
- `add(LUFp&, const LUFp&, const UFp&)`
- `add(LUFp&, const LUFp&, const LUFp&)`

Same pattern for `sub` and fused helpers.

An explicit `promote` helper may still exist for debug or convenience use, but it
should not be required for performance-sensitive code.

## Renorm Policy

`renorm` should use a compile-time policy parameter, not a runtime numeric bound.

Suggested shape:

- `renorm<LUFpRenormClass::Tight>(...)`
- `renorm<LUFpRenormClass::Medium>(...)`
- `renorm<LUFpRenormClass::Loose>(...)`

Rationale:

- callers often know the drift class statically
- compile-time policy keeps codegen regular
- `if constexpr` can select the minimum carry/repair logic needed
- avoids a single worst-case renorm for every call site

The policy should describe:

- maximum expected pre-renorm limb magnitude
- postcondition target after renorm
- intended use sites

Initial intent:

- `Tight`: very short combine chains
- `Medium`: common local Karatsuba regions
- `Loose`: near the conservative operation budget limit

## Bounds Strategy

Use static bounds plus debug validation. Do not store dynamic drift metadata in
every object.

Recommended process:

1. Document a conservative signed limb bound for `LUFp`.
2. Derive per-op growth bounds for:
   - `add`
   - `sub`
   - each fused helper
3. Define a conservative operation budget before `renorm`
4. Define `renorm<B>` postconditions
5. Validate with debug telemetry and adversarial tests

Suggested initial debug constants:

- `LUFP_LIMB_BOUND = 2^61`
- `LUFP_RENORM_TARGET = 2^50`
- `LUFP_OP_BUDGET = 8`

These are starting points only. They must be tuned by measurement and stress.

## Relationship to `UFp`

This design depends on keeping the roles separate:

- `UFp`
  - strict widened domain
  - used by existing `redc_core`
  - current semantics unchanged
- `LUFp`
  - lazy widened domain
  - cheap local combine domain
  - must be normalized back to `UFp` before strict interfaces

This means:

- delaying current `UFp` correction logic does break the current `UFp` contract
- `LUFp` exists specifically so that this delay is explicit and documented

## First Use Sites

The first target should be subtraction-heavy widened combine regions, not general
rewrites of the whole tower.

Best initial candidates:

- `Fp6Pack8x1::mul` in [Fp6_avx.cpp](/c/Users/peter/repo/blst/blst/fields/Fp6_avx.cpp)
- `Fp6Pack8x1::sqr` in [Fp6_avx.cpp](/c/Users/peter/repo/blst/blst/fields/Fp6_avx.cpp)
- `Fp12Pack8x1::mul` in [Fp12_avx.cpp](/c/Users/peter/repo/blst/blst/fields/Fp12_avx.cpp)
- GT GLV4 loop paths in [GT_avx.cpp](/c/Users/peter/repo/blst/blst/curve/GT_avx.cpp) only after lower-layer wins are proven

Why these first:

- they are combine-heavy
- profiling already shows `UFp add/sub` is significant there
- many formulas are subtraction-heavy Karatsuba combinations

## Expected Outcomes

If `LUFp` works as intended:

- local widened add/sub becomes substantially cheaper than current `UFp add/sub`
- combine-heavy Fp6/Fp12/GT code pays strict `UFp` normalization less often
- reduction remains explicit and provable at boundaries

Expected measurable outcomes:

- lower share for `UFpPack8x1::sub` and `UFpPack8x1::add` in GT / Fp12 profiles
- better throughput for packed `Fp12 mul`, `Fp12 cyclotomic_sqr`, and GT scalar-mult
- no correctness regressions and no invariant violations under stress

## What To Watch

- `LUFp` can become a correctness trap if bounds are vague.
- `renorm` can erase the benefit if used too aggressively or implemented as a
  universal worst-case pass.
- `normalize(LUFp -> UFp)` must be very clearly documented as a strict boundary op.
- It is easy to overfit fused helpers without a reusable abstraction.

## What Not To Do

- Do not make `redc_core` accept `LUFp` early unless there is a written proof.
- Do not add runtime per-object drift counters in release code.
- Do not convert every existing `UFp` helper into `LUFp` at once.
- Do not start with `k>1` renorm kernels.

The likely order is:

1. single-object `LUFp add/sub`
2. `renorm<B>`
3. `normalize(LUFp -> UFp)`
4. mixed overloads
5. fused helpers
6. only later, benchmark whether `renorm_k2` is worth adding

## Benchmark Plan

Before any tower integration, benchmark:

- `LUFp add`
- `LUFp sub`
- `renorm<Tight/Medium/Loose>`
- `normalize(LUFp -> UFp)`
- `UFp add` / `UFp sub` as the baseline comparator

Then benchmark tower kernels with small targeted substitutions:

- one Fp6 combine region
- one Fp12 combine region
- one GT loop path

The realistic question is not whether isolated `LUFp add/sub` is faster. The
real question is whether `LUFp` wins after including:

- any required `renorm<B>`
- the final `normalize(LUFp -> UFp)`
- surrounding mul/select/reduction costs

So the benchmark plan should always run at two levels:

### Level 1: Primitive A/B

Measure the direct cost of:

- `UFp add`
- `UFp sub`
- `LUFp add`
- `LUFp sub`
- `LUFp renorm<Tight>`
- `LUFp renorm<Medium>`
- `LUFp renorm<Loose>`
- `LUFp normalize -> UFp`

This establishes whether the lazy widened domain is plausibly competitive before
touching tower code.

### Level 2: Realistic Kernel A/B

Measure one real combine-heavy kernel where only the widened local-combine domain
changes and everything else stays fixed.

The best first target is:

- `Fp6Pack8x1::mul` in [Fp6_avx.cpp](/c/Users/peter/repo/blst/blst/fields/Fp6_avx.cpp)

Recommended first A/B:

- baseline: current `UFp2` / `UFp` combine path
- variant: use `LUFp` only inside the subtraction-heavy widened combine region
- keep the same algorithm
- keep the same reduction boundary
- keep mul and final reduction unchanged

That gives the cleanest answer to the real design question:

- does `LUFp` beat `UFp` in an actual tower use case after paying for `renorm`
  and `normalize`?

Only after that should the plan expand to:

- one Fp12 combine region
- one GT loop slice

Use the same profiling approach already used for GT:

- pinned core
- throughput benchmark
- `perf record` / `perf report`

## Summary

`LUFp` is not a replacement for `UFp`. It is a new lazy widened domain whose job
is to make local signed combine regions cheaper.

The design should preserve the current strict `UFp` boundary, use compile-time
`renorm<B>` policies, prefer mixed overloads over explicit promotion in hot code,
and start with narrow benchmarked integrations before broader rollout.
