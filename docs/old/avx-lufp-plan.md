# AVX-512 `LUFp` Plan

This document is the execution plan for the AVX-512 lazy widened field work.
It replaces the earlier append-only notes with a single current plan.

Companion methodology note:

- [lazy-renorm-strategy.md](/c/Users/peter/repo/blst/docs/lazy-renorm-strategy.md)

## Purpose

Introduce a lazy widened domain for the AVX-512 IFMA path:

- `UFp`: strict widened domain, unchanged
- `LUFp`: lazy widened domain for local linear-combination regions

The purpose of `LUFp` is to remove repeated strict widened repair from hot
combine regions while keeping the current `UFp` contract intact at explicit
boundaries.

## Current Status

Implemented today:

- `LUFpPack8x1` exists in the core library as a documented lazy widened
  companion to `UFpPack8x1`
- mixed `UFp` / `LUFp` `add/sub` support exists
- bounded `LUFp -> UFp` normalize support exists in the core library
- focused `LUFp` tests exist in [Fp_avx_tests.cpp](/c/Users/peter/repo/blst/test/Fp_avx_tests.cpp)
- benchmark-only `Fp6` experiments exist in
  [Fp12_tower_avx_bench.cpp](/c/Users/peter/repo/blst/benchmark/Fp12_tower_avx_bench.cpp)

What is not settled yet:

- the long-term renorm API
- the correct bound classes for schedule-preserving `Fp6` paths
- the final repair kernel shape for those classes
- whether `Fp6 mul` and `Fp6 sqr` can beat the current packed baseline with a
  real bounded implementation

## Goals

- preserve the current strict `UFp` contract
- make widened add/sub-style combine regions cheaper than strict `UFp`
- keep all repair and conversion boundaries explicit and documented
- derive exact source and target bounds for measured kernels
- preserve the baseline `Fp6` schedules as closely as possible

## Non-Goals

- do not weaken the existing `UFp` contract in place
- do not make `redc_core` silently accept broader widened inputs
- do not add runtime drift metadata to each object
- do not integrate benchmark-only mock paths into the library
- do not roll out a generic renorm family before the first real `Fp6` wins are
  proven

## What The Benchmark Work Has Established

The benchmark file
[Fp12_tower_avx_bench.cpp](/c/Users/peter/repo/blst/benchmark/Fp12_tower_avx_bench.cpp)
has established the following:

- strict widened `UFp` repair is a real cost center
- primitive add/sub microbench alone is not enough to predict tower wins
- blanket renorm after lazy operations is too expensive
- schedule preservation matters more than ad hoc `LUFp` rewrites
- the correct near-term ceiling is the baseline schedule with widened repair
  removed

These are the important benchmark categories:

- real `LUFp`
  - benchmark path that uses real lazy add/sub and real `LUFp -> UFp`
    normalize
- `BaselineNoRepair`
  - benchmark-only copy of the baseline widened schedule with strict widened
    repair deleted
  - this is the current ceiling/reference
- `LUFp_NoNormalize`
  - explicit parity/mock harness only
  - useful to verify that a lazy path collapses to the no-repair baseline when
    all repair is disabled
- `LUFp_RenormNoNormalize`
  - explicit negative control
  - showed that renorm-after-every-op is not a viable direction

These paths are not interchangeable. Only the real `LUFp` path is an
implementation candidate. The others are benchmark methodology tools.

## Representation Model

Model a widened value as:

- `X = lo + R * hi`
- `R = 2^384`

Track bounds separately:

- low-half slack in units of `R`
- high-half slack in units of `p`

Uniform naming:

- `LNegA_PosB`: `-A*R <= lo < B*R`
- `HNegA_PosB`: `-A*p <= hi < B*p`

Examples:

- narrow `Fp` embedded in widened form:
  - `LNeg0_Pos1_HNeg0_Pos0`
- strict widened `UFp`:
  - `LNeg0_Pos1_HNeg0_Pos1`

This notation should be used in:

- helper contracts
- debug assertions
- bound worksheets
- benchmark notes

## Why Low Slack Matters

If low limbs are forced canonical after every conceptual add/sub, a single low
carry or borrow can get rounded into an extra full `+p` or `-p` in the high
bound. That is too pessimistic.

Allowing low-half slack means:

- low carry/borrow can be deferred
- high bounds do not inflate prematurely
- one final bounded repair can often replace many intermediate repairs

## Mixed Wide And Narrow Inputs

Treat narrow as widened with exact zero high half:

- narrow embedded: `LNeg0_Pos1_HNeg0_Pos0`

So:

- wide + narrow first widens low slack
- wide - narrow first widens low slack
- high drift should only be charged when a repair boundary deliberately folds
  low slack back into the chosen target class

## Repair Framing

`renorm` and `normalize` are both repair kernels.

They differ by target class:

- `renorm`
  - repairs into a cheaper lazy class
- `normalize`
  - repairs into a stricter final consumer-legal class

Neither should be assumed to use repeated single-`p` correction. Depending on
the source and target bounds, both may need scaled-`p` repair.

The methodology for repair placement and repair-shape choice lives in
[lazy-renorm-strategy.md](/c/Users/peter/repo/blst/docs/lazy-renorm-strategy.md).

The current boundary design rule is:

- keep the hot path as raw limb-wise add/sub with no early physical bias
- statically track the lazy source class as before
- at the boundary, choose a precomputed bias so that after carry the value lands
  in a sign-friendly one-sided class such as `[-k p, p)`
- then repair using sign-driven stages:
  - linear `+p` if `k` is very small
  - binary `+2^i p` if `k` is larger
- fall back to upper-threshold lex compare only if a source class cannot be
  mapped efficiently into that sign-driven form

## Repair Implementation Options

The repair implementation is not fixed yet. For a given source and target class,
the candidates are:

- sign-driven one-sided ladder after a precomputed bias
- linear `+p` sign ladder
- binary `+2^i p` sign ladder
- lex-compare-based ladder only as a reference/fallback
- ternary symmetric ladder
- independent two-sided ladder

We are not choosing one globally yet. The source and target bound contract
should stay fixed; the chosen repair kernel for that contract should be driven
by benchmark results.

## Current Design Direction

The current best hypothesis is:

- `Fp6 mul`
  - preserve the baseline widened schedule exactly
  - try to avoid intermediate `renorm`
  - use one final repair if the derived bound allows it
- `Fp6 sqr`
  - preserve the baseline widened schedule exactly
  - may require one selective renorm, most likely in the long `c2` path

This is intentionally different from "generic renorm first."

## Immediate Execution Order

The next work should happen in this order:

1. exact bound worksheet for schedule-preserving `Fp6 mul`
2. final repair kernel for the actual `Fp6 mul` output class
3. schedule-preserving `Fp6 mul` implementation using that repair kernel
4. exact bound worksheet for schedule-preserving `Fp6 sqr`
5. if required, one selective renorm plus final repair for `Fp6 sqr`

The plan does not currently start with:

- generic renorm-family rollout
- broad tower integration
- lazy narrow implementation

## Phase Plan

### Phase 0: Benchmark Hygiene

Implement:

- keep `BaselineNoRepair` as the explicit performance ceiling
- keep `LUFp_NoNormalize` only as an explicitly labeled mock/parity harness
- keep `LUFp_RenormNoNormalize` only as an explicitly labeled negative control
- do not treat benchmark-only mock paths as integration candidates

Tests gating this phase:

- preserve explicit comments in the benchmark file identifying mock/reference
  paths
- keep correctness checking on real implementation candidates only

Performance requirement:

- none

### Phase 1: Bound Worksheet For Baseline `Fp6 Mul`

Implement:

- derive the bound of every widened intermediate in the current baseline `Fp6`
  `mul` schedule
- use `LNegA_PosB_HNegC_PosD` notation throughout
- include mixed wide/narrow transitions explicitly
- include `mul_by_u_plus_1` explicitly
- make the carry transform explicit:
  - the common raw source class is `LNeg4_Pos4_HNeg4_Pos4`
  - after low carry folding, the injected quotient is `q = floor(L / R)` with
    `q in {-4, -3, -2, -1, 0, 1, 2, 3}`
  - so the carried high class is `H' = H + q`, i.e. `[-4p - 4, 4p + 3)`
- determine the narrowest final repair target that lets `mul` avoid
  intermediate renorm
- maintain the working derivation in
  [fp6-mul-bound-worksheet.md](/c/Users/peter/repo/blst/docs/fp6-mul-bound-worksheet.md)

Tests gating this phase:

- every bound step must correspond to a concrete operation in the baseline
  widened schedule
- add debug assertions or helper comments for every precondition that will later
  be enforced in code
- distinguish clearly between:
  - the pre-carry grouped source class
  - the post-carry high class actually consumed by repair
- do not round the small carry-induced slack `[-4, 3]` into a whole extra
  multiple of `p` without an explicit justification

Performance requirement:

- none

### Phase 2: Final Repair Kernel For `Fp6 Mul`

Implement:

- one real repair kernel for the exact carried `Fp6 mul` output class
- treat the current `1/2/4-step` repair ladders as scenario probes only; they
  are not yet the proven `Fp6 mul` boundary
- benchmark candidate repair shapes for that class only:
  - sign-driven one-sided ladder after a precomputed bias
  - linear `+p` sign ladder
  - binary `+2^i p` sign ladder
  - lex-compare-based ladder only as a reference/fallback
- keep the source and target bound contract fixed while comparing kernels
- the base case must collapse to the same structure as the current `UFp add`
  kernel when the source class is the exact `UFp` add class

Tests gating this phase:

- exact equality with a scalar/reference repair for all values in the documented
  source class
- exact equality with a scalar/reference repair for the carried class implied by
  `LNeg4_Pos4_HNeg4_Pos4 -> [-4p - 4, 4p + 3)`
- edge coverage for the negative and positive extremes of the source class
- exact postcondition proof for the target class, not only modular equivalence

Performance requirement:

- must be materially faster than the current bounded normalize path for the same
  class
- for the exact one-add source class, the sign-driven boundary should be at or
  near parity with the existing `UFp add` structure
- if it is not materially faster, stop and revisit the bound class before
  integrating it into `Fp6 mul`

### Phase 3: Schedule-Preserving `Fp6 Mul`

Implement:

- rebuild `Fp6 mul` as a line-by-line preservation of the existing widened
  baseline schedule
- replace only the targeted widened combine repair with the new bound-driven
  lazy path
- avoid intermediate renorm unless the worksheet proves it is required

Tests gating this phase:

- exact equality against current `Fp6Pack8x1::mul`
- randomized and edge-biased corpus
- debug assertions that every intermediate entering the final repair kernel
  satisfies the documented source bound

Performance requirement:

- first checkpoint:
  - within `5%` of `BaselineNoRepair`
- continuation bar:
  - beat the current packed baseline by at least `10%`

### Phase 4: Bound Worksheet For Baseline `Fp6 Sqr`

Implement:

- derive the bound of every widened intermediate in the current baseline `Fp6`
  `sqr` schedule
- explicitly analyze the long `c2` chain
- explicitly analyze split and fanout points
- determine whether `sqr` can also use final-repair-only or whether one
  selective renorm is forced

Tests gating this phase:

- every bound step must correspond to a concrete operation in the baseline
  widened schedule
- every claimed forced renorm point must be justified by a documented
  precondition failure, not by guesswork

Performance requirement:

- none

### Phase 5: Minimum-Repair `Fp6 Sqr`

Implement:

- preserve the baseline widened `Fp6 sqr` schedule exactly
- if the worksheet allows it, use final repair only
- if not, introduce exactly one real selective renorm at the justified point and
  no earlier

Tests gating this phase:

- exact equality against current `Fp6Pack8x1::sqr`
- targeted edge cases around the deepest `c2` path
- debug assertions for the selective renorm precondition and the final repair
  precondition

Performance requirement:

- first checkpoint:
  - within `5%` of `BaselineNoRepair`
- continuation bar:
  - beat the current packed baseline by at least `8%`

### Phase 6: Cleanup And Generalization

Implement only after `Fp6 mul` and `Fp6 sqr` show real wins:

- replace temporary policy-style names in measured paths with exact bound-based
  kernel naming
- remove or quarantine obsolete benchmark-only paths
- consider broader tower integration
- consider lazy narrow based on the same topology-aware methodology

Tests gating this phase:

- every measured helper must have an exact bound contract
- benchmark-only mock paths must stay clearly labeled or be removed

Performance requirement:

- no regression relative to the winning `Fp6` implementations

## Testing Rules

For every real implementation phase:

- test arithmetic equality against the current implementation
- test exact postcondition, not just modular equality, for repair kernels
- include randomized inputs and edge-biased inputs
- use debug assertions for all documented preconditions

For benchmark-only mock paths:

- they must be explicitly labeled as mock/reference/negative-control paths
- they must not be confused with integration candidates

## Performance Rules

Primitive microbench is a sanity check, not the final decision maker.

The real decision points are:

- does the repair kernel beat the current bounded normalize for the same
  documented class
- does schedule-preserving `Fp6 mul` approach `BaselineNoRepair`
- does schedule-preserving `Fp6 sqr` approach `BaselineNoRepair`
- do the real paths beat the current packed baseline by the stated bars

## Rule For Future Experiments

Either:

- the path is fully implemented with its documented bound contracts

or:

- it must be called out explicitly as a mock experiment before it is built or
  benchmarked

This rule applies especially to:

- fake no-normalize paths
- parity harnesses
- negative-control renorm experiments

## Deferred Work

Explicitly deferred until after real `Fp6` wins:

- generic renorm-family rollout
- broad `Fp12` / GT integration
- multiply contracts for broader lazy widened inputs
- lazy narrow implementation

These remain interesting, but they are not the next execution steps.

One specific deferred opportunity is worth calling out explicitly:

- if a future lazy narrow type becomes multiply-legal under a documented bound,
  then `Fp6 mul` has a second optimization frontier before the widened combine
  region
  - specifically, the six narrow `Fp2` pre-adds that currently must be repaired
    before multiply
  - this does not change the current wide-only Phase 1/2/3 work, but it is a
    real follow-on opportunity once multiply input contracts are revisited
