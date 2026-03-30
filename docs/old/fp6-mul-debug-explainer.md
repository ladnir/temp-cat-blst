# `Fp6 mul` Debug Explainer

This note captures the current debugging state of the `LUFp v2` `Fp6 mul`
work in the benchmark-local v2 path.

The immediate goal is:

- recover a **correct** `LUFp v2` `Fp6 mul`
- with a final shape that is as close as possible to the existing fast path
- without losing the current performance advantage that the broken fast path
  demonstrated

There is also an explicit performance requirement:

- the final correct path should recover the current fast-path performance class
- the target is approximately `1700 ns/Fp6`
- more concretely, the frozen fast reference currently measures around
  `1680-1690 ns/Fp6`

This note is intentionally operational. It records what is known, what has been
proven experimentally, what is still unclear, and how the current candidate
paths relate to the bound worksheet.

## Current code paths

The relevant benchmark-local code lives in:

- [Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h)
- [Fp_avx_lufp_v2_fp6_bench.cpp](/c/Users/peter/repo/blst/benchmark/Fp_avx_lufp_v2_fp6_bench.cpp)

There are currently several important `Fp6` paths:

1. `Wide`
   - the **fast but broken** reference
   - this is the path that produced about `1680-1690 ns/Fp6`
   - it must stay frozen as a performance reference

2. `WideSafe`
   - the **correct but slow** reference
   - this uses the older safe `HiNeg4Pos4`-style renorm logic
   - it passes validation, but is much slower

3. `WideStepCandidate`
   - the current best **correct** step-style candidate
   - this renorms after every add/sub step
   - it is much faster than `WideSafe`, but still materially slower than `Wide`

4. `WideStepCandidateRelaxed2`
   - a correct deferred-checkpoint candidate
   - this now uses `step<2>, step<1>, step<1>` at each two-step checkpoint
   - it is correct, but slower than `WideStepCandidate`

5. `WideStepCandidateHybrid`
   - a mixed checkpoint experiment
   - also correct
   - slower than `WideStepCandidate`

## Current state matrix

Current short remote benchmark state:

| Path | Correct | Approx. ns/Fp6 | Notes |
|---|---:|---:|---|
| `Wide` | No | `~1686` | frozen fast reference; target performance class |
| `WideSafe` | Yes | `~3083` | trusted correctness anchor |
| `WideStepCandidate` | Yes | `~1936` | current best correct path |
| `WideStepCandidateRelaxed2` | Yes | `~2228` | two-step deferred checkpoints with `2,1,1` |
| `WideStepCandidateHybrid` | Yes | `~2149` | mixed checkpoint strategy |
| `WideStepCandidateFinalFast` | No | n/a | fails validation immediately |

This is the practical target statement:

- correctness must be preserved
- the final path should approach the frozen fast reference
- the explicit target is the current fast-path class, about `1700 ns/Fp6`

## Pass/fail status

Currently passing:

- `WideSafe`
- `WideStepCandidate`
- `WideStepCandidateRelaxed2`
- `WideStepCandidateHybrid`
- primitive v2 prechecks:
  - `lufp_v2_sub2`
  - `v2_signed_carry_range`
  - `v2_sign_step<1>`
  - `v2_sign_step<2>`
  - `v2_sign_step<4>`

Currently failing:

- `Wide`
- `WideStepCandidateFinalFast`
- earlier versions of `WideStepCandidateRelaxed2` that used weaker repair than
  `2,1,1`

## Core facts already established

### 1. The core fast-style primitive ops are not the problem

We added primitive-level benchmark-local prechecks for:

- `lufp_v2_sub2`
- `v2_signed_carry_range`
- `v2_sign_step<1>`
- `v2_sign_step<2>`
- `v2_sign_step<4>`

These pass.

So the issue is not that the low-level v2 AVX building blocks are fundamentally
wrong.

### 2. The raw `Fp6` combine formulas are not the problem

We proved:

- `prepare_lufp_v2_fp6_raw(...)` is correct
- `raw + old HiNeg4Pos4 normalize + redc` matches canonical `Fp6Pack8x1::mul`

So the bug is not in the algebraic `Fp6` schedule itself.

### 3. The main problem is the renorm contract

The fast path uses a final `bias + carry + 4,2,1,1` style repair.

That shape is fast, but currently not correct.

The working candidate paths show that correctness depends on **when** we repair
and **what repair class** we use at each checkpoint.

## What we think is wrong with the fast path

The current best hypothesis is:

- the fast path is assuming that some deferred intermediate or final states are
  already in the right bound class for the current `4,2,1,1` repair
- but that assumption is not actually established by the implementation

More concretely:

- the fast path defers repair aggressively
- the intermediate states drift more than the fast renorm expects
- the current final fast renorm is therefore being applied to states outside its
  true domain

This is why:

- the fast path is fast
- but fails validation

The problem no longer looks like:

- a typo in the core schedule
- or a basic failure of `v2_sign_step`

It looks like:

- a bound/repair-class mismatch

## The current working version

The strongest correct candidate today is:

- `WideStepCandidate`
- renorm after every add/sub step
- use `ExactAdd`-style repair at each checkpoint

This path is correct locally and remotely.

Recent short remote numbers were approximately:

- `Wide` broken fast: `~1686 ns/Fp6`
- `WideSafe`: `~3083 ns/Fp6`
- `WideStepCandidate`: `~1936 ns/Fp6`

So the working candidate is:

- far faster than the old safe path
- but still not as fast as the broken fast path

That makes it the best current correctness anchor for continued optimization.

## What the deferred-checkpoint experiments proved

We tested several relaxation strategies starting from the correct
`WideStepCandidate`.

### Renorm after every step

This passes.

That means the step-style renorm machinery can produce correct results for this
kernel if it is applied often enough.

### Renorm after every two steps with `ExactAdd`

This failed.

### Renorm after every two steps with repeated `ExactAdd`

This still failed.

That means:

- the two-step deferred checkpoint is no longer in the `ExactAdd` repair class

### Renorm after every two steps with `step<2>, step<1>, step<1>`

This passes.

That is a major result.

It means:

- after two deferred ops, the checkpoint needs stronger repair than `1` or `1,1`
- a `2,1,1`-style step repair is sufficient for correctness

It also means:

- the problem is not “we can never defer”
- the problem is “we must use the right repair class for the amount of drift”

### Hybrid relaxation

We also tried a hybrid:

- weaker repair where earlier experiments suggested it might be enough
- stronger repair only where needed

This was correct, but slower than the simple every-step candidate.

So a more complicated checkpoint schedule has not yet produced a win.

## Checkpoint map

The useful checkpoint-level picture so far is:

| Checkpoint style | Known result | Current interpretation |
|---|---|---|
| one deferred add/sub step | `ExactAdd` works | this remains in the `ExactAdd` repair class |
| two deferred steps on `c0`-style path | weaker repair can sometimes work | `c0` is more tolerant |
| two deferred steps on `c1.re` / `c2.re` | `ExactAdd` and repeated `ExactAdd` fail | these leave the `ExactAdd` class |
| two deferred steps with `2,1,1` | passes | this is currently the right repair class for deferred-two checkpoints |
| final-fast-only swap | fails immediately | the fast final renorm is not a valid drop-in boundary for the working candidate |

The key concrete deferred checkpoints are:

- `w0 - t1`
- `w0 - t1 - t2`
- `w1 - t0 - t1`
- `w2 - t0 - t2`
- post-`mul_by_u_plus_1` states

Current understanding:

- `w0 - t1` was the first observed bad state in the original failing step path
- once that was repaired with `ExactAdd`, the first major integration divergence
  moved later
- `c1.re` and `c2.re` are the checkpoints that force the stronger deferred
  repair

## Relationship to the bound worksheet

See:

- [fp6-mul-bound-worksheet.md](/c/Users/peter/repo/blst/docs/fp6-mul-bound-worksheet.md)

That worksheet derives the raw lazy-combine classes for the baseline `Fp6`
schedule and concludes that the narrowest common final source class is:

- `LNeg4_Pos4_HNeg4_Pos4`

This is important, but it is easy to misread.

What the worksheet proves:

- if we preserve the raw schedule
- and look at the final raw `c0/c1/c2` outputs
- then one final repair kernel should only need to accept that common class

What the worksheet does **not** prove:

- that every intermediate deferred checkpoint is in the `ExactAdd` class
- that the current fast implementation of the final repair is correct
- that the implementation details of low/high carry, sign mapping, and masked
  repair exactly match the abstract bound model

So the worksheet still appears useful and likely correct in spirit.

The current best interpretation is:

- the worksheet supports the existence of a fast final-repair design
- but the current fast path is not yet that correct design

One especially important consequence for finishing the work:

- the worksheet supports the belief that a final correct fast-shape kernel
  should exist
- therefore the current gap is more likely an implementation-level repair-class
  mistake than a fundamentally wrong optimization goal

In particular, the two-step deferred experiment strongly suggests:

- single-step checkpoints behave like `ExactAdd`
- two-step checkpoints behave like a wider class that needs `2,1,1`

That is consistent with the worksheet’s general message that the kernel leaves
the strict narrow world quickly and must be repaired according to actual drift.

## Why we want to keep the fast path frozen

The broken fast path is valuable even though it is wrong.

It gives us:

- a stable performance reference
- a known-good codegen shape
- a target end-state structure

If we keep mutating it directly, we lose:

- the ability to measure whether a “fix” is actually converging on the original
  fast shape
- the ability to compare object shape and performance against the known fast
  version

So the correct workflow is:

1. keep `Wide` frozen
2. keep `WideSafe` frozen as the correctness reference
3. improve candidate paths until one is both correct and close to `Wide`

This should be treated as a hard requirement:

- do not “fix” the problem by discarding the fast-shape target
- the job is not finished unless the final correct path is back in the frozen
  fast-path performance class, approximately `1700 ns/Fp6`

## Where we are currently stuck

The main unresolved problem is:

- we still do not have a precise implementation-level bound model for the
  deferred checkpoint states

In other words:

- we know some checkpoint schedules work
- we know some do not
- but we do not yet know a clean rule like:
  - “this checkpoint is in class X”
  - “therefore repair sequence Y is sufficient”

This is the current bottleneck.

We are no longer blocked by:

- core primitive correctness
- the outer `Fp6` algebra

We are blocked by:

- not knowing how much growth each deferred checkpoint actually admits
- not knowing the minimal correct repair class for each checkpoint pattern

In short:

- we now know enough to keep making progress
- but we do not yet have the bound model needed to finish the job mechanically
- the remaining work is checkpoint classification and repair selection

## Practical current conclusions

1. `Wide` is still the performance target.

2. `WideSafe` is still the correctness anchor.

3. `WideStepCandidate` is currently the best correct path.

4. A two-step checkpoint is **not** in the `ExactAdd` class.

5. A two-step checkpoint **is** repairable by `2,1,1`.

6. Correctness is therefore compatible with step-style repair, but the minimal
   cheap checkpoint schedule is still unknown.

## Most promising next directions

These are the most useful next steps.

### 1. Derive checkpoint-specific classes

Do not reason only from the final raw envelope.

Instead, classify the actual deferred checkpoints used by candidate schedules,
for example:

- `w0 - t1`
- `w0 - t1 - t2`
- `w1 - t0 - t1`
- `w2 - t0 - t2`
- post-`mul_by_u_plus_1` states

The goal is to identify which checkpoints are:

- `ExactAdd`
- `2,1,1`
- or something wider

### 2. Use the working candidate as the main optimization base

The best current base is not `WideSafe`.

It is:

- `WideStepCandidate`

That already proves the step-style repair model works.

### 3. Only relax where a checkpoint class supports it

The failed experiments show that blind batching is inefficient.

The next relaxation should be driven by checkpoint class, not by aesthetics or
 “every two everywhere”.

### 4. Keep the object/codegen viewpoint

The frozen fast path is useful not only for timing, but for code shape.

A correct candidate that converges toward the fast path should also converge in:

- live-range structure
- spill behavior
- helper boundaries
- renorm shape

### 5. Avoid work that is not in the intended end-state

We already learned this the hard way.

The end-state is step-style repair, not a return to the old safe normalizer.

So exploration should stay aligned with:

- carry
- sign
- step repairs
- deferred-checkpoint classes

## Short summary

The project is in a good debugging state now.

We know:

- the fast path is broken because its repair assumptions are too optimistic
- the safe path is correct but too slow
- the step-style repair model can be correct
- two-step deferred checkpoints require stronger repair than `ExactAdd`
- `2,1,1` works for those two-step checkpoints

We do **not** yet know:

- the minimal correct checkpoint schedule that reproduces the fast-path shape
- the exact bound class of each deferred intermediate

The best current finishing strategy is:

1. keep `Wide` frozen as the `~1700 ns/Fp6` requirement target
2. keep `WideSafe` frozen as the correctness anchor
3. keep using `WideStepCandidate` as the working optimization base
4. relax only one checkpoint class at a time
5. require both:
   - validator pass
   - convergence back toward the `~1700 ns/Fp6` fast-path class

So the remaining problem is no longer broad.

It is:

- derive or experimentally pin down the checkpoint growth classes
- then rebuild the fast path from the correct candidate using only justified
  relaxations
