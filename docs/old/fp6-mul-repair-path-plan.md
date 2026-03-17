# `Fp6 mul` Repair Paths Plan

This note describes the two implementation paths we want to compare for `Fp6`
multiply repair.

It is a planning document only.

Non-goal:

- do **not** modify the existing broken fast path that currently measures around
  `~1700 ns/Fp6`
- treat that path as a performance reference and a repair-shape hint only

Related notes:

- [fp6-mul-debug-explainer.md](/c/Users/peter/repo/blst/docs/fp6-mul-debug-explainer.md)
- [fp6-mul-bound-worksheet.md](/c/Users/peter/repo/blst/docs/fp6-mul-bound-worksheet.md)
- [fp6-mul-sparse5b-bounds.md](/c/Users/peter/repo/blst/docs/fp6-mul-sparse5b-bounds.md)

## Goal

Compare two repair strategies for the widened `Fp6 mul` combine phase:

1. `ExactAdd`-driven immediate repair
2. delayed repair with a larger admissible source class and a staged sign-step
   ladder

The comparison is about:

- correctness from first principles
- admissible source bounds
- likely cost shape
- implementation risk

## Shared First-Principles Model

Track each final `Fp` half independently.

Use:

- `W` for one strict widened `Fp` term
- `C(P,N)` for an exact signed sum of `P` positive `W` terms and `N` negative
  `W` terms

After `v2_carry_from_biased(...)`, a `C(P,N)` value has:

- low canonicalized into `[0, R)`
- high residue in `[-N*p, P*p)`

This gives the dual view we want:

- linear-combination provenance: `C(P,N)`
- repair-facing envelope: `[-N*p, P*p)`

The current no-intermediate-repair block envelopes are:

- `c0`: `[-4*p, 4*p)`
- `c1`: `[-3*p, 3*p)`
- `c2`: `[-2*p, 2*p)`

## Path 1: Immediate `ExactAdd`

### Idea

Repair as soon as a checkpoint enters the one-step signed class that
`ExactAdd` can handle.

At code level, `ExactAdd` is:

- low carry
- sign test on the top limb
- one masked add of `p` for negative lanes

So it is a one-step repair primitive, not a general-purpose renorm.

### Admissible Class

This path is valid only when the carried value is already in the one-step
domain:

- roughly `[-p, p)`

That means:

- subtraction from a repaired `W` state is admissible
- positive-only overflow such as `W + W` is **not** directly an `ExactAdd`
  checkpoint
- two unrepaired subtractions in a row are already outside the intended class

### Expected Schedule Shape

Under the current first-principles interpretation, the admissible checkpoint
schedule is:

- `c0_sub_t1`
- `c0_sub_t2`
- `c1_sub_t0`
- `c1_sub_t1`
- `c2_sub_t0`
- `c2_sub_t2`

This is the `0x333` mask in the current sparse-checkpoint encoding.

No checkpoint repair should be spent on:

- `mul_by_u_plus_1` outputs directly
- final positive adds

unless a separate proof shows those states still fit the one-step class.

### Advantages

- simplest proof story
- reuses existing trusted primitive
- easiest path to a fully correct checkpoint schedule

### Disadvantages

- likely too many repairs
- probably leaves a material performance gap versus the fast-path target
- still needs a final consumer-legal output contract

## Path 2: Delayed Large-Class Repair

### Idea

Do not repair at each checkpoint.

Instead:

1. let the block drift within a larger bounded source class
2. bias by `k*p` to center the carried high range inside `[-2^m*p, p)`
3. apply `O(m)` signed repair steps

This is the general delayed-repair strategy.

### Canonical Example

For a carried value in `[0, 2*p)`:

1. subtract `p`
2. now the value is in `[-p, p)`
3. perform one sign-based repair step

So positive overflow such as `W + W` is not a problem by itself. It just
belongs to a different repair class than plain `ExactAdd`.

### Block Targets

Using the current end-only per-`Fp` derivation:

- `c0` needs a final repair family for `[-4*p, 4*p)`
- `c1` needs a final repair family for `[-3*p, 3*p)`
- `c2` needs a final repair family for `[-2*p, 2*p)`

This suggests one of two designs:

- one common ladder for the largest class `[-4*p, 4*p)`
- or a family of ladders where `c2` uses a cheaper final repair than `c0`

### Admissible Repair Contract

For a carried source envelope `[-A*p, B*p)`, choose:

- a bias `k*p` such that the biased range is inside `[-2^m*p, p)`
- then a sign-step ladder that is proved correct for that range

This is the contract that must be proved against the actual renorm primitives in
[Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h).

### Advantages

- matches the intended lazy-repair philosophy
- naturally explains both negative drift and positive overflow
- offers a plausible route toward the fast-path cost profile

### Disadvantages

- proof burden is higher
- requires a new repair-family contract, not just checkpoint bookkeeping
- may need different ladders for `c0`, `c1`, and `c2`

## Comparison

`Path 1` is the low-risk baseline:

- easiest to prove
- easiest to implement correctly
- least likely to hit the performance target

`Path 2` is the high-upside path:

- more aligned with the observed fast-path performance ceiling
- more natural for the current end-only bounds
- requires a real repair-family proof before code changes are trustworthy

## Recommended Order

1. Keep the existing broken fast path unchanged.
2. Treat `Path 1` as the correctness baseline and checkpoint canary.
3. Treat `Path 2` as the real optimization target.
4. First prove a delayed final repair ladder for:
   - `[-2*p, 2*p)`
   - `[-3*p, 3*p)`
   - `[-4*p, 4*p)`
5. Only after that, decide whether some selective intermediate repairs are still
   cheaper than a single stronger final ladder.

## Concrete Next Steps

1. Write a code-level contract note for the repair primitives:
   - what source envelope each current primitive actually accepts
   - what target class it guarantees

2. Derive a candidate delayed final ladder for:
   - `c2` first on `[-2*p, 2*p)`
   - then `c1`
   - then `c0`

3. Benchmark the delayed final ladder in isolation before integrating it into
   `Fp6 mul`.

4. Compare against the immediate-`ExactAdd` schedule and the untouched broken
   fast path.

## Decision Rule

Choose the delayed strategy only if both are true:

- the repair-family proof is explicit and tied to the actual code-level steps
- the benchmark win over the immediate-`ExactAdd` path is material

Otherwise, keep the immediate path as the correctness baseline and continue
using the broken fast path only as a performance reference.
