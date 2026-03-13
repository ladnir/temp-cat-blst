# Fp12 LUFp V2 Proposal

## Audience

This document is for someone new to this repository who needs enough context to work effectively on the next `Fp12` multiplication design.

It is meant to do four jobs:

1. explain the relevant field tower structure
2. explain how the successful `Fp6` design works
3. explain what is and is not done yet for `Fp12`
4. give a concrete plan for a first correct `Fp12` implementation

## Executive Summary

The only completed optimization result so far is at the `Fp6` layer.

In benchmark `LUFp v2` code, `Fp6` multiplication was improved from roughly `2200 ns/Fp6` to roughly `1660 ns/Fp6` by deriving the final one-sided block bounds and implementing the fused strict-negative repair.

That result was then improved again by upstreaming the seeded first-step entry
into the live fused helper. The current live measurement is about
`1611 ns/Fp6`, with the broken historical reference still around
`1567 ns/Fp6`.

Nothing analogous has been implemented for `Fp12` yet.

This document is therefore not a report on an existing `Fp12` implementation. It is a design brief for how to start one, using the same reasoning pattern that worked for `Fp6`.

## Goal

Reuse the lazy widened design ideas that won for `Fp6`, but apply them correctly at the `Fp12` layer.

The intended pattern is:

1. compute raw widened values with cheap linear operations
2. stop at a well-defined post-carry representation
3. derive explicit bounds for the next block
4. apply one delayed repair at the latest safe point

Status today:

- `Fp6` has a real measured win
- `Fp12` does not yet have a new lazy implementation
- no new `Fp12` bounds proof exists yet
- the main design question is what lazy `Fp6` interface `Fp12` should consume
- the likely answer is a new post-carry lazy `Fp6` interface, not the current terminal repaired `Fp6` winner output

This distinction matters:

- `Fp6` is done enough to benchmark and use as a reference design
- `Fp12` is still at the proposal stage

## What Delayed Repair Means

“Delayed repair” is the optimization strategy used throughout this work.

The ordinary approach to modular arithmetic is:

1. do a small amount of arithmetic
2. immediately repair or reduce the result back toward `[0, p)` or another tight class
3. continue

Delayed repair does something different:

1. do several arithmetic operations in a widened representation
2. allow the value to drift away from canonical range by a bounded number of multiples of `p`
3. keep track of that drift symbolically
4. repair only when a consumer actually needs a tighter bound

So “delayed” does not mean “ignore correctness until the end.” It means:

- we permit temporary drift
- but only inside a proved envelope
- and we choose the latest safe point where one repair can cover several earlier operations

Why this can be faster:

- repair is not free
- if you repair after every subtract or add, you spend a lot of instructions re-tightening the same value repeatedly
- if you can prove that a block of operations only grows the value to, for example, `[-4p, 4p)`, then you may be able to do one larger repair at the end instead of several smaller repairs in the middle

In this project, delayed repair usually looks like:

```text
raw linear operations
-> post-carry signed interval like [-A*p, B*p)
-> one bias into a strict-negative class
-> a short keep-negative repair ladder
-> final lazy or canonical class
```

Example at the `Fp6` layer:

- the winning `Fp6` kernel does not repair after every intermediate step
- it forms the whole `c0`, `c1`, and `c2` block first
- then it repairs each final block once, using a schedule matched to that block’s proved interval

That is the pattern we want to repeat for `Fp12`, but using a composable post-carry `Fp6` interface.

## Field Tower Primer

This code is working over the usual BLS12-381 extension tower:

```text
Fp2  = Fp[u]   / (u^2 + 1)
Fp6  = Fp2[v]  / (v^3 - (u + 1))
Fp12 = Fp6[w]  / (w^2 - v)
```

Operationally:

- `Fp2` elements are pairs `(re, im)`
- `Fp6` elements are triples `(c0, c1, c2)` of `Fp2`
- `Fp12` elements are pairs `(c0, c1)` of `Fp6`

Two very important linear maps appear repeatedly:

```text
(u + 1) * (x.re, x.im) = (x.re - x.im, x.re + x.im)
v * (x0, x1, x2)       = ((u + 1) * x2, x0, x1)
```

You do not need deep algebraic background to work on this optimization, but you do need to recognize that:

- `Fp12` multiplication is mostly three `Fp6` multiplications
- the remaining work is linear
- delayed repair is most valuable precisely in that situation

## Terminology

### `UFp`

Widened field element. Arithmetic has been performed, but the value is not yet fully reduced to canonical `[0, p)`.

### `LUFp`

Lazy widened field element. Same mathematical value, but we intentionally allow bounded drift so we can do fewer repairs.

### carry

In this document, “carry” means converting the widened limb representation into a stable signed post-carry form that is easy to reason about.

The important split is:

- pre-carry: easy for raw arithmetic, awkward for proofs
- post-carry: right abstraction for repair contracts

### repair / renorm

A repair step adds or subtracts multiples of `p` so the value lands in a bounded class acceptable to the next consumer.

### terminal lazy class

A lazy representation that is valid for immediate final reduction, but not yet proved valid for more widened arithmetic.

### composable lazy class

A lazy representation whose bounds are explicit enough that further adds, subtracts, and nonresidue maps can be reasoned about safely.

This distinction is the core design issue for `Fp12`.

## Code Map

These are the files that matter for this work.

### Canonical tower code

- [Fp12_avx.cpp](/c/Users/peter/repo/blst/blst/fields/Fp12_avx.cpp)
  - canonical packed `Fp12` multiply
  - current widened `Fp12` entry point is at [Fp12_avx.cpp](/c/Users/peter/repo/blst/blst/fields/Fp12_avx.cpp#L312)
- [Fp6_avx.cpp](/c/Users/peter/repo/blst/blst/fields/Fp6_avx.cpp)
  - canonical packed `Fp6` multiply
  - widened `Fp6` multiply is at [Fp6_avx.cpp](/c/Users/peter/repo/blst/blst/fields/Fp6_avx.cpp#L200)

### Benchmark-only `LUFp v2` code

- [Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h)
  - the promoted `Fp6` winner lives here
  - current live `Fp6` winner starts at [Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L1268)
  - the generic fused strict-negative renorm helper is at [Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L830)
- [Fp_avx_lufp_v2_fp6_bench.cpp](/c/Users/peter/repo/blst/benchmark/Fp_avx_lufp_v2_fp6_bench.cpp)
  - benchmark and validation harness for the `Fp6` work
- [Fp12_tower_avx_bench.cpp](/c/Users/peter/repo/blst/benchmark/Fp12_tower_avx_bench.cpp)
  - `Fp12` tower benchmarks
  - a good place for future benchmark-only `Fp12` experiments

### Existing documentation

- [fp6-mul-delayed-block-final-proof.md](/c/Users/peter/repo/blst/docs/fp6-mul-delayed-block-final-proof.md)
  - proof-style writeup for the current `Fp6` winner
- [fp6-mul-sparse5b-bounds.md](/c/Users/peter/repo/blst/docs/fp6-mul-sparse5b-bounds.md)
  - earlier bounds exploration and counterexample history

## How Canonical Fp12 Mul Works Today

The canonical widened packed `Fp12` multiply is in [Fp12_avx.cpp](/c/Users/peter/repo/blst/blst/fields/Fp12_avx.cpp#L312).

Its structure is:

```text
t0 = a0 * b0
t1 = a1 * b1
m  = (a0 + a1) * (b0 + b1)

r1 = m - t0 - t1
r0 = t0 + v * t1
```

The corresponding code does exactly that:

- `t0 = a0*b0`
- `t1 = a1*b1`
- `r.c[1] = (a0+a1)(b0+b1) - t0 - t1`
- `r.c[0] = t0 + v*t1`

The important engineering fact is that `Fp12` multiply is mostly three `Fp6` multiplies plus a small amount of linear glue.

That is why `Fp6` performance directly matters for `Fp12`.

## Reference: Current Fp6 Winner

The promoted `Fp6` kernel in [Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L1268) computes the three Karatsuba blocks and then applies the fused strict-negative repair independently to each final `Fp` half.

The live source classes are:

- `c0.re in [-3p, 4p)`
- `c0.im in [-4p, 3p)`
- `c1.re in [-3p, 2p)`
- `c1.im in [-2p, 3p)`
- `c2.re in [-2p, 2p)`
- `c2.im in [-2p, 2p)`

At a high level, the winner does this:

1. compute the three base `Fp2` products `t0`, `t1`, `t2`
2. compute the three Karatsuba cross-products
3. form the raw `c0`, `c1`, `c2` blocks using only linear ops
4. bias each `Fp` half into a strict-negative class
5. run the fused repair ladder `4 -> 2 -> 1 -> <1,1>` or `2 -> 1 -> <1,1>`,
   depending on the proved per-half source class

One implementation detail matters for future work: the live generic helper no
longer always materializes `x - B*p` and then applies the first repair step
separately. When `B <= bit_floor(A + B)`, it starts from the equivalent seeded
state:

- `x + (bit_floor(A + B) - B) * p`

That "seed4" idea was first benchmarked as a separate `Fp6` candidate and then
upstreamed into the generic `lufp_v2_renorm_negative_fold_fused<A, B>` helper.
The proof contract stayed the same; only the entry state changed.

Operationally, the current winner returns a widened `UFp6` value with all six `Fp` halves in:

- `[0, p)`

This is a successful design for `Fp6`.

But there is a subtlety:

- these envelopes are enough for immediate `UFp6Pack8x1::redc`
- they are not yet enough to justify arbitrary later widened arithmetic

That is exactly the question `Fp12` has to answer before implementation starts.

## Why Fp6 Success Does Not Automatically Lift To Fp12

It is tempting to say:

1. the `Fp6` winner is correct
2. its output is bounded
3. therefore `Fp12` can consume that output lazily

That reasoning is incomplete.

What is actually proved today is weaker:

- the current repaired `Fp6` winner produces a widened value acceptable to `UFp6Pack8x1::redc`

What is not yet proved:

- that the same output class is stable under more widened `UFp6` adds, subs, and `mul_by_nonresidue`

So the current `Fp6` winner output is a terminal lazy class, not a composable lazy class.

That single distinction explains why `Fp6` success does not automatically give an `Fp12` design.

## Why The Current Fp6 Winner Is Still The Wrong Fp12 Interface

The current `Fp6` winner returns all six `Fp` halves in:

```text
[0, p)
```

That is excellent for standalone `Fp6` multiplication because it means:

- the widened result is already in a canonical terminal class
- `UFp6Pack8x1::redc` can run immediately
- the kernel is easy to validate against the ordinary `Fp6` product

But for `Fp12`, that same fact is a warning sign.

If one tries to build `Fp12` directly from that terminal interface, the bound
story collapses back toward ordinary canonical arithmetic:

- `r1 = m - t0 - t1` becomes a difference of three `[0, p)` terms
- `r0 = t0 + v * t1` becomes a sum of `[0, p)` terms and small linear maps of
  `[0, p)`

Those bounds are valid, but they do not expose the larger lazy drift class that
actually made the `Fp6` kernel fast. In other words:

- the current `Fp6` winner output is too repaired to be an interesting `Fp12`
  building block
- it is a terminal class, not a composable lazy class

That is why the correct `Fp12` starting point is still a new post-carry
interface:

- after the raw `Fp6` combine
- after signed carry
- before the final strict-negative repair

That post-carry boundary is where the useful `Fp12` proof work belongs.

## Current State

What is real today:

- the `Fp6` fused strict-negative design is implemented
- the `Fp6` fused strict-negative design is benchmarked
- the `Fp6` fused strict-negative design is documented and proved at the kernel level
- the practical result is roughly `2200 ns/Fp6 -> 1660 ns/Fp6 -> 1611 ns/Fp6`
- the old delayed-centered path is kept only as a broken benchmark reference

What is not real today:

- a new lazy `Fp12` implementation
- a proved post-carry lazy `Fp6` interface for `Fp12`
- a validated `Fp12` delayed-final repair schedule
- any measured `Fp12` win from this design

So this document should be read as a proposal, not as implementation history.

## The Correct Next Interface

The next attempt should be organized around a new benchmark-only helper:

```text
mul_fp6_pack_lufpv2_wide_postcarry(...)
```

Its contract should be:

- compute the same raw `Fp6` result as the current winner
- stop after post-carry normalization
- do not apply the terminal fused strict-negative repair used for immediate `UFp6Pack8x1::redc`
- return explicit per-half post-carry envelopes

This is the interface `Fp12` likely needs.

Why:

- post-carry is the right proof abstraction
- post-carry classes compose through linear `Fp12` operations
- terminal repaired classes only guarantee immediate reduction safety

## Implementation Roadmap

This is the order of work I would recommend to a new contributor.

### 1. Write down the post-carry `Fp6` interface explicitly

Create a short note that states, per `Fp` half:

- the exact raw `C(P,N)` class
- the post-carry interval `[-A*p, B*p)`
- which operations are allowed before another repair

Do not start from the repaired `Fp6` winner output. Start from the post-carry state before its terminal repair.

### 2. Add a benchmark-only `Fp6` post-carry kernel

Add a helper in [Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h) that:

- computes the current raw `Fp6` combine
- runs post-carry
- stops there

This helper should not replace the current `Fp6` winner. It is an internal building block for `Fp12` experiments.

### 3. Derive the `Fp12` envelopes from that interface

Repeat exactly the proof pattern that worked for `Fp6`:

- write `r0` and `r1` as linear combinations of post-carry `Fp6` terms
- derive per-half interval bounds
- choose the smallest delayed repair schedule that covers those bounds

This is where the proof work belongs.

### 4. Implement only after the bounds are written down

Once the intervals are explicit:

- add a benchmark-only `mul_fp12_pack_lufpv2_wide(...)`
- validate it against canonical `Fp12Pack8x1::mul`
- only then benchmark it remotely

### 5. Compare three things, not two

The useful comparison set is:

- current canonical `Fp12`
- legacy `Fp12`
- new `LUFp v2` `Fp12`

That tells you whether the new path is actually beating the current baseline or only beating older code.

## Expected Comparison

The baseline `UFp` `Fp12` multiply in [Fp12_avx.cpp](/c/Users/peter/repo/blst/blst/fields/Fp12_avx.cpp#L312) pays for three canonical widened `Fp6` multiplies and then does the `Fp12` linear combine.

The intended `LUFp v2` path would replace those three `Fp6` multiplies with a composable post-carry `Fp6` kernel, keep the `Fp12` combine widened, and repair only once per final `Fp6` block at the `Fp12` level.

If this works, the expected win comes from:

- faster widened `Fp6` multiplication
- reusing lazy widened `Fp6` drift classes instead of canonicalizing early
- moving the next repair step to the latest safe point

As of the latest remote measurements, the relevant current baselines are:

- `Fp6 LUFp v2 Wide`: about `1611 ns/Fp6`
- `Fp6 LUFp v2 WideBroken`: about `1567 ns/Fp6`
- current canonical `Fp12`: about `7165 ns/Fp12`
- legacy `Fp12`: about `7224 ns/Fp12`

So there is room, but `Fp12` will only benefit if the lazy `Fp6` interface is actually composable.

## What We Learned After Promoting The Fp6 Winner

Two implementation findings are relevant for the eventual `Fp12` work.

### 1. The real win was the seeded fused entry, not a new tower formula

The measured `Fp6` improvement after the first fused strict-negative design came
from shortening the renorm entry path:

- raw `Fp6` combine stayed algebraically unchanged
- the helper started directly from the post-first-step candidate when the first
  repair scale covered the strict-negative bias
- that reduced raw renorm cost and survived into full `Fp6` cost

So the next `Fp12` work should keep the same mindset:

- do not change the tower algebra without a strong reason
- look first at the delayed repair boundary and its entry state

### 2. The sign-mask codegen suspicion did not produce a new win

Object inspection of the live `Fp6` winner showed GCC scalarizing parts of the
sign-mask path inside the very large inlined multiply symbol. A native-mask
helper experiment was tried specifically to recover `vpmovq2m`-style codegen in
that hot symbol.

That experiment did not change:

- object size
- sign-mask instruction shape
- or benchmark numbers in a meaningful way

So there is no new promoted optimization result from that line of work today.
For `Fp12`, the actionable lesson is simply:

- prefer concrete proof/benchmark wins over codegen speculation
- verify any assembly-motivated change with both object inspection and end-to-end timing

## Validation Plan

1. Random packed validation against `Fp12Pack8x1::mul`.
2. Benchmark split:
   - full reduced `Fp12`
   - wide-only `Fp12`
   - redc-only `Fp12`
3. Compare:
   - canonical `Fp12`
   - legacy `Fp12`
   - new `LUFp v2` `Fp12`

There is no new `Fp12` path to validate yet. The purpose of this section is to define the validation bar for the first implementation.

## Final Guidance

If you are new to the project, the most important thing to internalize is this:

- `Fp6` was won by delaying repair to the end of the `Fp6` block
- the current concrete win is roughly `2200 ns/Fp6 -> 1660 ns/Fp6 -> 1611 ns/Fp6`
- `Fp12` has not been redesigned yet
- the current task is to build that redesign from the same proof pattern, starting from a composable post-carry interface

That is the accurate state of the project today.
