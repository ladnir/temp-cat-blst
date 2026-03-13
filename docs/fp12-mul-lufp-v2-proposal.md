# Fp12 LUFp V2 Proposal

## Audience

This document is for someone new to this repository who needs enough context to work effectively on the next `Fp12` multiplication design.

It is meant to do four jobs:

1. explain the relevant field tower structure
2. explain how the successful `Fp6` design works
3. explain what is and is not done yet for `Fp12`
4. give a concrete plan for a first correct `Fp12` implementation

## Executive Summary

This document started as a proposal. It is now a mixed design note and status report.

The `Fp6` `LUFp v2` winner is real and remains the main success story:

- core-library wide `Fp6` (`UFp`, type 2) is about `2203 ns/Fp6`
- `LUFp v2` with terminal `Fp6` renorm (`mul_fp6_pack_lufpv2_wide<true>`) is about `1599 ns/Fp6`
- the raw lazy type-3 kernel (`mul_fp6_pack_lufpv2_wide<false>`) is about `1328 ns/Fp6`

There is now also a correct benchmark-only `Fp12` `LUFp v2` path.

At the `Fp12` layer there are currently three relevant variants:

1. `basic`
   - the core-library packed `Fp12` multiply
   - internally this is already type-2 / `UFp` wide arithmetic
2. `Fp6<true>`
   - the `Fp12` formula with the three inner `Fp6` multiplies replaced by the live `LUFp v2` `Fp6` winner
   - each inner `Fp6` multiply ends with its own terminal renorm
3. `Fp6<false>`
   - the true delayed-final-repair `Fp12` path
   - the three inner `Fp6` multiplies use the same kernel with renorm disabled
   - `Fp12` then performs lazy type-3 linear ops and renorms only the final `Fp12` halves

The important current result is:

- both `Fp6<true>` and `Fp6<false>` beat `basic`
- `Fp6<false>` is correct
- but `Fp6<false>` is not yet the fastest of the three
- today, `Fp6<true>` is still slightly faster than `Fp6<false>`

So `Fp12` is no longer at the purely hypothetical stage. The open problem has changed:

- not "can we make a correct delayed-final-repair `Fp12` path at all?"
- but "why is the correct delayed-final-repair path still a little slower than the `Fp6<true>` hybrid path?"

## Goal

Reuse the lazy widened design ideas that won for `Fp6`, but apply them correctly at the `Fp12` layer and understand the remaining gap between the current correct delayed-final-repair path and the current fastest hybrid path.

The intended pattern is:

1. compute raw widened values with cheap linear operations
2. keep them in a provable lazy representation
3. derive explicit bounds for the next block
4. apply one delayed repair at the latest safe point

Status today:

- `Fp6` has a real measured win
- `Fp12` now has two benchmark-only `LUFp v2` implementations
- the delayed-final-repair `Fp12` path is correct against canonical `Fp12Pack8x1::mul`
- explicit final-half `Fp12` bounds exist for that path
- the remaining problem is performance, not first correctness

This distinction matters:

- `Fp6` is done enough to benchmark and use as a reference design
- `Fp12` is no longer purely a proposal
- but the final delayed-final-repair design is not yet the performance winner

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

That is the pattern we want to repeat for `Fp12`, but now in the concrete `Fp6<false>` form that keeps the current `Fp6` kernel lazy through the `Fp12` linear stage and repairs only at the end.

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

There are three distinct notions of "reduction" in this work. Keeping them separate is essential.

### `Fp` (type 1)

Ordinary reduced field element.

- represented in the normal 6-limb radix-48 form
- logically reduced modulo `p`
- acceptable as an ordinary `Fp` input or output

### `UFp` (type 2)

Wide field element.

- double-wide limb storage
- still logically reduced modulo `p`
- after each linear operation, the code re-establishes the type-2 contract with carry and conditional `p` shifts
- this is the model used by the core-library wide `Fp6`/`Fp12` code

### `LUFp` (type 3)

Lazy wide field element.

- also double-wide limb storage
- may spill beyond the 48-bit limb contract before later repair
- not yet logically reduced modulo `p`
- the whole point is to defer the type-2 style repair work and batch it later

### carry

In this document, “carry” means radix carry normalization, not full modular repair.

For a double-wide radix-48 pack:

- after carry, limbs `0..14` are back in the normal 48-bit range
- the top limb still carries the signed drift
- this is a representation cleanup step
- it is not the same as restoring a type-1 or type-2 modular contract

### repair / renorm

A repair step adds or subtracts multiples of `p` so the value lands in a proved target class.

Examples:

- type-2 repair after a `UFp` add/sub
- strict-negative fused renorm in `LUFp v2`
- final `UFp12Pack8x1::redc`

### wide-only

Benchmark term used throughout this file.

- benchmark only the widened multiply kernel
- do not include final reduction to ordinary `Fp12`

### full multiply

Benchmark term used throughout this file.

- benchmark widened multiply
- then benchmark the final reduction to ordinary `Fp12`

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
  - current live `Fp6` winner template starts at [Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L1379)
  - `mul_fp6_pack_lufpv2_wide<true>` is the live `Fp6` winner
  - `mul_fp6_pack_lufpv2_wide<false>` is the same kernel with terminal renorm disabled
  - current `Fp12` delayed-final-repair kernel is at [Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L1436)
  - `lufp6_v2_mul_by_nonresidue` is now a one-parameter in-place helper at [Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L1183)
  - the generic fused strict-negative renorm helper is at [Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L830)
- [Fp_avx_lufp_v2_fp6_bench.cpp](/c/Users/peter/repo/blst/benchmark/Fp_avx_lufp_v2_fp6_bench.cpp)
  - benchmark and validation harness for the `Fp6` work
  - includes both `BM_Fp6MulLufpWide` and `BM_Fp6MulLufpWideLazy`
- [Fp12_tower_avx_bench.cpp](/c/Users/peter/repo/blst/benchmark/Fp12_tower_avx_bench.cpp)
  - `Fp12` tower benchmarks
  - now contains benchmark-only `Fp12` comparisons for `basic`, `Fp6<true>`, and `Fp6<false>`

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

One correction is important here: the core-library `Fp12` path is not doing ordinary fully reduced `Fp` arithmetic internally. Its inner wide multiplies and wide linear ops are already living in the type-2 `UFp` world. So the right comparison is:

- `basic`: type-2 `UFp` wide arithmetic throughout
- `Fp6<true>`: use the live type-3 `Fp6` winner, then re-enter the type-2 `UFp6` `Fp12` glue
- `Fp6<false>`: stay in type-3 lazy form through the `Fp12` linear stage and repair only at the end

## Reference: Current Fp6 Winner

The promoted `Fp6` kernel in [Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L1379) computes the three Karatsuba blocks and then applies the fused strict-negative repair independently to each final `Fp` half.

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

Operationally, `mul_fp6_pack_lufpv2_wide<true>` returns a widened `UFp6` value with all six `Fp` halves in:

- `[0, p)`

This is a successful design for `Fp6`.

The same code body can now also be instantiated as `mul_fp6_pack_lufpv2_wide<false>`, which returns the raw type-3 lazy output before the terminal `Fp6` renorm.

That split is the entire basis of the current `Fp12` experiments.

## Current State

What is real today:

- the `Fp6` fused strict-negative design is implemented
- the `Fp6` fused strict-negative design is benchmarked
- the `Fp6` fused strict-negative design is documented and proved at the kernel level
- the practical `Fp6` result is roughly `2203 ns/Fp6 -> 1599 ns/Fp6`
- the raw lazy type-3 `Fp6` kernel is roughly `1328 ns/Fp6`
- a correct benchmark-only `Fp12` delayed-final-repair path exists
- a correct benchmark-only `Fp12` hybrid `Fp6<true>` path exists
- both `Fp12` `LUFp v2` variants beat the current `basic` baseline

What was tried and then discarded:

- a separate post-carry `Fp6` interface
- that experiment did not establish a useful new composable contract
- the current code now keeps exactly one `mul_fp6_pack_lufpv2_wide` implementation and selects the terminal renorm with a template parameter

What is still not done:

- a proof-quality explanation of why `Fp6<false>` is still slightly slower than `Fp6<true>`
- the final performance cleanup needed to make the true delayed-final-repair path the clear winner
- any promotion of the `Fp12` benchmark-only code into the library

So this document should now be read as:

- historical proposal
- current implementation status
- and a list of the specific remaining performance questions

## Current Fp12 Variants

There are now three concrete `Fp12` multiplication variants worth comparing.

### 1. `basic`

This is the core-library packed `Fp12` wide multiply in [Fp12_avx.cpp](/c/Users/peter/repo/blst/blst/fields/Fp12_avx.cpp#L312).

Properties:

- type-2 `UFp` arithmetic throughout
- inner `Fp6` multiplies are the core-library wide `UFp6` multiplies
- linear `Fp12` glue uses the existing `UFp6` operators

### 2. `Fp6<true>`

This is the benchmark-only hybrid path in [Fp12_tower_avx_bench.cpp](/c/Users/peter/repo/blst/benchmark/Fp12_tower_avx_bench.cpp#L957).

Properties:

- the three inner `Fp6` multiplies call `mul_fp6_pack_lufpv2_wide<true>`
- each inner `Fp6` multiply ends with the live fused strict-negative renorm
- after that, the path returns to the ordinary `UFp6` `Fp12` glue

This path is useful because it isolates the value of swapping in the live `Fp6` winner without redesigning the rest of `Fp12`.

### 3. `Fp6<false>`

This is the true delayed-final-repair benchmark-only path in [Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L1436).

Properties:

- the three inner `Fp6` multiplies call `mul_fp6_pack_lufpv2_wide<false>`
- the `Fp12` linear stage uses explicit lazy helpers:
  - `lufp6_v2_add`
  - `lufp6_v2_sub`
  - `lufp6_v2_mul_by_nonresidue`
- only the final `Fp12` halves are renormed

This is the version that best matches the original optimization goal.

## Current Fp12<false> Bound Contract

The current delayed-final-repair `Fp12` path uses the canonical tower formulas:

```text
t0 = a0 * b0
t1 = a1 * b1
m  = (a0 + a1) * (b0 + b1)

r1 = m - t0 - t1
r0 = t0 + v * t1
```

Here `t0`, `t1`, and `m` are produced by `mul_fp6_pack_lufpv2_wide<false>`, so they are type-3 lazy values.

The current proved final per-half source classes used by `mul_fp12_pack_lufpv2_wide` are:

For `r1`:

- `r1.c0.re in [-11p, 10p)`
- `r1.c0.im in [-10p, 11p)`
- `r1.c1.re in [-7p, 8p)`
- `r1.c1.im in [-8p, 7p)`
- `r1.c2.re in [-6p, 6p)`
- `r1.c2.im in [-6p, 6p)`

For `r0`:

- `r0.c0.re in [-7p, 8p)`
- `r0.c0.im in [-8p, 7p)`
- `r0.c1.re in [-6p, 6p)`
- `r0.c1.im in [-6p, 6p)`
- `r0.c2.re in [-5p, 4p)`
- `r0.c2.im in [-4p, 5p)`

Those are implemented as:

- `<11,10>`
- `<10,11>`
- `<7,8>`
- `<8,7>`
- `<6,6>`
- `<5,4>`
- `<4,5>`

at the final renorm call sites in [Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L1449) and [Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L1458).

## Correctness History And Debugging Lessons

The current `Fp12<false>` path is correct, but it did not become correct on the first try.

### Invalid linear-op contract

An earlier version tried to compose lazy values with `UFp6Pack8x1::add/sub/mul_by_nonresidue`.

That was wrong.

Reason:

- those helpers are type-2 / `UFp` operators
- they do carry and conditional `p`-repair after each op
- they are not pure type-3 lazy lane-wise add/sub operators

The fix was to use explicit lazy helpers for the `Fp12` linear stage.

### `mul_by_nonresidue` alias bug

An earlier version of `lufp6_v2_mul_by_nonresidue` was used in-place without being alias-safe.

That corrupted `r.c[0]` in the `Fp12` combine and failed validation against canonical `Fp12Pack8x1::mul`.

The fix was:

- make `lufp6_v2_mul_by_nonresidue` always in-place
- make the API one-parameter only
- preserve the old `c0` component before rotating the tuple

That helper now lives at [Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L1183).

### Removal of decomposed scaled-`p` additions in live helpers

At one point the live renorm helpers built larger multiples of `p` at runtime, for example by expressing `10p` as `8p + 2p`.

That is now removed from the live helpers in [Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h).

The current state is:

- direct precomputed scales exist for `1..16`
- live scaled-`p` helpers no longer decompose large scales at runtime
- correctness still passed after this cleanup

## Implementation Roadmap

Most of the original roadmap has now been executed. The remaining work is narrower.

### 1. Keep the current three-way comparison intact

Do not collapse the comparison set. Keep benchmarking:

- `basic`
- `Fp6<true>`
- `Fp6<false>`

That separation is now essential for reasoning about where wins and losses come from.

### 2. Preserve the single-source `Fp6` kernel split

The current arrangement is the right one:

- one `mul_fp6_pack_lufpv2_wide`
- template parameter chooses whether the terminal `Fp6` renorm runs

That keeps the algebra identical between the `Fp6<true>` and `Fp6<false>` comparisons.

### 3. Explain the remaining `Fp12<false>` gap

This is now the main open problem.

The abstract work count says `Fp6<false>` should be the ceiling:

- it does the least repair work at the algorithm level
- it avoids type-2 repair after inner `Fp6` multiplies
- it has fewer total renorm ladder steps than `Fp6<true>`

But it still loses slightly in measured time.

So the next work should focus on:

- instruction count in the lazy `Fp12` glue
- register pressure and live-range length
- codegen quality of the final renorm block
- any avoidable moves or spills in the current `Fp12<false>` implementation

### 4. Do not go back to runtime-decomposed scaled `p`

The live helper cleanup that replaced runtime `8p + 2p` style construction with direct precomputed scales should stay.

That was both a performance and code-quality correction.

## Expected Comparison

The expected ordering on paper was:

1. `Fp6<false>` fastest
2. `Fp6<true>` next
3. `basic` slowest

The reason was straightforward:

- `basic` pays type-2 repair throughout
- `Fp6<true>` replaces the three inner `Fp6` multiplies with the fast live `LUFp` kernel, but still pays a terminal renorm inside each of those multiplies
- `Fp6<false>` keeps the same raw `Fp6` kernel and defers all `Fp12` repair to the very end

The measured result today is different:

- `Fp6<true>` is fastest
- `Fp6<false>` is second
- `basic` is slowest

That means the current delayed-final-repair implementation is correct and beneficial, but not yet fully optimized.

## Current Measurements

All numbers below are from the remote Linux benchmark server `18.188.231.170`.

### Fp6 wide-only comparison

This is useful background, because it shows the raw headroom available at the `Fp6` layer.

| Variant | Meaning | Time |
| --- | --- | --- |
| `basic` | core-library type-2 wide `Fp6` | `563983 ns / 256 = 2203.1 ns/Fp6` |
| `LUFp<true>` | `mul_fp6_pack_lufpv2_wide<true>` | `409402 ns / 256 = 1599.2 ns/Fp6` |
| `LUFp<false>` | `mul_fp6_pack_lufpv2_wide<false>` raw kernel only | `339919 ns / 256 = 1327.8 ns/Fp6` |

Relative to `basic`:

- `LUFp<true>` is about `1.377x` faster
- `LUFp<false>` is about `1.659x` faster

So the raw type-3 `Fp6` kernel has substantial headroom.

### Fp12 wide-only comparison

This is the main comparison to watch because it isolates the widened multiply kernel and excludes final `UFp12Pack8x1::redc`.

Latest run:

| Variant | Meaning | Time |
| --- | --- | --- |
| `basic` | core-library `Fp12` wide kernel | `793017 ns / 128 = 6195.4 ns/Fp12` |
| `Fp6<true>` | inner `Fp6` uses `mul_fp6_pack_lufpv2_wide<true>` | `701288 ns / 128 = 5478.8 ns/Fp12` |
| `Fp6<false>` | inner `Fp6` uses `mul_fp6_pack_lufpv2_wide<false>` and repairs at end of `Fp12` | `718963 ns / 128 = 5616.9 ns/Fp12` |

Relative to `basic`:

- `Fp6<true>` is about `1.130x` faster, roughly `13.0%`
- `Fp6<false>` is about `1.103x` faster, roughly `10.3%`

Relative `Fp6<true>` vs `Fp6<false>`:

- `Fp6<true>` is about `1.025x` faster, roughly `2.4%`

An earlier run immediately before this one gave the same ordering:

- `basic`: `794859 ns / 128`
- `Fp6<true>`: `704313 ns / 128`
- `Fp6<false>`: `720877 ns / 128`

So the ranking is stable.

### Fp12 full multiply comparison

This includes the final reduction back to ordinary `Fp12`.

| Variant | Meaning | Time |
| --- | --- | --- |
| `basic` | core-library full multiply | `894637 ns / 128 = 6989.4 ns/Fp12` |
| `Fp6<true>` | hybrid path | `800799 ns / 128 = 6256.2 ns/Fp12` |
| `Fp6<false>` | delayed-final-repair path | `809626 ns / 128 = 6325.2 ns/Fp12` |

Relative to `basic`:

- `Fp6<true>` is about `11.7%` faster
- `Fp6<false>` is about `10.5%` faster

So both `LUFp v2` `Fp12` paths are wins over `basic`, but the true delayed-final-repair path is not yet the best one.

## What We Have Learned

### 1. The type-3 delayed-final-repair path is correct

This is the most important status change from the original proposal.

The code in [Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L1436) now validates against canonical `Fp12Pack8x1::mul`.

So the question is no longer "can this be made correct at all?"

### 2. The type-3 path is not losing because it performs more nominal renorm ladder work

Counting renorm ladder steps:

- `Fp6<true>` inside `Fp12`:
  - 3 inner `Fp6` muls
  - 6 halves per mul
  - each half uses `4,2,1,1`
  - total seeded ladder steps: `72`
- `Fp6<false>` inside `Fp12`:
  - 12 final halves
  - 2 halves use `16,8,4,2,1,1`
  - 10 halves use `8,4,2,1,1`
  - total seeded ladder steps: `62`

So by this count, `Fp6<false>` does less nominal repair-ladder work than `Fp6<true>`.

That means "it is slower because it repairs more" is not a valid explanation.

### 3. The remaining gap is therefore implementation-level

The current evidence points to machine-level cost, not a flaw in the high-level delayed-repair idea.

The most plausible remaining sources are:

- longer live ranges in the type-3 `Fp12` glue
- higher register pressure
- worse scheduling or spilling around the final renorm block
- extra moves in the current lazy helper organization

This is not yet proved down to assembly, but it is the right category of problem now.

### 4. The raw type-3 `Fp6` speedup does not transfer linearly into `Fp12`

Standalone `Fp6<false>` is dramatically faster than standalone `basic` `Fp6`.

But `Fp12` is not just `3 * Fp6`.

`Fp12` also pays for:

- the `Fp12` linear combine
- `mul_by_nonresidue`
- final per-half renorms
- and, in the full benchmark, final `UFp12Pack8x1::redc`

So the current `Fp12` gap is about how well the type-3 `Fp12` glue is implemented, not about whether the underlying `Fp6` kernel is fundamentally good.

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
   - `basic`
   - `Fp6<true>`
   - `Fp6<false>`

This validation plan is no longer hypothetical. It is the exact benchmark structure currently present in [Fp12_tower_avx_bench.cpp](/c/Users/peter/repo/blst/benchmark/Fp12_tower_avx_bench.cpp).

## Final Guidance

If you are new to the project, the most important thing to internalize is this:

- `Fp6` was won by delaying repair to the end of the `Fp6` block
- the current concrete `Fp6` numbers are roughly `2203 ns/Fp6 -> 1599 ns/Fp6 -> 1328 ns/Fp6`
- `Fp12` has now been redesigned in benchmark-only form
- both current `LUFp v2` `Fp12` paths beat the `basic` baseline
- the true delayed-final-repair `Fp12` path is correct but still slightly slower than the `Fp6<true>` hybrid path
- the remaining work is therefore low-level performance explanation and cleanup, not first-principles design anymore

That is the accurate state of the project today.
