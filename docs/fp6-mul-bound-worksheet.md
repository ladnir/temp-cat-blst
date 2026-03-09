# `Fp6 mul` Bound Worksheet

This note is the Phase 1 worksheet for the baseline widened
`Fp6Pack8x1::mul` schedule in
[Fp6_avx.cpp](/c/Users/peter/repo/blst/blst/fields/Fp6_avx.cpp).

It derives the raw widened bound of every intermediate under a
schedule-preserving lazy-combine model.

Companion documents:

- [avx-lufp-plan.md](/c/Users/peter/repo/blst/docs/avx-lufp-plan.md)
- [lazy-renorm-strategy.md](/c/Users/peter/repo/blst/docs/lazy-renorm-strategy.md)

## Scope

This worksheet covers:

- the current baseline widened `Fp6 mul` schedule
- exact raw add/sub propagation in `LNegA_PosB_HNegC_PosD` notation
- the narrowest common wide source class for a final repair-first design
- the effect of the internal `mul_by_u_plus_1` split
- what changes conceptually if a lazy narrow type exists

This worksheet does not define the final repair kernel. It only derives the
source classes that such a kernel must accept.

## Baseline Schedule

The widened baseline implementation is:

```cpp
UFp2Pack8x1x1::mul(t0, a.c[0], b.c[0]);
UFp2Pack8x1x1::mul(t1, a.c[1], b.c[1]);
UFp2Pack8x1x1::mul(t2, a.c[2], b.c[2]);

// c0
Fp2Pack8x1x1::add(aa, a.c[1], a.c[2]);
Fp2Pack8x1x1::add(bb, b.c[1], b.c[2]);
UFp2Pack8x1x1::mul(r.c[0], aa, bb);
UFp2Pack8x1x1::sub(r.c[0], r.c[0], t1);
UFp2Pack8x1x1::sub(r.c[0], r.c[0], t2);
UFp2Pack8x1x1::mul_by_u_plus_1(t3, r.c[0]);
UFp2Pack8x1x1::add(r.c[0], t3, t0);

// c1
Fp2Pack8x1x1::add(aa, a.c[0], a.c[1]);
Fp2Pack8x1x1::add(bb, b.c[0], b.c[1]);
UFp2Pack8x1x1::mul(r.c[1], aa, bb);
UFp2Pack8x1x1::sub(r.c[1], r.c[1], t0);
UFp2Pack8x1x1::sub(r.c[1], r.c[1], t1);
UFp2Pack8x1x1::mul_by_u_plus_1(t3, t2);
UFp2Pack8x1x1::add(r.c[1], r.c[1], t3);

// c2
Fp2Pack8x1x1::add(aa, a.c[0], a.c[2]);
Fp2Pack8x1x1::add(bb, b.c[0], b.c[2]);
UFp2Pack8x1x1::mul(r.c[2], aa, bb);
UFp2Pack8x1x1::sub(r.c[2], r.c[2], t0);
UFp2Pack8x1x1::sub(r.c[2], r.c[2], t2);
UFp2Pack8x1x1::add(r.c[2], r.c[2], t1);
```

The worksheet below preserves that schedule exactly.

## Notation

For one widened `Fp` component:

- `X = lo + R * hi`
- `R = 2^384`

Track bounds as:

- `LNegA_PosB`: `-A*R <= lo < B*R`
- `HNegC_PosD`: `-C*p <= hi < D*p`

Combined form:

- `LNegA_PosB_HNegC_PosD`

## Starting Classes

Reduced narrow input `Fp` embedded in widened form:

- `N = LNeg0_Pos1_HNeg0_Pos0`

Strict widened `UFp` multiply output:

- `W = LNeg0_Pos1_HNeg0_Pos1`

For this worksheet:

- every component of `t0`, `t1`, and `t2` starts in `W`
- the narrow pre-adds `aa` and `bb` remain narrow and repaired in the current
  baseline schedule

## Raw Propagation Rules

These are the raw lazy-combine rules used in this worksheet. They do not fold
low carry or borrow into the high half yet.

If:

- `A = LNegA_PosB_HNegC_PosD`
- `B = LNegE_PosF_HNegG_PosH`

then raw add gives:

- `A + B = LNeg(A+E)_Pos(B+F)_HNeg(C+G)_Pos(D+H)`

and raw sub gives:

- `A - B = LNeg(A+F)_Pos(B+E)_HNeg(C+H)_Pos(D+G)`

These rules are intentionally more precise than the old "immediately charge a
whole extra `p`" model because low carry remains deferred.

## `mul_by_u_plus_1`

For an `Fp2` value `(x0, x1)`:

- `y0 = x0 - x1`
- `y1 = x0 + x1`

This is a split point. A single input bound feeds two downstream outputs.

That split matters because:

- it can amplify drift if the input is already wide
- it is one of the few places where an earlier renorm could be justified later
  by fanout

For `Fp6 mul`, this split occurs only once on a truly lazy value: the `c0` path.

## Derived Bounds: `c0`

Let each component of:

- `t0`, `t1`, `t2`
- the fresh product `(a1+a2)*(b1+b2)`

start in:

- `W = LNeg0_Pos1_HNeg0_Pos1`

Then for each `Fp` component:

1. `m0 - t1`
   - `W - W`
   - `LNeg1_Pos1_HNeg1_Pos1`

2. `(m0 - t1) - t2`
   - `LNeg1_Pos1_HNeg1_Pos1 - W`
   - `LNeg2_Pos1_HNeg2_Pos1`

Call that shared pre-split class:

- `C0_pre_split = LNeg2_Pos1_HNeg2_Pos1`

Now apply `mul_by_u_plus_1`:

3. low component: `x0 - x1`
   - `C0_pre_split - C0_pre_split`
   - `LNeg3_Pos3_HNeg3_Pos3`

4. high component: `x0 + x1`
   - `C0_pre_split + C0_pre_split`
   - `LNeg4_Pos2_HNeg4_Pos2`

Now add `t0`, whose components are in `W`:

5. final `c0[0]`
   - `LNeg3_Pos3_HNeg3_Pos3 + W`
   - `LNeg3_Pos4_HNeg3_Pos4`

6. final `c0[1]`
   - `LNeg4_Pos2_HNeg4_Pos2 + W`
   - `LNeg4_Pos3_HNeg4_Pos3`

So the final component classes are:

- `c0[0] = LNeg3_Pos4_HNeg3_Pos4`
- `c0[1] = LNeg4_Pos3_HNeg4_Pos3`

## Derived Bounds: `c1`

The main lazy chain is:

- fresh product
- subtract `t0`
- subtract `t1`
- add `(u+1)*t2`

The first two subtracts match the `c0` pre-split chain:

1. `m1 - t0`
   - `W - W`
   - `LNeg1_Pos1_HNeg1_Pos1`

2. `(m1 - t0) - t1`
   - `LNeg1_Pos1_HNeg1_Pos1 - W`
   - `LNeg2_Pos1_HNeg2_Pos1`

Call that:

- `C1_main = LNeg2_Pos1_HNeg2_Pos1`

Now split `t2` through `mul_by_u_plus_1`, where each component of `t2` is still
strict widened `W`:

3. low component: `t2[0] - t2[1]`
   - `W - W`
   - `LNeg1_Pos1_HNeg1_Pos1`

4. high component: `t2[0] + t2[1]`
   - `W + W`
   - `LNeg0_Pos2_HNeg0_Pos2`

Now add those to `C1_main`:

5. final `c1[0]`
   - `LNeg2_Pos1_HNeg2_Pos1 + LNeg1_Pos1_HNeg1_Pos1`
   - `LNeg3_Pos2_HNeg3_Pos2`

6. final `c1[1]`
   - `LNeg2_Pos1_HNeg2_Pos1 + LNeg0_Pos2_HNeg0_Pos2`
   - `LNeg2_Pos3_HNeg2_Pos3`

So the final component classes are:

- `c1[0] = LNeg3_Pos2_HNeg3_Pos2`
- `c1[1] = LNeg2_Pos3_HNeg2_Pos3`

## Derived Bounds: `c2`

The `c2` path is linear and has no internal split after the fresh product:

1. `m2 - t0`
   - `W - W`
   - `LNeg1_Pos1_HNeg1_Pos1`

2. `(m2 - t0) - t2`
   - `LNeg1_Pos1_HNeg1_Pos1 - W`
   - `LNeg2_Pos1_HNeg2_Pos1`

3. final add with `t1`
   - `LNeg2_Pos1_HNeg2_Pos1 + W`
   - `LNeg2_Pos2_HNeg2_Pos2`

So:

- `c2[0] = c2[1] = LNeg2_Pos2_HNeg2_Pos2`

## Common Final Wide Source Class

If we want one common final repair-first class for all `Fp6 mul` output
components, the least common upper envelope from the derived classes above is:

- `LNeg4_Pos4_HNeg4_Pos4`

That is the narrowest single common class proved by this worksheet for the raw
schedule-preserving lazy-combine version of `Fp6 mul`.

This is the main output of Phase 1.

## Important Observation: `c0` Is The Driver

`c0` determines the common class.

Specifically:

- the second component of `c0` reaches `LNeg4_Pos3_HNeg4_Pos3`
- that comes from the positive branch of `mul_by_u_plus_1`
- this is the only place in `Fp6 mul` where the pre-split lazy value fans out
  into two outputs before the final add

So if later benchmarking shows that one final repair from
`LNeg4_Pos4_HNeg4_Pos4` is too expensive, the first place to inspect is:

- the `c0` pre-split class
- and whether a selective repair before `mul_by_u_plus_1` is cheaper than
  repairing the larger post-split class

## What This Implies For Phase 2

The first repair kernel candidate for `Fp6 mul` should target:

- source class: `LNeg4_Pos4_HNeg4_Pos4`

and repair into:

- strict widened `UFp`

or, if we introduce an intermediate target later:

- a tighter lazy class that is immediately followed by a cheap strict boundary

But the important point is that Phase 2 should start from the actual
schedule-derived class above, not from placeholder names such as `Loose`.

## Relationship To Low Carry Folding

This worksheet intentionally keeps low-half slack explicit and does not charge
it into the high half at every step.

If a final repair kernel first canonicalizes low from `LNegA_PosB` back into
`LNeg0_Pos1`, that low carry changes the high half by only a small integer, not
by a whole extra multiple of `p`.

So for this worksheet:

- the `HNegA_PosB` bounds above are the right coarse source classes
- if Phase 2 needs more precision, it can refine the source class with an extra
  small-integer high slack term rather than inflating to a looser whole-`p`
  class

## Hypothetical Lazy Narrow Type

This phase is about the widened schedule, but it is still useful to ask where a
hypothetical lazy narrow type would matter.

### Where It Would Appear

In the current `Fp6 mul` schedule, the narrow-only operations are the six `Fp2`
pre-adds:

- `a1 + a2`
- `b1 + b2`
- `a0 + a1`
- `b0 + b1`
- `a0 + a2`
- `b0 + b2`

Each of those feeds exactly one multiply and has no fanout.

### What It Would Not Change

A lazy narrow type would not change:

- the wide classes of `t0`, `t1`, `t2`
- the wide combine topology after the three multiplies
- the fact that `c0` is still the bound driver on the wide side

### When It Could Help

There are two different lazy-narrow scenarios:

1. lazy narrow exists, but multiply still requires strict narrow
   - then the six pre-adds still need an immediate repair before multiply
   - because each sum is single-use and has no fanout, the topology does not
     strongly favor this as an early optimization target

2. lazy narrow exists and multiply accepts a bounded lazy narrow input
   - then the six pre-add repairs might also be removable
   - that would be a real additional optimization opportunity
   - but it changes the multiply input contract and is therefore a separate
     project from the current wide-only Phase 1/2/3 plan

### Conclusion For Phase 1

The presence of a hypothetical lazy narrow type does not change the Phase 1
wide bound worksheet.

It does change the broader optimization space:

- if lazy narrow becomes multiply-legal later, `Fp6 mul` has an additional
  pre-multiply optimization frontier
- but that frontier is separate from the current problem of replacing strict
  widened repair in the post-multiply combine region

## Summary

For the current baseline widened `Fp6 mul` schedule:

- the exact raw lazy-combine classes are derived above
- `c0` is the bound-driving path
- the narrowest common final wide source class is:
  - `LNeg4_Pos4_HNeg4_Pos4`

So the next step should be:

- benchmark and choose a real repair kernel for
  `LNeg4_Pos4_HNeg4_Pos4 -> strict UFp`
