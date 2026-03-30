# `Fp12 mul` LUFp v2 First-Pass Contract

This note records the first benchmark-only composable `Fp6` interface for the
`Fp12` follow-up and the resulting `Fp12` final-half bounds.

It is intentionally narrow in scope:

- benchmark-only
- multiply only
- linear `Fp12` combine only
- no claim that the interface is a general-purpose lazy `Fp6` type

## 1. Post-carry `Fp6` interface

The new helper is:

```text
mul_fp6_pack_lufpv2_wide_postcarry(...)
```

Its contract is:

1. compute the same raw `Fp6` Karatsuba result as the live
   `mul_fp6_pack_lufpv2_wide(...)` winner
2. perform a full signed carry on each final `Fp` half
3. stop before the terminal fused strict-negative repair

So the returned state is:

- exact modulo `p`
- low half canonicalized by the signed carry pass
- high half still allowed to drift in the proved one-sided interval

Per final `Fp` half, the post-carry source classes are exactly the current live
winner source classes:

- `c0.re in [-3p, 4p)`
- `c0.im in [-4p, 3p)`
- `c1.re in [-3p, 2p)`
- `c1.im in [-2p, 3p)`
- `c2.re in [-2p, 2p)`
- `c2.im in [-2p, 2p)`

This is the composable contract consumed by the first `Fp12` path.

## 2. Valid operations on the post-carry interface

Two different statements matter here.

### 2.1 Candidate composable contract

The first intended `Fp12` design wanted the post-carry `Fp6` result to be valid
for:

- component-wise widened add
- component-wise widened subtract
- `Fp2` multiply by `(u + 1)`
- `Fp6` multiply by the nonresidue `v`
- one final per-half strict-negative repair

That is enough for the canonical `Fp12` tower formulas:

```text
t0 = a0 * b0
t1 = a1 * b1
m  = (a0 + a1) * (b0 + b1)

r1 = m - t0 - t1
r0 = t0 + v * t1
```

### 2.2 Current promoted contract

That stronger composable contract did **not** validate end-to-end in the first
`Fp12` experiment.

So the current benchmark-only code only relies on the weaker contract:

- the post-carry `Fp6` result computes the same raw product as the live winner
- it is valid input to a slower obvious per-product repair (`renorm_any`)
- after that per-product repair, the ordinary `Fp12` combine is correct

This note therefore does **not** claim the current post-carry interface is
already valid for direct composable `Fp12` widened add/sub/nonresidue use in
promoted code.

It also does **not** claim the post-carry interface is directly valid for:

- another `Fp6` multiply
- direct `UFp6Pack8x1::redc` without a final repair
- arbitrary future lazy schedules

## 3. `Fp12` bounds from the post-carry interface

Every `t0`, `t1`, and `m` product uses the same post-carry `Fp6` bounds from
Section 1.

### 3.1 `r1 = m - t0 - t1`

Per half:

- `r1.c0.re in [-11p, 10p)`
- `r1.c0.im in [-10p, 11p)`
- `r1.c1.re in [-7p, 8p)`
- `r1.c1.im in [-8p, 7p)`
- `r1.c2.re in [-6p, 6p)`
- `r1.c2.im in [-6p, 6p)`

Derivation is direct interval propagation:

- subtracting `[-A*p, B*p)` contributes `+B` to the negative side and `+A` to
  the positive side
- so `m - t0 - t1` is just the component-wise accumulation of those terms

### 3.2 `r0 = t0 + v * t1`

For `x = t1.c2`, `mul_by_u_plus_1(x)` maps:

- `x.re, x.im in [-2p, 2p)` to
- `x.re - x.im in [-4p, 4p)`
- `x.re + x.im in [-4p, 4p)`

So `v * t1 = ((u + 1) * t1.c2, t1.c0, t1.c1)` gives:

- first `Fp2` half: `([-4p, 4p), [-4p, 4p))`
- second `Fp2` half: `t1.c0`
- third `Fp2` half: `t1.c1`

Adding `t0` then gives:

- `r0.c0.re in [-7p, 8p)`
- `r0.c0.im in [-8p, 7p)`
- `r0.c1.re in [-6p, 6p)`
- `r0.c1.im in [-6p, 6p)`
- `r0.c2.re in [-5p, 4p)`
- `r0.c2.im in [-4p, 5p)`

## 4. Candidate delayed repair schedule

The repair uses the same one-sided strict-negative framework as the live `Fp6`
winner.

For a source class `[-A*p, B*p)`, the first repair scale is:

```text
bit_floor(A + B)
```

That gives the smallest descending power-of-two ladder that can still shrink
the post-bias class `[-(A+B)*p, 0)` to `[-p, 0)`.

For the first delayed-final candidate:

- `r1.c0.re` and `r1.c0.im` need `16 -> 8 -> 4 -> 2 -> 1 -> <1,1>`
- every other final half needs `8 -> 4 -> 2 -> 1 -> <1,1>`

Reason:

- `r1.c0.*` has width `A + B = 21`
- the other final halves have widths `15`, `12`, or `9`
- `8 + 4 + 2 + 1 = 15`, so width `21` requires a `16` step
- width `9..15` requires an `8` step but no `16` step

## 5. Current implementation status

The current correct benchmark-only code computes:

1. `t0`, `t1`, and `m` with `mul_fp6_pack_lufpv2_wide_postcarry(...)`
2. repair those three `Fp6` products with `renorm_any`
3. run the ordinary widened `Fp12` linear combine from those repaired products

The first attempted fast terminal schedule was the per-half strict-negative
repair from Section 4.

That schedule is still the right candidate to investigate next, but it is **not
promoted in code yet** because end-to-end validation against `renorm_any`
showed a remaining proof gap on several final halves.

So the current benchmark-only `mul_fp12_pack_lufpv2_wide(...)` uses the slower
obvious product-level repair:

```text
renorm_any on t0, t1, and m
```

This gives a correct first `Fp12` path built on the new post-carry `Fp6`
interface while preserving the unresolved composable-lazy and fast-final-repair
work as a clearly identified next step.
