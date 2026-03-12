# `Fp6 mul` Fused Strict-Negative Proof Note

This note documents the arithmetic and range contract of the current benchmark
`Fp6` winner:

- `mul_fp6_pack_lufpv2_wide(...)`

implemented in
[Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L1209).

The file name is historical. This note no longer describes the older delayed
block-final winner. That older centered path is still kept in the benchmark
file as:

- `mul_fp6_pack_lufpv2_wide_broken(...)`

for A/B comparison only.

What is proved here:

1. the raw unrepaired `c0`, `c1`, `c2` formulas are the correct `Fp6`
   Karatsuba formulas
2. each final `Fp` half enters the documented one-sided source envelope
3. the implemented fused strict-negative repair maps those envelopes into the
   canonical widened target class `[0, p)` per `Fp` half
4. every repair step preserves the represented value modulo `p`

This note does **not** re-prove:

- the correctness of the widened multiply/add/sub primitives
- the correctness of the Montgomery reduction path
- machine-level absence of signed overflow in the SIMD lane arithmetic

So the final theorem is:

- assuming the existing widened arithmetic primitives are exact on the stated
  integer values, and assuming `UFp6Pack8x1::redc` is correct, the current
  fused `Wide` kernel computes the correct `Fp6` product

Related notes:

- [fp6-mul-sparse5b-bounds.md](/c/Users/peter/repo/blst/docs/fp6-mul-sparse5b-bounds.md)
- [fp6-mul-debug-explainer.md](/c/Users/peter/repo/blst/docs/fp6-mul-debug-explainer.md)

## 1. Model

Fix one SIMD lane and one `Fp` half. Work entirely in ordinary integers.

Write a widened value as:

```text
X = lo + R * hi
```

where:

- `R = 2^384`
- `lo` is the low radix-48 half
- `hi` is the high radix-48 half

Use the strict widened source class:

```text
W := { lo + R*hi : 0 <= lo < R and 0 <= hi < p }.
```

Also use:

```text
C(P,N) := exact sum of P positive W-terms and N negative W-terms.
```

So:

- `W = C(1,0)`
- `W - W = C(1,1)`
- `W + W = C(2,0)`

## 1.5 Methodology

The current methodology is:

1. derive an exact one-sided source interval for each final `Fp` half
2. bias by `B*p` so the source class `[-A*p, B*p)` becomes a strict-negative
   class `[-(A+B)*p, 0)`
3. walk that strict-negative class back toward `0` with descending repair
   scales
4. fold the final `+p` canonicalization into the same step machine

This is the core idea behind the current winner.

It is not a symmetric "repair around zero" design. It is a one-sided design:

- first force the state into a purely negative class
- then only keep repaired candidates that are still negative
- stop when the remaining class is `[-p, 0)`
- then add the final `+p` to land in `[0, p)`

For a source interval:

```text
X in [-A*p, B*p)
```

the algorithm starts with:

```text
X <- X - B*p
```

so:

```text
X in [-(A+B)*p, 0).
```

Then it applies descending repair scales:

```text
4, 2, 1
```

for the current `Fp6` classes. At each repair scale `a`, the reference
operation is:

```text
X' = X + a*p
if X' < 0:
    X = X'
```

So the step does **not** mean:

- "if `X < 0`, add `a*p`"

It means:

- "if `X < -a*p`, replace `X` by `X + a*p`"

That difference is the entire reason the earlier fast substitution attempts
failed.

### Why bias to strict negative?

Because it turns an asymmetric interval problem into a one-sided shrink
problem.

Once the state is in `[-M*p, 0)`, every repair step has a simple meaning:

- try adding `a*p`
- keep the candidate only if it is still negative

That means each step shrinks the remaining interval toward `0` without ever
crossing into the wrong positive representative class.

For the current `Fp6` source classes:

- `c0.re in [-3p, 4p)` becomes `[-7p, 0)` after bias by `4p`
- `c0.im in [-4p, 3p)` becomes `[-7p, 0)` after bias by `3p`
- `c1.re in [-3p, 2p)` becomes `[-5p, 0)` after bias by `2p`
- `c1.im in [-2p, 3p)` becomes `[-5p, 0)` after bias by `3p`
- `c2.re`, `c2.im in [-2p, 2p)` become `[-4p, 0)` after bias by `2p`

Those are exactly the negative classes that the `4 -> 2 -> 1` ladder reduces
to `[-p, 0)`.

### Why does the fused version work?

The fused implementation is an optimization of the same algorithm, not a
different algorithm.

The reference step machine is:

```text
X' = X + a_i*p
if X' < 0:
    keep X'
else:
    keep X
```

The fused machine carries the previous step decision forward as a sign bit and
updates the next candidate by:

```text
X <- X + a_i*p - (!sign) * a_{i-1} * p
```

This works because:

- if the previous candidate was kept, there is nothing to undo
- if the previous candidate was rejected, undoing the prior `a_{i-1}*p` is
  exactly what returns to the accepted pre-step state

So the fused recurrence is just the reference keep-negative machine with the
"undo previous rejected step" logic folded into the next update.

The terminal `<1,1>` step is the same idea again:

- if the final candidate is still negative, add the final `+p`
- otherwise leave it alone

That is why `<1,1>` is not a special ad hoc case. It is the final
canonicalization folded into the same recurrence.

### What failed before

Several earlier ideas looked plausible but were wrong.

1. Symmetric block-level schedules were too coarse.

The old broken winner treated the final blocks as:

- `c0 in [-4p, 4p)`
- `c1 in [-3p, 3p)`
- `c2 in [-2p, 2p)`

and applied one centered block repair per `Fp2`. That was fast, but it was not
actually a proved match for the concrete per-half source classes.

2. Replacing keep-negative steps with plain sign steps was wrong.

`v2_sign_step<a>` means:

- add `a*p` on lanes already known negative

but the strict-negative algorithm needs:

- add `a*p` only when the value is below `-a*p`

Those are not equivalent. The interval `[-a*p, 0)` is the difference.

3. The old hand-tuned `c1` schedule was insufficient.

The restored historical `c1` schedule:

- `<2;3,1>` for `c1.re`
- `<3;3,1>` for `c1.im`

looks tighter, but it does not actually shrink `[-5p, 0)` to `[-p, 0)`. That
is why the old restored `wide` failed debug validation even though `c0` and
`c2` looked fine.

4. Generic correctness without a debug oracle was not enough.

The project only stabilized once `renorm_any` was used as the explicit oracle
for:

- raw combine correctness
- candidate renorm equivalence
- fast-path postconditions

That is the practical lesson: the bounds proof and the executable postcondition
have to agree.

## 2. Raw `Fp6` Formulas

The kernel computes:

```text
t0 = a0*b0
t1 = a1*b1
t2 = a2*b2
m0 = (a1 + a2)(b1 + b2)
m1 = (a0 + a1)(b0 + b1)
m2 = (a0 + a2)(b0 + b2)
```

and then forms:

```text
c0 = (u+1)(m0 - t1 - t2) + t0
c1 = (m1 - t0 - t1) + (u+1)t2
c2 = (m2 - t0 - t2) + t1
```

This is exactly what the code does in
[Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L1219),
[Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L1235),
and
[Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L1249).

### Lemma 2.1

For `x = x_re + u*x_im`,

```text
(u+1)x = (x_re - x_im) + u(x_re + x_im).
```

### Proof

By distributivity and `u^2 = -1` in `Fp2`:

```text
(u+1)(x_re + u*x_im)
= x_re*u + x_re + u^2*x_im + u*x_im
= (x_re - x_im) + u(x_re + x_im).
```

## 3. Raw Per-`Fp` Expansion

Each component of `t0`, `t1`, `t2`, `m0`, `m1`, `m2` is one strict widened
term in `W`.

Using Lemma 2.1:

### `c0`

Let:

```text
x = m0 - t1 - t2.
```

Then:

```text
c0 = (u+1)x + t0.
```

So:

```text
c0.re = x.re - x.im + t0.re
      = m0.re - t1.re - t2.re - m0.im + t1.im + t2.im + t0.re
      = C(4,3)

c0.im = x.re + x.im + t0.im
      = m0.re - t1.re - t2.re + m0.im - t1.im - t2.im + t0.im
      = C(3,4)
```

### `c1`

```text
c1 = (m1 - t0 - t1) + (u+1)t2
```

So:

```text
c1.re = m1.re - t0.re - t1.re + t2.re - t2.im
      = C(2,3)

c1.im = m1.im - t0.im - t1.im + t2.re + t2.im
      = C(3,2)
```

### `c2`

```text
c2 = (m2 - t0 - t2) + t1
```

So:

```text
c2.re = m2.re - t0.re - t2.re + t1.re
      = C(2,2)

c2.im = m2.im - t0.im - t2.im + t1.im
      = C(2,2)
```

## 4. Envelope From `C(P,N)`

### Lemma 4.1

After a full signed carry of a `C(P,N)` value:

- the low half is canonicalized into `[0, R)`
- the carried high half lies in `[-N*p, P*p)`

### Proof

Each positive `W` term contributes:

- low in `[0, R)`
- high in `[0, p)`

Each negative `W` term contributes:

- low in `(-R, 0]`
- high in `(-p, 0]`

Summing `P` positive and `N` negative terms gives a raw value whose lower part
lies in `(-N*R, P*R)` and whose upper part lies in `(-N*p, P*p)`.

The function `v2_signed_carry_range(...)` in
[Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L325)
propagates signed radix carries from the low half into the high half without
changing the represented integer. After this carry, the low half is canonical
in `[0, R)`, and the remaining signed drift is in the high half.

So the high interval is exactly `[-N*p, P*p)`.

### Corollary 4.2

Applying Lemma 4.1 to Section 3 gives:

```text
c0.re in [-3p, 4p)
c0.im in [-4p, 3p)

c1.re in [-3p, 2p)
c1.im in [-2p, 3p)

c2.re in [-2p, 2p)
c2.im in [-2p, 2p)
```

These are exactly the live kernel instantiations:

- `c0.re -> <A, B> = <3, 4>`
- `c0.im -> <4, 3>`
- `c1.re -> <3, 2>`
- `c1.im -> <2, 3>`
- `c2.re`, `c2.im -> <2, 2>`

See
[Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L1231),
[Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L1245),
and
[Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L1255).

## 5. Reference Strict-Negative Repair

For a source interval `[-A*p, B*p)`, define:

```text
Bias_B(X) := X - B*p
```

and for `a >= 1`:

```text
Keep_a(X) :=
    X + a*p   if X + a*p < 0
    X         otherwise.
```

Let:

```text
a0 := 2^floor(log2(A + B)).
```

The reference one-sided repair ladder is:

```text
RefFold_{A,B}(X)
= (+p) after Keep_1(...Keep_2(Keep_4(Bias_B(X)))...)
```

with descending powers of two from `a0` to `1`.

For the current `Fp6` source classes:

- `A + B = 7` for `c0.re` and `c0.im`
- `A + B = 5` for `c1.re` and `c1.im`
- `A + B = 4` for `c2.re` and `c2.im`

So every current half uses the same descending ladder:

```text
4, 2, 1
```

followed by the final `+p`.

### Lemma 5.1

`Bias_B` preserves the represented value modulo `p`.

### Proof

Immediate, because it subtracts a multiple of `p`.

### Lemma 5.2

`Keep_a` preserves the represented value modulo `p`.

### Proof

It changes the value by either `0` or `a*p`.

### Lemma 5.3

For the current source classes, after bias we have:

- `[-3p, 4p) - 4p = [-7p, 0)`
- `[-4p, 3p) - 3p = [-7p, 0)`
- `[-3p, 2p) - 2p = [-5p, 0)`
- `[-2p, 3p) - 3p = [-5p, 0)`
- `[-2p, 2p) - 2p = [-4p, 0)`

Applying `Keep_4`, `Keep_2`, `Keep_1` maps each of these to `[-p, 0)`.

### Proof

For `[-7p, 0)`:

- after `Keep_4`: `[-4p, 0)`
- after `Keep_2`: `[-2p, 0)`
- after `Keep_1`: `[-p, 0)`

For `[-5p, 0)`:

- after `Keep_4`: `[-4p, 0)`
- after `Keep_2`: `[-2p, 0)`
- after `Keep_1`: `[-p, 0)`

For `[-4p, 0)`:

- after `Keep_4`: `[-4p, 0)`
- after `Keep_2`: `[-2p, 0)`
- after `Keep_1`: `[-p, 0)`

So every current half reaches `[-p, 0)`, and the final `+p` maps that to
`[0, p)`.

## 6. Fused Strict-Negative Repair

The live code does not execute the reference ladder literally. It uses the
fused helper
[v2_keep_negative_fuse_step](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L725)
through
[lufp_v2_renorm_negative_fold_fused](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L789).

The first fused step is the seed:

```text
Candidate_0 = Bias_B(X) + a0*p
```

and returns:

```text
s0 = 1 if Candidate_0 < 0
     0 otherwise.
```

For later steps the fused recurrence is:

```text
Candidate_i
= Candidate_{i-1} + a_i*p - (1 - s_{i-1}) * a_{i-1} * p
```

and it returns:

```text
s_i = 1 if Candidate_i < 0
      0 otherwise.
```

The terminal `<1,1>` step folds the final `+p` into the same recurrence:

```text
Final
= Candidate_last + p - (1 - s_last) * 1*p.
```

### Lemma 6.1

The fused recurrence produces exactly the same candidate sequence as the
reference `Keep_a` ladder.

### Proof

Assume `Candidate_{i-1}` is the candidate produced at the previous reference
step, and `s_{i-1}` records whether that candidate was still negative.

There are two cases.

If `s_{i-1} = 1`, the previous candidate was kept by the reference ladder, so
the accepted state entering the next step is exactly `Candidate_{i-1}`. The
next reference candidate is therefore:

```text
Candidate_{i-1} + a_i*p.
```

The fused recurrence gives the same value because the subtract term vanishes.

If `s_{i-1} = 0`, the previous candidate was rejected by the reference ladder.
Since the previous candidate was formed by adding `a_{i-1}*p`, the accepted
state is:

```text
Candidate_{i-1} - a_{i-1}*p.
```

So the next reference candidate is:

```text
Candidate_{i-1} - a_{i-1}*p + a_i*p,
```

which is exactly the fused recurrence when `s_{i-1} = 0`.

So the fused step produces the same candidate as the reference ladder in both
cases.

### Lemma 6.2

The terminal fused `<1,1>` step produces the same final state as the explicit
final `+p` in the reference ladder.

### Proof

Again split into two cases.

If `s_last = 1`, the last candidate was kept, so the reference ladder adds
`+p` to that candidate. The fused terminal step becomes:

```text
Candidate_last + p.
```

If `s_last = 0`, the last candidate was rejected. The accepted state is:

```text
Candidate_last - p.
```

Then the explicit final `+p` returns:

```text
Candidate_last.
```

The fused terminal step with `s_last = 0` is:

```text
Candidate_last + p - p = Candidate_last.
```

So the fused terminal step matches the reference finalization in both cases.

### Corollary 6.3

`lufp_v2_renorm_negative_fold_fused<A,B>(...)` is value-equivalent to the
reference strict-negative ladder for the same `<A, B>` source class.

## 7. Final Range Contract

By Corollary 4.2, every half of the raw `Fp6` result lies in its stated
`[-A*p, B*p)` source class.

By Lemma 5.3, the reference strict-negative ladder maps every current source
class to `[-p, 0)`, and the final `+p` maps that to `[0, p)`.

By Corollary 6.3, the fused implementation is value-equivalent to that
reference ladder.

Therefore the live fused repair maps every current `Fp6` half to:

```text
[0, p).
```

Since every repair stage changes the value only by multiples of `p`, the final
canonical widened outputs remain congruent modulo `p` to the correct raw
product.

## 8. Final Theorem

### Theorem 8.1

For each lane, `mul_fp6_pack_lufpv2_wide(...)` computes three widened `Fp2`
outputs `c0`, `c1`, `c2` such that:

1. each output is congruent modulo `p` to the correct raw `Fp6` product block
2. every final `Fp` half lies in `[0, p)`

### Proof

Section 2 and Section 3 give the correct raw algebraic formulas.

Corollary 4.2 gives the source envelopes:

- `c0.re in [-3p, 4p)`
- `c0.im in [-4p, 3p)`
- `c1.re in [-3p, 2p)`
- `c1.im in [-2p, 3p)`
- `c2.re in [-2p, 2p)`
- `c2.im in [-2p, 2p)`

Section 5 proves that the reference strict-negative ladder maps those current
source classes into `[0, p)`.

Section 6 proves that the fused implementation is equivalent to that reference
ladder.

Therefore the live fused repair produces canonical widened outputs in `[0, p)`
while preserving congruence modulo `p`.

### Corollary 8.2

Assuming `UFp6Pack8x1::redc` correctly reduces widened values from that target
class, the reduced output of `mul_fp6_pack_lufpv2_wide(...)` is the correct
`Fp6` product.

## 9. Historical Note

The older delayed-centered path:

- `mul_fp6_pack_lufpv2_wide_broken(...)`

is still useful as a performance reference, but it is no longer the current
baseline and is not the subject of this proof note.
