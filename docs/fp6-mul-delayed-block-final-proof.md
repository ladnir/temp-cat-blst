# `WideDelayedBlockFinal` Proof Note

This note gives an explicit elementary proof for the arithmetic and range
contract of:

- `mul_fp6_pack_lufpv2_wide_delayed_block_final(...)`

implemented in
[Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L1049).

It proves:

1. the raw unrepaired `c0`, `c1`, `c2` formulas are the correct `Fp6` Karatsuba
   formulas
2. each final `Fp` half enters the documented delayed-repair source envelope
3. the implemented delayed repair maps those envelopes into explicit nonnegative
   target envelopes
4. every repair step preserves the represented value modulo `p`

This note does **not** re-prove:

- the correctness of the library Montgomery reduction itself
- machine-level absence of signed overflow in the SIMD lane arithmetic

So the final theorem is:

- assuming the existing widened arithmetic primitives are exact on the stated
  integer values, and assuming `UFpPack8x1::redc` is correct on the resulting
  widened outputs, `WideDelayedBlockFinal` computes the correct `Fp6`
  product

Related notes:

- [fp6-mul-sparse5b-bounds.md](/c/Users/peter/repo/blst/docs/fp6-mul-sparse5b-bounds.md)
- [fp6-mul-repair-path-plan.md](/c/Users/peter/repo/blst/docs/fp6-mul-repair-path-plan.md)

## 1. Model

Fix one SIMD lane and one `Fp` half. We work entirely in ordinary integers.

Write a widened value as:

```text
X = lo + R * hi
```

where:

- `R = 2^384`
- `lo` is the lower 384-bit radix value
- `hi` is the upper 384-bit radix value

We use the standard strict widened source class:

```text
W := { lo + R*hi : 0 <= lo < R and 0 <= hi < p }.
```

We also use:

```text
C(P,N) := exact sum of P positive W-terms and N negative W-terms.
```

So:

- `W = C(1,0)`
- `W - W = C(1,1)`
- `W + W = C(2,0)`

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
[Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L1055),
[Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L1059),
[Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L1074),
and
[Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L1088).

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
[Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L311)
propagates signed radix carries from the low half into the high half without
changing the represented integer. After this carry, the low half is canonical
in `[0, R)`, and the entire remaining signed drift is in the high half.

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

Therefore:

```text
c0 fits symmetric envelope [-4p, 4p)
c1 fits symmetric envelope [-3p, 3p)
c2 fits symmetric envelope [-2p, 2p)
```

These are exactly the source envelopes used by the kernel comments at
[Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L1071),
[Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L1085),
and
[Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L1099).

## 5. Delayed Repair Operator

The helper
[lufp_v2_renorm_delayed_centered](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L539)
does three things:

```text
1. subtract BiasScale * p from the high half
2. perform full signed carry
3. apply a sequence of sign steps
```

Formally, define:

```text
Bias_k(X) := X - k*p
```

and for `s >= 1`:

```text
Step_s(X) :=
    X + s*p   if X < 0
    X         if X >= 0
```

Then the integer effect of the helper is:

```text
Repair_{k; s1,...,sr}(X)
= Step_sr(...Step_s2(Step_s1(Bias_k(X)))).
```

### Lemma 5.1

`Bias_k` changes the represented value by `-k*p`, hence preserves the value
modulo `p`.

### Proof

Immediate, because subtracting `k*p` changes the integer by a multiple of `p`.

### Lemma 5.2

Each `Step_s` changes the represented value by either `0` or `s*p`, hence
preserves the value modulo `p`.

### Proof

By definition:

- if `X >= 0`, `Step_s(X) = X`
- if `X < 0`, `Step_s(X) = X + s*p`

In both cases the difference is a multiple of `p`.

### Lemma 5.3

If `X` lies in `[A, B)`, then:

```text
Step_s(X) lies in [min(A + s*p, 0), max(B, s*p)).
```

### Proof

Split into two disjoint cases.

If `X < 0`, then:

```text
X in [A, 0)
=> X + s*p in [A + s*p, s*p).
```

If `X >= 0`, then:

```text
X in [0, B).
```

Taking the union gives the stated interval.

## 6. Explicit Block Repairs

The kernel uses:

- `c0`: `Repair_{3; 4,2,1}`
- `c1`: `Repair_{2; 4,1}`
- `c2`: `Repair_{1; 2,1}`

See:

- [Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L1072)
- [Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L1086)
- [Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L1100)

We now compute each target interval explicitly.

### 6.1 `c0`

From Corollary 4.2:

```text
c0 in [-4p, 4p).
```

After bias by `3p`:

```text
[-4p, 4p) - 3p = [-7p, p).
```

Apply `Step_4`. By Lemma 5.3:

```text
[-7p, p) -> [-3p, 4p).
```

Apply `Step_2`:

```text
[-3p, 4p) -> [-p, 4p).
```

Apply `Step_1`:

```text
[-p, 4p) -> [0, 4p).
```

So:

```text
Repair_{3;4,2,1}(c0) in [0, 4p).
```

### 6.2 `c1`

From Corollary 4.2:

```text
c1 in [-3p, 3p).
```

After bias by `2p`:

```text
[-3p, 3p) - 2p = [-5p, p).
```

Apply `Step_4`:

```text
[-5p, p) -> [-p, 4p).
```

Apply `Step_1`:

```text
[-p, 4p) -> [0, 4p).
```

So:

```text
Repair_{2;4,1}(c1) in [0, 4p).
```

### 6.3 `c2`

From Corollary 4.2:

```text
c2 in [-2p, 2p).
```

After bias by `1p`:

```text
[-2p, 2p) - p = [-3p, p).
```

Apply `Step_2`:

```text
[-3p, p) -> [-p, 2p).
```

Apply `Step_1`:

```text
[-p, 2p) -> [0, 2p).
```

So:

```text
Repair_{1;2,1}(c2) in [0, 2p).
```

## 7. Congruence Preservation

### Lemma 7.1

Each block repair preserves the raw block value modulo `p`.

### Proof

Each block repair is a composition of:

- one bias by a multiple of `p`
- several conditional additions by multiples of `p`

By Lemma 5.1 and Lemma 5.2, every stage preserves the represented value modulo
`p`. Therefore the full composition preserves the represented value modulo `p`.

## 8. Final Theorem

### Theorem 8.1

For each lane, `mul_fp6_pack_lufpv2_wide_delayed_block_final(...)` computes
three widened `Fp2` outputs `c0`, `c1`, `c2` such that:

1. each output is congruent modulo `p` to the correct raw `Fp6` product block
2. `c0` ends in `[0, 4p)` per `Fp` half
3. `c1` ends in `[0, 4p)` per `Fp` half
4. `c2` ends in `[0, 2p)` per `Fp` half

### Proof

Section 2 and Section 3 give the correct raw algebraic formulas.

Corollary 4.2 gives the source intervals:

- `c0 in [-4p, 4p)`
- `c1 in [-3p, 3p)`
- `c2 in [-2p, 2p)`

Section 6 proves that the implemented repairs map these intervals to:

- `c0 -> [0, 4p)`
- `c1 -> [0, 4p)`
- `c2 -> [0, 2p)`

Section 7 proves that each repair preserves congruence modulo `p`.

Therefore the repaired block outputs are congruent to the correct raw `Fp6`
product and lie in the claimed explicit target intervals.

### Corollary 8.2

Assuming the existing reduction path

- `UFp6Pack8x1::redc`

correctly reduces widened values from those target intervals, the reduced output
of `WideDelayedBlockFinal` is the correct `Fp6` product.

## 9. What Is Still Assumed

This note is a proof of the kernel-level arithmetic and range contract. It
still relies on the following already-existing components:

1. the widened multiply/add/sub primitives compute the exact intended integers
2. the carry primitives preserve the represented integer exactly
3. `UFpPack8x1::redc` correctly reduces the resulting widened outputs

Those are shared trust assumptions for the rest of the widened arithmetic code
too. They are not unique to `WideDelayedBlockFinal`.
