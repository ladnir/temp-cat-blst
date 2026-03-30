# `Fp6 mul` Sparse5B Bounds Worksheet

This note analyzes the former sparse winner:

- `mul_fp6_pack_lufpv2_wide_step_candidate_sparse<0x3C2u>(...)`
- implemented in
  [Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h)

The goal is to do the bound work that the earlier debug note said was still
missing:

- compute symbolic ranges for each intermediate
- propagate those through the sparse add/sub schedule
- identify the exact checkpoints that still require a proof obligation

This document is intentionally conservative.

It gives:

- the range propagation you can justify from the current code shape
- the checkpoints that are symbolically safe under the current model
- the checkpoints that are still only empirically justified
- the end-only per-`Fp` bounds after logically inlining `mul_by_u_plus_1`

It does **not** claim a full proof for all valid inputs. In fact, `Sparse5B`
now has a concrete valid-input counterexample, so the final sections explain
where the original proof attempt stopped and what the first-principles end-only
model says instead.

Companion notes:

- [fp6-mul-debug-explainer.md](/c/Users/peter/repo/blst/docs/fp6-mul-debug-explainer.md)
- [fp6-mul-bound-worksheet.md](/c/Users/peter/repo/blst/docs/fp6-mul-bound-worksheet.md)

## Scope

This worksheet covers only `Sparse5B`:

- `Sparse5B`
- mask `0x3C2`

Its active repairs are:

- `c0_sub_t2`
- `c1_uplus1_t2`
- `c1_add_t3`
- `c2_sub_t0`
- `c2_sub_t2`

Its skipped repairs are:

- `c0_sub_t1`
- `c0_uplus1`
- `c0_add_t0`
- `c1_sub_t0`
- `c1_sub_t1`
- `c2_add_t1`

The reordering of the implementation into `c0`, `c1`, `c2` blocks does not
change any bounds. It only changes live ranges and code shape.

This note also records the later result that `Sparse5B` is not correct for all
valid inputs: the targeted adversarial search found a concrete counterexample on
the `c0` path.

## Notation

Reuse the notation from
[fp6-mul-bound-worksheet.md](/c/Users/peter/repo/blst/docs/fp6-mul-bound-worksheet.md):

- `X = lo + R * hi`
- `R = 2^384`

Track bounds as:

- `LNegA_PosB`: `-A*R <= lo < B*R`
- `HNegC_PosD`: `-C*p <= hi < D*p`
- combined: `LNegA_PosB_HNegC_PosD`

Use the same strict widened product class:

- `W = LNeg0_Pos1_HNeg0_Pos1`

So:

- each component of `t0`, `t1`, `t2` starts in `W`
- each fresh Karatsuba side product also starts in `W`

## Raw Rules

For:

- `A = LNegA_PosB_HNegC_PosD`
- `B = LNegE_PosF_HNegG_PosH`

raw add gives:

- `A + B = LNeg(A+E)_Pos(B+F)_HNeg(C+G)_Pos(D+H)`

raw sub gives:

- `A - B = LNeg(A+F)_Pos(B+E)_HNeg(C+H)_Pos(D+G)`

For `mul_by_u_plus_1(x0, x1)`:

- low output is `x0 - x1`
- high output is `x0 + x1`

## What `ExactAdd` Needs To Prove

The sparse kernel uses `lufp2_v2_renorm_exact_add(...)` as its only repair.

At the code level, that primitive is:

- low carry through the whole widened value
- sign test on the top limb
- one masked add of `p` if the carried value is negative

See:

- [Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L327)
- [Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L349)

That means one `ExactAdd` repair is only obviously sufficient when the carried
high half is already known to lie in a one-step neighborhood of `[0, p)`.

The current symbolic model can therefore justify `ExactAdd` cleanly only for
states that are still "one deferred add/sub away" from `W`.

Once a checkpoint drifts to a two-step class such as:

- `LNeg2_Pos1_HNeg2_Pos1`
- `LNeg2_Pos2_HNeg2_Pos2`

the old independent low/high interval model no longer proves that one sign
repair is enough. At that point, the proof would need either:

- tighter low/high correlation, or
- a direct proof of the `ExactAdd` implementation on that wider source class

That distinction is the main result of this worksheet.

## Sparse5B Schedule

The active sparse path is:

```cpp
// c0
m0 = (a1 + a2) * (b1 + b2);
c0 = m0 - t1;
c0 = c0 - t2;
repair(c0);              // active
t3 = mul_by_u_plus_1(c0);
c0 = t3 + t0;

// c1
m1 = (a0 + a1) * (b0 + b1);
c1 = m1 - t0;
c1 = c1 - t1;
t3 = mul_by_u_plus_1(t2);
repair(t3);              // active
c1 = c1 + t3;
repair(c1);              // active

// c2
m2 = (a0 + a2) * (b0 + b2);
c2 = m2 - t0;
repair(c2);              // active
c2 = c2 - t2;
repair(c2);              // active
c2 = c2 + t1;
```

## Symbolic Propagation

### `c0`

Start:

- `m0 = W`
- `t1 = W`
- `t2 = W`
- `t0 = W`

Step 1:

- `c0_1 = m0 - t1`
- `W - W`
- `c0_1 = LNeg1_Pos1_HNeg1_Pos1`

Step 2:

- `c0_2 = c0_1 - t2`
- `LNeg1_Pos1_HNeg1_Pos1 - W`
- `c0_2 = LNeg2_Pos1_HNeg2_Pos1`

Step 3:

- active repair at `c0_sub_t2`
- if that repair succeeds, it returns the checkpoint to strict widened class:
- `c0_2_repaired = W`

Step 4: split through `mul_by_u_plus_1`

From `W`:

- low branch: `W - W = LNeg1_Pos1_HNeg1_Pos1`
- high branch: `W + W = LNeg0_Pos2_HNeg0_Pos2`

So:

- `t3[0] = LNeg1_Pos1_HNeg1_Pos1`
- `t3[1] = LNeg0_Pos2_HNeg0_Pos2`

No repair is applied here in `Sparse5B`.

Step 5: final add with `t0`

Final `c0.re`:

- `t3[0] + t0`
- `LNeg1_Pos1_HNeg1_Pos1 + W`
- `c0.re = LNeg1_Pos2_HNeg1_Pos2`

Final `c0.im`:

- `t3[1] + t0`
- `LNeg0_Pos2_HNeg0_Pos2 + W`
- `c0.im = LNeg0_Pos3_HNeg0_Pos3`

So the `c0` outputs leave the kernel without a final repair in classes:

- `c0.re = LNeg1_Pos2_HNeg1_Pos2`
- `c0.im = LNeg0_Pos3_HNeg0_Pos3`

These are both strictly narrower than the raw baseline worksheet's final common
envelope `LNeg4_Pos4_HNeg4_Pos4`.

### `c1`

Start:

- `m1 = W`
- `t0 = W`
- `t1 = W`

Step 1:

- `c1_1 = m1 - t0`
- `W - W`
- `c1_1 = LNeg1_Pos1_HNeg1_Pos1`

Step 2:

- `c1_2 = c1_1 - t1`
- `LNeg1_Pos1_HNeg1_Pos1 - W`
- `c1_2 = LNeg2_Pos1_HNeg2_Pos1`

No repair is applied yet.

Step 3: split `t2` through `mul_by_u_plus_1`

Since `t2 = W`:

- `t3_raw.re = W - W = LNeg1_Pos1_HNeg1_Pos1`
- `t3_raw.im = W + W = LNeg0_Pos2_HNeg0_Pos2`

Step 4:

- active repair at `c1_uplus1_t2`
- if that repair succeeds on both halves:
- `t3 = W`

Step 5:

- `c1_3 = c1_2 + t3`
- `LNeg2_Pos1_HNeg2_Pos1 + W`
- `c1_3 = LNeg2_Pos2_HNeg2_Pos2`

Step 6:

- active repair at `c1_add_t3`
- if that repair succeeds:
- final `c1 = W`

So the whole `c1` path hinges on two proof obligations:

- `ExactAdd` must succeed on the split `t2` class
- `ExactAdd` must also succeed on `LNeg2_Pos2_HNeg2_Pos2`

### `c2`

Start:

- `m2 = W`
- `t0 = W`
- `t2 = W`
- `t1 = W`

Step 1:

- `c2_1 = m2 - t0`
- `W - W`
- `c2_1 = LNeg1_Pos1_HNeg1_Pos1`

Step 2:

- active repair at `c2_sub_t0`
- if that repair succeeds:
- `c2_1_repaired = W`

Step 3:

- `c2_2 = c2_1_repaired - t2`
- `W - W`
- `c2_2 = LNeg1_Pos1_HNeg1_Pos1`

Step 4:

- active repair at `c2_sub_t2`
- if that repair succeeds:
- `c2_2_repaired = W`

Step 5:

- final `c2 = c2_2_repaired + t1`
- `W + W`
- `c2 = LNeg0_Pos2_HNeg0_Pos2`

So the final unrepaired `c2` output class is:

- `c2.re = c2.im = LNeg0_Pos2_HNeg0_Pos2`

Again this is much narrower than the raw baseline `LNeg4_Pos4_HNeg4_Pos4`
envelope.

## Checkpoint Table

| Checkpoint | Symbolic class | Repair applied? | Current model verdict |
|---|---|---:|---|
| `c0_sub_t1` | `LNeg1_Pos1_HNeg1_Pos1` | No | safe to defer one more step |
| `c0_sub_t2` | `LNeg2_Pos1_HNeg2_Pos1` | Yes | not proved by current coarse model |
| `c0_uplus1` low | `LNeg1_Pos1_HNeg1_Pos1` | No | safe as intermediate |
| `c0_uplus1` high | `LNeg0_Pos2_HNeg0_Pos2` | No | safe as intermediate |
| final `c0.re` | `LNeg1_Pos2_HNeg1_Pos2` | No | narrower than raw worksheet envelope |
| final `c0.im` | `LNeg0_Pos3_HNeg0_Pos3` | No | narrower than raw worksheet envelope |
| `c1_sub_t0` | `LNeg1_Pos1_HNeg1_Pos1` | No | safe to defer one more step |
| `c1_sub_t1` | `LNeg2_Pos1_HNeg2_Pos1` | No | becomes proof obligation at next add |
| `c1_uplus1_t2.re` | `LNeg1_Pos1_HNeg1_Pos1` | Yes | plausibly in one-step domain |
| `c1_uplus1_t2.im` | `LNeg0_Pos2_HNeg0_Pos2` | Yes | not proved by current coarse model |
| `c1_add_t3` | `LNeg2_Pos2_HNeg2_Pos2` | Yes | not proved by current coarse model |
| `c2_sub_t0` | `LNeg1_Pos1_HNeg1_Pos1` | Yes | symbolically consistent with one-step repair |
| `c2_sub_t2` | `LNeg1_Pos1_HNeg1_Pos1` | Yes | symbolically consistent with one-step repair |
| final `c2` | `LNeg0_Pos2_HNeg0_Pos2` | No | narrower than raw worksheet envelope |

## What The Current Bounds Actually Prove

The sparse schedule is clearly *structurally better* than the raw no-repair
baseline:

- `c0` and `c2` final classes are much narrower than the raw worksheet classes
- the `c2` path is fully aligned with a one-step repair story

More concretely, the current symbolic worksheet proves:

1. the reordering of the code does not change the algebraic schedule
2. `c2` uses `ExactAdd` only on one-step subtraction classes
3. `c0` and `c2` final unrepaired outputs are both much narrower than the
   original raw `LNeg4_Pos4_HNeg4_Pos4` envelope

## Where The Proof Still Stopped

The current worksheet did **not** fully prove Sparse5B correct for all valid
inputs.

The unresolved obligations are:

1. `c0_sub_t2`
   - input class: `LNeg2_Pos1_HNeg2_Pos1`
   - current independent interval model does not prove that one `ExactAdd`
     repair always returns this class to `W`

2. `c1_uplus1_t2.im`
   - input class: `LNeg0_Pos2_HNeg0_Pos2`
   - current coarse model does not prove that the existing `ExactAdd`
     implementation always normalizes this positive-side split class

3. `c1_add_t3`
   - input class: `LNeg2_Pos2_HNeg2_Pos2`
   - this is the strongest unresolved checkpoint in the schedule

That gap turned out to be real rather than just a proof artifact:

- the structured adversarial search found a concrete valid-input failure on the
  `c0` path
- so `Sparse5B` is now known incorrect, not merely unproved

So the answer to "are the repair steps enough?" is now:

- **no for `Sparse5B`**
- **symbolically yes for the `c2` subpath**
- **insufficient on `c0` and `c1` under both the coarse model and actual
  execution**

## Per-`Fp` End-Only Model

The earlier sections tracked the actual sparse checkpoint schedule in
`Fp2`-sized chunks. For delayed repair design, the cleaner first-principles
model is to:

- track each `Fp` half independently
- logically inline `mul_by_u_plus_1`
- assume no repairs until the final end-of-block renorm

This strips away the checkpoint placement question and asks a simpler one:

- what final repair envelope does each unrepaired `Fp` half actually need?

### Contract Notation

Let:

- `W` be one strict widened `Fp` term
- `C(P,N)` be an exact signed sum of `P` positive `W` terms and `N` negative
  `W` terms

So:

- `W = C(1,0)`
- `W - W = C(1,1)`
- `W + W = C(2,0)`

After `v2_carry_from_biased(...)`, the low half is canonicalized into
`[0, R)`, and the carried high half lies in:

- `C(P,N) -> hi in [-N*p, P*p)`

That statement follows directly from the code-level contract of
`v2_carry_from_biased(...)` in
[Fp_avx_lufp_v2_impl.h](/c/Users/peter/repo/blst/benchmark/internal/Fp_avx_lufp_v2_impl.h#L327):

- the low limbs are carried and masked into radix range
- only the carried signed high residue remains

### `c0` With `mul_by_u_plus_1` Inlined

Start from:

```text
c0_pre = m0 - t1 - t2
c0 = (u+1)c0_pre + t0
```

with:

- `m0.re`, `m0.im`, `t0.re`, `t0.im`, `t1.re`, `t1.im`, `t2.re`, `t2.im` all
  in `W`

Inline the `Fp2` multiply by `(u+1)`:

```text
(u+1)(x_re + u x_im) = (x_re - x_im) + u(x_re + x_im)
```

So:

```text
c0.re =  m0.re - t1.re - t2.re
       - m0.im + t1.im + t2.im
       + t0.re
       = C(4,3)

c0.im =  m0.re - t1.re - t2.re
       + m0.im - t1.im - t2.im
       + t0.im
       = C(3,4)
```

After end-only carry:

- `c0.re -> hi in [-3*p, 4*p)`
- `c0.im -> hi in [-4*p, 3*p)`

So an end-only `c0` repair needs to handle the full envelope:

- `c0 in [-4*p, 4p)` up to branch asymmetry between real and imaginary halves

### `c1` With `mul_by_u_plus_1` Inlined

Start from:

```text
c1 = m1 - t0 - t1 + (u+1)t2
```

Inline:

```text
c1.re = m1.re - t0.re - t1.re + t2.re - t2.im
      = C(2,3)

c1.im = m1.im - t0.im - t1.im + t2.re + t2.im
      = C(3,2)
```

After end-only carry:

- `c1.re -> hi in [-3*p, 2*p)`
- `c1.im -> hi in [-2*p, 3*p)`

So an end-only `c1` repair needs to handle:

- `c1 in [-3*p, 3p)` up to real/imag asymmetry

### `c2` Without Intermediate Repairs

Start from:

```text
c2 = m2 - t0 - t2 + t1
```

Inline per `Fp`:

```text
c2.re = m2.re - t0.re - t2.re + t1.re
      = C(2,2)

c2.im = m2.im - t0.im - t2.im + t1.im
      = C(2,2)
```

After end-only carry:

- `c2.re -> hi in [-2*p, 2*p)`
- `c2.im -> hi in [-2*p, 2*p)`

So an end-only `c2` repair needs only:

- `c2 in [-2*p, 2p)`

### End-Only Repair Implications

This end-only per-`Fp` model exposes the repair problem much more cleanly than
the earlier coarse `Fp2` interval model.

The final unrepaired block requirements are:

- `c0`: repair family for `[-4*p, 4*p)`
- `c1`: repair family for `[-3*p, 3*p)`
- `c2`: repair family for `[-2*p, 2*p)`

That immediately explains why single-step checkpoint repair was too weak as a
general delayed-repair story:

- `ExactAdd` is only a one-step negative repair
- it is not an end-only repair for `c0`, `c1`, or even symmetric `c2`
  without more structure

It also suggests the correct synthesis order for future kernels:

1. derive a final repair ladder with a proven admissible source envelope
2. classify each block by the smallest final envelope it needs
3. only then decide whether mid-block checkpoint repairs are cheaper than a
   stronger final repair

## Practical Conclusion

Sparse5B is not "mystically correct." It is now known incorrect on a valid
input. But the bound story is still useful because it shows exactly why the
naive sparse-5 argument was too weak:

- the checkpoint model already isolated `c0` and `c1` as the unsafe parts
- the end-only per-`Fp` model shows the true target repair envelopes:
  - `c0`: `[-4*p, 4*p)`
  - `c1`: `[-3*p, 3*p)`
  - `c2`: `[-2*p, 2*p)`

But the current symbolic model is still missing exactly the same thing the
debug explainer identified:

- a repair-family proof tied directly to the actual code-level renorm steps,
  not just the old coarse `LNegA_PosB_HNegC_PosD` bookkeeping

For this kernel, the next real proof step is therefore not more sparse-mask
search. It is one of:

1. prove a final repair ladder on:
   - `[-2*p, 2*p)`
   - `[-3*p, 3*p)`
   - `[-4*p, 4*p)`

2. derive the cheapest checkpoint schedule whose intermediate states fit that
   ladder's admissible source classes

Until one of those is done, Sparse5B should be treated as:

- disproved by counterexample
- still useful as a bounds case study
- not a correct kernel candidate
