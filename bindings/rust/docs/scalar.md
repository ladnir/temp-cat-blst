# Lean scalar curve API

The `catblst::scalar` module is an allocation-free Rust interface over
category-blst's scalar backend. Its type names and operator conventions are
inspired by `blstrs`, which keeps migration of typed Rust curve code direct.
It is not a source-compatibility fork of `blstrs`.

The Cargo package is temporarily named `blst` for native-link compatibility.
New Category consumers should alias it as `catblst`; examples and new source use
the intended namespace.

## Contract model

Let `r` be the BLS12-381 scalar-field and prime-subgroup order. The safe types
carry these invariants:

- `catblst::Scalar` is a reduced element of `Fr = Z/rZ`, stored in the native
  Montgomery representation. Zero is valid.
- G1 and G2 values are on the corresponding curve and in its order-`r`
  subgroup. The group identity is valid.
- `catblst::scalar::MillerLoopResult` is a nonzero product of Miller-loop
  contributions evaluated on validated subgroup points. It is not a target
  group value before final exponentiation.
- `catblst::scalar::Gt` is in the order-`r` pairing target group. Its identity
  is valid.

Safe constructors establish these invariants, and safe arithmetic preserves
them. Consequently, ordinary arithmetic does not repeat curve, subgroup, or
target-group validation. Validation occurs at decoding and other explicit
trust boundaries. Observation methods such as `is_on_curve` and
`is_in_subgroup` are redundant for safe values; they are useful as defensive
assertions around unsafe native interoperation.

Native projective coordinates, native affine coordinates, Montgomery scalar
limbs, and native Fp12 limbs are arithmetic storage. They are not wire formats
and need not be unique representations of a mathematical value. Equality on
point wrappers compares mathematical points rather than coordinate limbs.

## Types and operator semantics

`catblst::Scalar` is a transparent wrapper over `blst_fr`. It is shared by
scalar and AVX protocol implementations. Scalar `+`, `-`, `*`, and unary `-`
are arithmetic in `Fr`. `Default` and an empty iterator sum return zero.

The scalar backend owns these distinct transparent wrappers:

- `catblst::scalar::G1Projective` over `blst_p1`
- `catblst::scalar::G1Affine` over `blst_p1_affine`
- `catblst::scalar::G2Projective` over `blst_p2`
- `catblst::scalar::G2Affine` over `blst_p2_affine`
- `catblst::scalar::MillerLoopResult` over a general `blst_fp12`
- `catblst::scalar::Gt` over a final-exponentiated target-group `blst_fp12`

For projective G1 and G2 values, `+`, `-`, and unary `-` are the curve group
operations, and `point * scalar` is scalar multiplication. `Default` and an
empty iterator sum return the identity.

`Gt` follows the additive convention used by `blstrs`, although its underlying
Fp12 group law is multiplicative:

- `a + b` multiplies their Fp12 representatives.
- `a - b` multiplies `a` by the inverse of `b`.
- `-a` computes the group inverse.
- `a * s` exponentiates the Fp12 representative by scalar `s`.
- `Default` and an empty iterator sum return the target-group identity.

Multiplication of two `MillerLoopResult` values combines their Miller-loop
products. It does not perform final exponentiation.

The corresponding AVX point and target-group types remain under
`catblst::avx`. Scalar and AVX point representations have different sizes and
must not be reinterpreted or implicitly converted.

## Preconditions, validation, and failure

| Operation | Preconditions | Validation and failure |
| --- | --- | --- |
| Scalar arithmetic | Inputs satisfy the `Scalar` invariant | No validation; safe operations preserve reduction |
| Strict scalar decoding | Exactly 32 bytes in the selected byte order | Accepts the unsigned integer exactly when it is in `[0, r)`; returns `None` otherwise; does not reduce |
| Wide scalar reduction | Exactly 64 bytes in the selected byte order | Every input succeeds and is reduced modulo `r` |
| Scalar inversion | Any valid scalar | Returns `None` for zero; otherwise returns the unique inverse |
| Point arithmetic | Safe subgroup-point inputs | No repeated curve or subgroup validation; identity inputs are valid |
| Constant-time point multiplication | Any safe point and scalar | Uses a canonical fixed-width 256-bit scalar encoding |
| Variable-time point multiplication | Scalar is public | Result is algebraically identical to constant-time multiplication |
| Point decoding | Exact type-level compressed or uncompressed size | Native decoding checks flags, canonical coordinates, and the curve equation; the wrapper then checks subgroup membership; the identity is accepted; errors are returned as `BLST_ERROR` |
| Single Miller loop | Safe G1 and G2 affine inputs | No repeated point validation; identity input contributes the empty product |
| Batched Miller loop | G1 and G2 slices have equal length | A mismatch returns `PairingError::LengthMismatch` before the native call; empty slices return the empty product |
| Final exponentiation | A safe, therefore nonzero, `MillerLoopResult` | No repeated provenance or nonzero check; returns `Gt` |
| GT decoding | Exactly 576 bytes | Rejects noncanonical base-field coefficients with `BLST_BAD_ENCODING`; rejects values outside the target group with `BLST_POINT_NOT_IN_GROUP`; accepts the identity |

A failed decoder exposes no partially constructed safe value. Decoder timing may
depend on the input, so encoded inputs must be treated as public. Validation
results establish only the documented algebraic invariant; protocol-level
authorization and proof verification remain the caller's responsibility.

## Encodings

Scalar encoders return the unique 32-byte integer in `[0, r)` in the selected
byte order. Strict decoders are suitable for parsing canonical protocol fields.
Wide constructors perform deterministic modular reduction; they are suitable
only when a protocol's scalar-derivation rule specifies wide reduction. They do
not implement strict decoding or rejection sampling.

G1 compressed and uncompressed encodings are 48 and 96 bytes. G2 compressed
and uncompressed encodings are 96 and 192 bytes. Each safe encoder produces the
canonical category-blst encoding, including for the identity.

GT encoding consists of twelve canonical 48-byte big-endian base-field
coefficients, for a total of 576 bytes. Category-blst orders these coefficients
by Fp2 index, then Fp6 half, then base-field component. This order is part of
the CatBLST wire contract and is independent of how `blst_fp12` is traversed in
memory.

Random-byte generation remains application-owned. A caller that requires an
exactly uniform scalar uses rejection sampling with the strict 32-byte decoder.
A protocol that specifies wide reduction passes 64 bytes to the corresponding
constructor. This keeps entropy policy and RNG dependencies outside the
bottom-of-stack curve layer.

## Timing and allocation

Scalar-field addition, subtraction, multiplication, and negation use native
constant-time entry points. For a nonzero scalar, the native inversion is also
intended to be constant-time; the safe `invert` method nevertheless exposes
whether its input is zero through its return value and initial zero test.
Unsuffixed point and GT multiplication encode a canonical scalar into a fixed
256-bit stack buffer and use the native constant-time entry point. A method
ending in `_vartime` permits input-dependent execution and therefore requires a
public scalar.

No constant-time guarantee is made for decoding, validation, encoding,
equality, identity tests, affine conversion, pairing, final exponentiation, or
any other operation unless its API documentation makes that guarantee. In
particular, a method's lack of a `_vartime` suffix is not by itself a timing
promise.

Cryptographic methods in this module use stack or caller-owned storage and do
not allocate. Batched Miller loops pass the contiguous transparent-wrapper
slices to one native batch call through two fixed-size stack descriptor arrays;
they do not construct pointer arrays proportional to the batch size. The types
contain no pointers and own no external resources. Independent operations on
disjoint mutable values may run concurrently under Rust's ordinary aliasing
rules. Safe methods express their aliasing requirements through Rust borrows.
An unsafe caller passing raw pointers onward to a native kernel remains
responsible for the native kernel's aliasing contract.

## Raw boundary

`as_raw` and `into_raw` expose native storage without conversion.
`as_raw_mut` and `from_raw_unchecked` are unsafe because the compiler cannot
verify their algebraic preconditions:

- A raw scalar must be reduced and in Montgomery form. Canonical integer limbs
  or a `blst_scalar` value are not valid `blst_fr` storage.
- A raw point must be a valid native representation of an on-curve point in the
  correct prime-order subgroup.
- A raw Miller-loop result must be nonzero and must come from valid subgroup
  point contributions. An arbitrary nonzero Fp12 value is insufficient.
- A raw GT value must be in the order-`r` target subgroup. An unfinalized
  Miller-loop result is insufficient.

A caller using `as_raw_mut` must restore the relevant invariant before the
borrow ends. An invalid wrapper can violate the preconditions of later native
operations. Native raw values should cross a serialization, hashing, or backend
boundary only through the canonical encoders. `Scalar`'s `Debug` implementation
redacts its value; debugging output is not a serialization format for any type.

## Supported surface

The lean API provides:

- scalar identities, arithmetic, inversion, strict decoding, wide reduction,
  canonical encoding, and native access;
- G1 and G2 identities, generators, complete group arithmetic, constant-time
  and variable-time scalar multiplication, affine/projective conversion,
  canonical encoding, strict decoding, subgroup validation, and native access;
- single and allocation-free batched Miller loops, explicit final
  exponentiation, target-group arithmetic, strict GT encoding and decoding, and
  native access.

## Deliberate omissions

This surface does not provide:

- `ff`, `group`, or `pairing` trait implementations;
- Rayon, worker threads, thread pools, or asynchronous work;
- hidden heap allocation;
- Serde integration;
- generic hash-to-curve abstractions;
- implicit scalar/AVX point or GT conversions;
- a promise of complete `blstrs` API compatibility.

## Backend boundary

A protocol instance selects either scalar or AVX native point types and retains
that representation through its hot path. Both backends must produce identical
canonical encodings. Conversion through canonical bytes is appropriate for
serialization, hashing, conformance tests, and explicit backend selection; it
is not an internal arithmetic strategy.
