# Lean scalar curve API

The `catblst::scalar` module is an allocation-free Rust interface over
category-blst's scalar backend. Its type names and ordinary arithmetic surface
are inspired by `blstrs`, which makes migration from typed Rust curve code
direct. This is not a source-compatibility fork of `blstrs`.

The Cargo package is temporarily named `blst` for native-link compatibility.
New Category consumers should alias it as `catblst`; examples and new source use
the intended namespace.

## Types

`catblst::Scalar` is a transparent wrapper over the native Montgomery-form
`blst_fr`. It is shared by scalar and AVX protocol implementations.

The scalar backend owns these distinct transparent wrappers:

- `catblst::scalar::G1Projective` over `blst_p1`
- `catblst::scalar::G1Affine` over `blst_p1_affine`
- `catblst::scalar::G2Projective` over `blst_p2`
- `catblst::scalar::G2Affine` over `blst_p2_affine`
- `catblst::scalar::MillerLoopResult` over a general `blst_fp12`
- `catblst::scalar::Gt` over a final-exponentiated target-group `blst_fp12`

The corresponding native AVX point and target-group types remain under
`catblst::avx`. Scalar and AVX point representations have different sizes and
must not be reinterpreted or implicitly converted.

## Supported surface

`Scalar` provides:

- `ZERO`, `ONE`, and conversion from `u64`
- addition, subtraction, multiplication, negation, assignment, and sums
- inversion with explicit zero rejection
- strict canonical 32-byte big-endian and little-endian decoding
- canonical 32-byte big-endian and little-endian encoding
- 64-byte big-endian and little-endian reduction
- direct `blst_fr` borrows and an explicitly unsafe raw constructor

G1 and G2 types provide:

- identity and generator constructors
- complete addition, subtraction, doubling, and negation
- fixed-width constant-time multiplication for possibly secret scalars
- explicitly named variable-time multiplication for public scalars
- projective/affine conversion
- canonical compressed and uncompressed encoding
- strict decoding followed by prime-order subgroup validation
- identity, curve, and subgroup observations
- direct native borrows and explicitly unsafe raw constructors

Pairing types provide:

- single and allocation-free batched Miller loops
- a separate `MillerLoopResult` domain before final exponentiation
- final exponentiation into `Gt`
- GT identity, generator, addition, subtraction, negation, doubling, and scalar
  multiplication under the additive conventions used by `blstrs`
- strict canonical GT encoding and decoding with coefficient and subgroup checks
- direct native borrows and explicitly unsafe raw constructors

## Deliberate omissions

This initial surface does not provide:

- `ff`, `group`, or `pairing` trait implementations
- Rayon, worker threads, thread pools, or asynchronous work
- hidden heap allocation
- Serde integration
- generic hash-to-curve abstractions
- implicit scalar/AVX point or GT conversions
- a promise of complete `blstrs` API compatibility

Random-byte generation remains application-owned. Callers that require an
exactly uniform scalar use rejection sampling with the strict 32-byte decoder.
Protocols that specify wide reduction pass 64 bytes to the corresponding
constructor. This keeps entropy policy and RNG dependencies outside the
bottom-of-stack curve layer.

## Backend boundary

A protocol instance selects either scalar or AVX native point types and retains
that representation through its hot path. Both backends must produce identical
canonical encodings. Conversion through canonical bytes is appropriate for
serialization, hashing, conformance tests, and explicit backend selection; it
is not an internal arithmetic strategy.

Raw accessors are zero-conversion borrows. Safe constructors establish the
wrapper's mathematical invariant. `from_raw_unchecked` is reserved for outputs
of native kernels whose contracts already establish that invariant.
