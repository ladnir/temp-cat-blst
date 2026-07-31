# Phase 11 Conformance Vectors and Release Language

## Record

- Date: 2026-07-31
- State: complete
- Upstream rows: `ff66582`, `67fc54b`
- Scope: retained C++ and Rust conformance tests plus exact README language
- Benchmarks: excluded and not run
- CUDA: experimental and untouched
- BTx/BTE: frozen and source-untouched

## Decision

The useful intent of `ff66582` is independent conformance evidence, not its
Python runner. The draft-13 BLS12-381 generator encodings and pairing
coefficients were therefore translated into retained C++ and Rust tests. The
literal `67fc54b` README edit was not imported: Zig remains absent and a
floating draft label was replaced with a narrow version-pinned statement.

## Pairing evidence

The generator and target-group coefficients come from
[`draft-irtf-cfrg-pairing-friendly-curves-13`, Appendix B.1](https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-pairing-friendly-curves-13#appendix-B.1).
The draft lists `e_0` through `e_11` in tower order; the tests map them into
BLST's serialized order and round-trip all 576 bytes through the public field
API.

The draft also specifies that production libraries may return
`e_lib(P,Q) = e_pseudocode(P,Q)^k`, with `k = 3` for BLS12-381 and BLST named
as an example. The tests deserialize the literal draft value, cube it with the
public Fp12 API, and compare that result with scalar and AVX production
pairings. They do not incorrectly demand byte equality between the raw draft
pseudocode value and BLST's production convention.

## Signature evidence

The BLS document is pinned to
[`draft-irtf-cfrg-bls-signature-07`](https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-bls-signature-07).
Its test-vector appendix is unfinished, so the fixed secret key and message
were cross-checked independently with
[`@noble/curves` 2.2.0](https://github.com/paulmillr/noble-curves/tree/043905dfc90e75fe77bbbd60c791baa1e5177801),
npm `gitHead` `043905dfc90e75fe77bbbd60c791baa1e5177801`. The resulting
compressed signatures match CatBLST in both basic-scheme orientations:

- signature in G1 with public key in G2 and the G1 `_NUL_` domain;
- signature in G2 with public key in G1 and the G2 `_NUL_` domain.

Both C++ and the feature-gated high-level Rust APIs reproduce the independent
bytes and successfully verify them. This establishes the tested ciphersuite
behavior; it does not turn an active Internet-Draft into an IETF standard or
claim that every optional scheme has an official draft vector.

## Focused validation

| Surface | Evidence | Result |
| --- | --- | --- |
| C++ scalar | Release scalar `Conformance.*` filter | Generator encodings, cubed production pairing, and both signature orientations passed; AVX test skipped by design |
| C++ AVX API | Debug AVX-portable `Conformance.*` filter | Scalar, AVX conversion/pairing, and both signature orientations passed |
| Rust scalar | `cargo test --test conformance` | Draft generator and cubed production pairing passed |
| Rust high-level signatures | `cargo test --test conformance --features full-rust-bindings` | Both independently generated signature orientations passed and verified |
| Rust AVX-portable | `cargo test --test conformance --features avx` on Windows | Generator encodings and production pairing passed through the enabled AVX-portable backend |
| Rust formatting | `cargo fmt -- --check` | Passed |

The checks were deliberately limited to the affected conformance targets.
Production arithmetic did not change, so a full build/test sweep was not
warranted in this change loop. An uncached MSVC release build proved too slow
for focused iteration and was replaced with debug behavior checks. This host
does not expose AVX-512, so native AVX-512 execution remains a release-matrix
item rather than being inferred from the portable backend.

## Residual risk

These vectors cover the canonical generators, one pairing, and one fixed
basic-signature input per group orientation. They complement rather than
replace broader property, malformed-input, constant-time, scalar/AVX
differential, and platform release testing. Native AVX-512 hardware and Apple
Silicon execution remain CI evidence items, and the final supported release
matrix remains a milestone gate.
