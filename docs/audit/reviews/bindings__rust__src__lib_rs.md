# File Review: high-level Rust root API

## Record

| Field | Value |
| --- | --- |
| Path | `bindings/rust/src/lib.rs` |
| Roadmap phase / queue | Phase 18 safe Rust root and signature surface |
| Reviewed commit | `778f2f60d1fd4ba19569ba0dd1d4e3f5f1a50b88` |
| Reviewer | Codex audit team |
| Date | 2026-08-12 |
| Disposition | in-scope production |
| Subsystem | crate root, BLS signature wrappers, aggregation, and multipoint adapters |
| Profiles | no-std/alloc or std Rust; scalar or native AVX implementation beneath the scalar ABI |
| Generation/provenance | Handwritten safe layer with upstream BLST lineage; raw declarations included separately |
| Upstream references | Common base `f48500c1fdbefa7c0bf9800bccd65d28236799c1`; compatibility retained while contracts and ownership are hardened locally |
| Related files | `pairing.rs`, `pippenger.rs`, `primitives.rs`, `fft.rs`, `avx.rs`, full-binding guide, C scalar header |

## Purpose and trust boundary

The root exposes raw POD/FFI names and, behind `full-rust-bindings`, the two BLS
orientations, secret/public/signature wrappers, aggregates, verification, and
batch point adapters. Serialized adversarial bytes enter through decoding;
secret material enters through key derivation or validated scalar decoding.
Wrapper types deliberately do not carry a runtime validation-state tag.

## Contracts and invariants

- Decode/uncompress establishes canonical on-curve storage but not subgroup
  membership. `key_validate`/`sig_validate` establish the documented subgroup
  and infinity postconditions. A false validation flag means prior validation.
- Raw-to-wrapper `From` conversions perform no mathematical validation.
- Ordinary aggregation, verification, and `MultiPoint` operations require
  prime-order subgroup operands and are variable-time over public data.
  Scalar layout, length, allocation, error, and panic behavior is documented.
- Fast aggregate verification requires prior proof of possession. Batch
  verification requires independent unpredictable nonzero public weights.
- `SecretKey` is always canonical and nonzero, has no `Default`, constructs FFI
  outputs in private uninitialized raw storage, validates legacy v3's negligible
  zero case, and zeroizes every owned copy on drop.
- Slice representation conversions are valid because `PublicKey`/`Signature`
  are transparent over their affine carriers; each unsafe conversion records
  that layout argument locally.
- Pairing/aggregation is synchronous on the calling thread. No global pool or
  dynamic type dispatch is introduced.

## Checklist evidence

1-3. Identity/contracts/lineage: public methods and carriers are documented
with shared vocabulary; upstream names remain compatible.
4. Correctness: key derivation/serialization, signature verification,
aggregation, non-subgroup decoding, and multipoint behavior have focused tests.
5-6. Memory/UB: transparent conversions have explicit layout proofs;
`MaybeUninit` values are published only after native output/validation.
7-8. Side channels/secrets: secret constructors preserve zeroization; public
MSM/batch paths are labelled variable-time; no blanket zeroing was added.
9. API/ABI: scalar FFI types remain the boundary; validation flags and raw
conversions state exactly what they do not check.
10. Concurrency: all work is caller-thread synchronous.
11. Resources: owned pairing/MSM allocations and panic/error boundaries are
documented; no hidden scheduler exists.
12. Portability: no-std alloc plus Serde-secret check passes; custom cfg is
declared rather than globally suppressed.
13. Performance: no Release validation, allocation, or dynamic abstraction was
added to hot native paths.
14. Tests: rustfmt, rustdoc with warnings denied, structural contract gate, and
focused high-level tests pass.

## Findings

`AUD-085`, `AUD-088`, `AUD-089`, `AUD-094`, `AUD-098`, `AUD-149`, `AUD-155`,
`AUD-159`, `AUD-160`, `AUD-161`, `AUD-171`, `AUD-321`, `AUD-330`, `AUD-331`,
`AUD-332`.

## Validation

- `cargo fmt --check`: passed.
- `RUSTDOCFLAGS='-D warnings' cargo doc --no-deps --features
  full-rust-bindings`: passed.
- `check_rust_api_contracts.ps1`: documented public surface, non-default
  `SecretKey`, and six justified slice conversions passed.
- Secret-key invariant test: 2/2 orientations, covering v3/v4/v4.5/v5,
  EIP-2333 master/child, round trip, zero, noncanonical, and wrong length.
- `BLST_TEST_NO_STD=1 cargo check --features serde-secret` with warnings denied:
  passed.

## Automated evidence

| Tool/profile | Result | Blind spots |
| --- | --- | --- |
| rustdoc/structural contract gate | passed | Semantic quality still depends on this hand review |
| Rust focused tests | passed | Not a full release/platform sweep |
| Clippy baseline | Prior Phase 12 findings triaged | Does not inspect native C++ reached by FFI |
| Sanitizer/fuzzing | Owning C boundaries have focused evidence | Broader adversarial signature/state fuzzing remains a milestone campaign |

## Residual risk and decision

Validation state remains documentary by deliberate compatibility decision; a
caller can misuse a false validation flag without violating Rust memory safety.
`AUD-085` retains the future checked/unchecked deserialization and approved
batch-subgroup-validation design work. Those are explicit API evolution items,
not undocumented behavior in this reviewed file. The current root API is
complete at the reviewed commit.
