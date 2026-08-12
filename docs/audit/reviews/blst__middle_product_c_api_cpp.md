# File Review: `blst/middle_product_c_api.cpp`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/middle_product_c_api.cpp` |
| Roadmap phase / queue | Phase 17 Packet 17J |
| Reviewed implementation commit | `9e4ec59ff3ecf77a11c73688f47e3c9cc91c8ed9` (owner/header changed; adapter arithmetic remains `bdefea3`) |
| Reviewer | Codex audit team |
| Date | 2026-08-11 |
| Disposition | `in_scope_production` |
| Subsystem | Opaque C ownership and exception boundary for pairing middle products |
| Profiles | scalar boundary ABI with scalar or AVX native implementation |
| Generation/provenance | Locally composed C++ adapter; Apache-2.0 and Category Labs |
| Upstream references | `not_applicable: CatBLST-specific C API` |
| Related files | `bindings/blst_middle_product.h`, `blst/fft/middle_product.h`, safe Rust owner/guide, `test/MiddleProduct_c_api_tests.cpp` |

## Purpose and ABI boundary

The adapter owns one opaque heap-allocated B precompute and exposes creation,
evaluation, constant-term update, GT combination, and destruction. C storage is
always the scalar/public ABI layout. Static assertions prove exact identity,
standard layout, and trivial copyability for boundary G1, G2, and Fp12 types.
AVX-enabled builds convert input points explicitly and return GT through the
audited direct field bridge; no serialization/decode round trip remains.

## Safety and exception review

- Precompute validates `n` before constructing the `2n` span, translates every
  native/allocation exception to null, and publishes the opaque handle only
  after complete construction.
- Evaluation and update reject null or mismatched public inputs, compute into a
  private GT vector, and copy to C output only after native success. Failure is
  output-atomic.
- `combine` checks pointers and nonzero length before a nonallocating,
  nonthrowing elementwise Fp12 multiplication representing the GT group law.
  Exact aliasing with either input is supported and tested; partial overlap is
  a documented precondition violation.
- Destruction accepts null. A non-null pointer must be a live handle created by
  this API and must not race with evaluation or another destruction. Arbitrary
  or stale pointers cannot be validated by a C ABI and are now documented as
  outside the contract.
- Catch-all boundaries prevent C++ exceptions from crossing C or Rust. No
  pointer, span, reference, or exception object escapes the call lifetime.

## Cryptographic and concurrency contracts

Every point must already be a reduced member of G1/G2; the boundary performs no
Release membership validation. Outputs and combine inputs semantically belong
to GT despite using `blst_fp12` storage. Points and lengths are public and the
operations are variable-time. A completed handle is immutable and can be shared
by concurrent evaluations with distinct outputs; destruction must be excluded.
The adapter owns no global mutable state or thread pool.

## Performance and allocation

The C API intentionally stages output to guarantee failure atomicity and creates
a fresh native workspace per evaluation. Its documentation therefore continues
to state that evaluation/update may allocate on every call. The selected native
C++ API offers reusable workspaces and now avoids same-backend conversion
vectors; this does not alter the stable C ABI or its ownership model.

## Findings

`AUD-085`, `AUD-283`, `AUD-286`, `AUD-288`, `AUD-289`, `AUD-290`, `AUD-311`,
and `AUD-312`.

## Validation

- Scalar-default plus portable-AVX Debug focused middle-product suite: 22/22.
- Scalar-only Release focused middle-product suite: 14/14.
- AVX-default portable Debug selected middle-product checks: 2/2.
- Strict standalone C11 and C++20 public-header parsing at `/W4 /WX`: passed.
- Hostile-length reject-before-read, mismatched-length failure publication,
  scalar/AVX reconstruction, infinity, update, and both exact-alias combine
  directions are covered.
- No benchmark, CUDA, portable-AVX Release, full suite, sanitizer, native
  AVX-512 hardware, remote, Apple, or BTx/BTE audit/change was performed.

## Residual risk and decision

The C API cannot diagnose forged, stale, concurrently freed, overlapping, or
undersized pointers; these remain precise caller preconditions. Sanitizer,
fuzz, native AVX-512, and Apple evidence remain milestone gates. The containing
file is complete for line-by-line hand review.
