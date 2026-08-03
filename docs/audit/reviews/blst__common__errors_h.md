# File Review: `blst/common/errors.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/common/errors.h` |
| Roadmap phase / queue | Phase 13/18, common error ABI foundation |
| Reviewed commit | `80ae0f9a23da1c30415498241e204297e15157af` |
| Reviewer | Codex audit team |
| Date | 2026-08-03 |
| Disposition | `in_scope_production` |
| Subsystem | Canonical public error type and internal namespace exposure |
| Profiles | C, C++, Rust; scalar and AVX public adapters |
| Generation/provenance | Hand-maintained compatibility header derived from upstream BLST error codes |
| Upstream references | Upstream head `54e6e55674722fc2797ebb4bbb71b26d881eb4b8`, public `blst.h` error codes |
| Related files | `bindings/blst_types.h`, public C headers, `blst/c_api.cpp`, `blst/avx_c_api.cpp`, Rust bindings |

## Purpose and trust boundary

This header makes the public `BLST_ERROR` enum the single error type used by
internal C++ code. Error values cross C and Rust FFI and must never drift by
language or backend.

## Contracts and invariants

- The canonical public values are explicit integers zero through seven.
- `blst::BLST_ERROR` is the exact global `::BLST_ERROR` type, not a duplicate.
- Internal enumerator names alias the public enumerators.
- CUDA-specific batch errors remain a separate experimental enum.
- The header owns no state, allocation, exception, or control flow.

## Dependency and call-path map

Scalar point decoding, aggregation, verification, public C adapters, AVX
adapters, and Rust FFI all consume the canonical public type. Experimental
CUDA batch status is separate and excluded from the release claim.

## Upstream/specification comparison

Values zero through six retain upstream BLST ordering. CatBLST's public
`BLST_BAD_SCALAR = 7` is now explicit and visible internally, eliminating the
previous duplicate-type drift.

## Checklist evidence

### 1. Identity and purpose

Both former enum declarations and every supported error-returning boundary
were inventoried.

### 2. Contracts and invariants

Exact type identity and stable numeric values are compile-time assertions.

### 3. Upstream and specification lineage

The upstream ordering is preserved; the local scalar error extension is fixed
at value seven.

### 4. Correctness

There is no conversion table or duplicated declaration left to drift. Existing
switch/comparison code sees the same enumerators.

### 5. Memory and lifetime safety

N/A: enum aliases own no memory and retain no pointer.

### 6. Undefined and implementation-defined behavior

Cross-declaration C++ enum casts are no longer required for scalar internal
code. The C enum's implementation representation remains the supported ABI's
ordinary integer enum representation.

### 7. Side channels and cryptographic hygiene

Error codes are public control results. This file introduces no secret-
dependent computation.

### 8. Randomness and secret lifecycle

N/A.

### 9. API and ABI

Explicit values preserve the existing C/Rust ABI. Type identity is proven with
`std::is_same_v`; success and bad-scalar endpoints are value-asserted.

### 10. Concurrency

N/A: immutable compile-time declarations.

### 11. Resource exhaustion and denial of service

N/A.

### 12. Portability and build behavior

The public C-compatible enum and C++ aliases compile under the supported MSVC
profile. No platform-specific branch exists in the production portion.

### 13. Performance-sensitive structure

Compile-time aliases generate no code or storage.

### 14. Tests and evidence

Boundary static assertions and the complete MSVC scalar/AVX-default build pass.
Existing decoder/error-path suites continue to compile against the exact type.

### Specialized additions

Rust already exposes the same numeric enum values. CUDA batch errors are
explicitly outside this file's production disposition.

## Findings

`AUD-064` - validated.

## Validation

- `std::is_same_v<blst::BLST_ERROR, ::BLST_ERROR>`: passed.
- Explicit `BLST_SUCCESS == 0` and `BLST_BAD_SCALAR == 7`: passed.
- MSVC single-worker Debug target build: passed.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Phase 12 databases | No remaining mapped issue | Type duplication was an ABI design issue found by hand |
| Clang Static Analyzer / clang-tidy | Retained profiles | No mapped diagnostic | Cannot enforce cross-language numeric policy alone |
| Clippy | Retained Rust bindings | Numeric parity inspected | Rust compile refresh remains in the next Rust packet |
| ASan/UBSan | N/A | N/A | Compile-time declarations only |
| Fuzzing | N/A | N/A | No input parser or runtime operation |

## Residual risk and decision

Complete at the recorded commit. Future public error additions must receive an
explicit value and update every generated language binding in the same change.
