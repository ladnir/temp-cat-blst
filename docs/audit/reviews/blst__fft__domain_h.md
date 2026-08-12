# File Review: `blst/fft/domain.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/fft/domain.h` |
| Roadmap phase / queue | Phase 17 FFT support-header packet |
| Reviewed commit | `0fc26d2617f74b6ef82345fee4f4400979b1534d` |
| Reviewer | Codex audit team |
| Date | 2026-08-11 |
| Disposition | `in_scope_production` |
| Subsystem | Generic root construction and reusable prepared FFT domains |
| Profiles | scalar and selected AVX native C++; C/Rust Fr adapters |
| Generation/provenance | Locally composed handwritten C++ templates; Apache-2.0 and Category Labs |
| Upstream references | `not_applicable: upstream BLST has no corresponding FFT domain helper` |
| Related files | `blst/fft/ctFft.h`, `blst/fft/primitveRoot.h`, `blst/fft/context.h`, `blst/fft_c_api.cpp`, `bindings/blst_fft.h`, `bindings/rust/src/fft.rs`, `test/Fft_api_tests.cpp` |

## Purpose and trust boundary

This header constructs ordered Fr root tables and the retained prepared-domain
layout: `n` forward roots, `1/n`, and an Fr-encoded `n` tag. Transform sizes,
direction, and domain selection are public. Values transformed through a domain
may be secret-derived, but this file never branches on them.

The templates are narrow native helpers. Supported C and safe Rust adapters
validate untrusted lengths, exact storage, and the size tag before entry and
translate exceptions at their ABI boundaries.

## Contracts and invariants

- `n` is a nonzero supported power of two within the signed FFT index range.
- `make_roots` rejects an invalid `n` before allocation and may throw
  `std::invalid_argument` or an allocation exception.
- `fill_roots` receives exactly one supported root-table length and writes
  consecutive powers beginning at one. Its direction is selected by `inverse`.
- `initialize_domain` receives exactly `n + 2` elements. A failure from a
  hypothetical generic arithmetic context may leave partial preparation, so a
  caller publishes the domain only after return.
- `transform_with_domain` receives a domain initialized for exactly the value
  count. Values and domain do not overlap; domain storage remains immutable.
- Debug builds diagnose native shape violations. Hot Release helpers do not
  revalidate shapes, contents, the encoded tag, mathematical root order, or
  overlap. The C/Rust boundary performs its documented Release validation.

## Dependency and call-path map

`primRootOfUnity<Fr>` supplies the known BLS12-381 power-of-two root table.
`FrCoeffCtx` supplies field multiplication, inversion, and integer conversion.
`fft_with_forward_roots` consumes the stored roots and scale. Production callers
are the Fr C API, interpolation, and middle-product root construction; Rust calls
through the C adapter. Frozen BTx/BTE compiled as an unchanged dependency and
was not audited or modified.

## Upstream/specification comparison

Upstream BLST has no generic FFT domain facility. The layout is local CatBLST
design already exposed by `bindings/blst_fft.h`; this review preserves the exact
ABI-visible `n + 2` representation and root order.

## Checklist evidence

### 1. Identity and purpose

The complete file, history, production callers, supported adapters, tests, and
read-only frozen-source references were inventoried.

### 2. Contracts and invariants

Size, direction, layout, mutation, publication, overlap, exception, and domain-
content preconditions are stated above and in source.

### 3. Upstream and specification lineage

The helper is local and its Apache-2.0/Category Labs provenance is explicit.

### 4. Correctness

Forward prepared transforms match ordinary transforms and invert back to their
input. Root construction reuses the audited primitive-root table and shared FFT
size predicate.

### 5. Memory and lifetime safety

The allocating builder rejects invalid shapes before vector construction.
Exact Debug span checks precede every prepared-domain index. Returned vectors
own their storage; executors retain no pointers.

### 6. Undefined and implementation-defined behavior

The supported-size predicate prevents invalid signed FFT indices. Native Release
misuse remains a documented narrow-contract violation; public adapters validate
before span construction and indexing.

### 7. Side channels and cryptographic hygiene

Loops and addresses depend only on public domain length and direction. Fr
arithmetic supplies coefficient-value timing behavior. No data-dependent root
selection was added.

### 8. Randomness and secret lifecycle

N/A: the Fr specialization uses fixed roots and no entropy. The helper owns no
secret copy beyond caller-provided transformed data.

### 9. API and ABI

Templates export no binary ABI. C and Rust retain the same prepared-domain size,
initialization, transform, error, alias, and publication contracts.

### 10. Concurrency

There is no mutable state in this header. Immutable completed roots/domains may
be shared; callers synchronize mutable transform buffers.

### 11. Resource exhaustion and denial of service

Valid large root vectors can consume substantial memory. Invalid sizes fail
before allocation; supported public adapters report allocation failure.

### 12. Portability and build behavior

The header compiles standalone under MSVC 19.50. Scalar-default and portable-
AVX-selected Debug configurations instantiate its retained uses. CUDA is out of
scope.

### 13. Performance-sensitive structure

Root preparation remains one linear multiplication loop. Transform execution
adds only `BLST_ASSERT` shape checks, which disappear under `NDEBUG`; no Release
check, dispatch, allocation, or wrapper was added to a butterfly path.

### 14. Tests and evidence

Invalid allocation-bearing sizes throw, malformed hot shapes terminate in Debug,
and prepared forward/inverse transforms agree with the ordinary engine. The
selected scalar and portable-AVX Debug objects compile.

### Specialized additions

The encoded tag is a shape guard at supported ABI boundaries, not authentication
of arbitrary domain contents. Validating every root would be linear extra work
and remains outside the narrow native contract.

## Findings

`AUD-293`, `AUD-300`, and `AUD-304`.

## Validation

- MSVC 19.50 scalar-default Debug focused FFT build and 9-test context/domain
  run: passed.
- Scalar-default and portable-AVX-selected Debug four-test domain/context run:
  4/4 passed in each profile.
- Standalone `domain.h` compile with `/std:c++latest /EHsc /c`: passed.
- No benchmark, CUDA, full suite, sanitizer, native AVX-512, remote, or Apple
  run was performed.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Earlier Phase 12 database | `AUD-293` context; no new diagnostic | Does not establish mathematical root order |
| Clang Static Analyzer / clang-tidy | Earlier retained profiles | No mapped result | Current exact commit used focused MSVC validation |
| Clippy | N/A | Rust adapter reviewed separately | Native template header |
| ASan/UBSan | Earlier Phase 17 focused run | No mapped prepared-domain report | Not rerun at this exact commit |
| Fuzzing | Not yet run | Phase 17 residual | Malformed public shapes have direct tests; domain-content fuzzing remains |

## Residual risk and decision

Arbitrary native domain contents and overlap are preconditions rather than hot
Release checks. Native AVX-512 hardware, sanitizer, fuzzing, and Apple milestone
evidence remain pending. The file is complete for line-by-line hand review.
