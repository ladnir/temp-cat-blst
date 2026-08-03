# File Review: `blst/common/bytes.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/common/bytes.h` |
| Roadmap phase / queue | Phase 13, common representation and parsing foundation |
| Reviewed commit | `80ae0f9a23da1c30415498241e204297e15157af` |
| Reviewer | Codex audit team |
| Date | 2026-08-03 |
| Disposition | `in_scope_production` |
| Subsystem | Limb/byte endian conversion and hexadecimal conversion |
| Profiles | scalar, native C++, C API, Rust FFI; AVX callers through scalar boundaries |
| Generation/provenance | Hand-maintained C++ adaptation of upstream BLST `src/bytes.h` |
| Upstream references | Upstream head `54e6e55674722fc2797ebb4bbb71b26d881eb4b8`, `src/bytes.h` |
| Related files | `blst/exports.cpp`, scalar `Fp`/`Fr`, `bindings/blst_scalar.h`, `bindings/rust/src/bindings.rs`, `test/BindingHeader_tests.cpp` |

## Purpose and trust boundary

This header converts fixed-width limbs to/from byte and hexadecimal forms.
Hexadecimal text crosses public C/C++/Rust boundaries. Callers require exact
endianness, no overflow, strict whole-input parsing, and unchanged output on
failure.

## Contracts and invariants

- Supported hosts declare byte order at compile time; unknown byte order is a
  configuration error.
- Little-endian exact in-place conversion remains a no-op.
- Hex input may have an optional `0x`/`0X` prefix and must then contain one or
  more hexadecimal digits, consume the complete string, and fit the output.
- Parse helpers return `bool` and do not modify output on failure.
- Valid output is little-endian/canonical-limb form as required by its caller.
- Functions allocate nothing and do not throw.

## Dependency and call-path map

Scalar fields and serialization call the limb/byte helpers. The retained C
scalar/Fr/Fp hex APIs call the strict helpers and translate to Montgomery form
only after successful parsing. Rust exposes the same Boolean C ABI.

## Upstream/specification comparison

Core nibble and endian arithmetic matches upstream. CatBLST intentionally
replaces C/C++ union-punning endianness detection and changes test-facilitator
hex parsing into a strict supported API contract with a failure channel.

## Checklist evidence

### 1. Identity and purpose

All helper definitions and direct supported callers were inspected; benchmark,
BTx, and experimental CUDA consumers were excluded by policy.

### 2. Contracts and invariants

Endianness, prefix, digit, width, null, success, and output-publication rules
are listed above and tested at the public boundary.

### 3. Upstream and specification lineage

Compared against upstream `src/bytes.h`; strict whole-string parsing and
defined compile-time byte order are documented local hardening divergences.

### 4. Correctness

Digit count is bounded before multiplication by two, over-width strings are
rejected by the terminating-character check, and odd/even digit placement
preserves the existing little-endian representation.

### 5. Memory and lifetime safety

No allocation or retained pointer. Null inputs/outputs fail. The parser first
validates the complete string, then writes output, providing failure atomicity.

### 6. Undefined and implementation-defined behavior

The inactive-union-member read is gone. `defines.h` requires
`std::endian::little`, matching every supported target, so the exact in-place
path needs no runtime or preprocessor byte-order probe. Size doubling is
guarded against overflow.

### 7. Side channels and cryptographic hygiene

Text parsing is variable-time in public text length/content and is not a
constant-time secret-decoding API. Fixed-width binary field operations retain
their separate constant-time requirements.

### 8. Randomness and secret lifecycle

N/A: no entropy. Temporary parsed field values are stack objects; broader
secret-zeroization policy is reviewed at owning APIs.

### 9. API and ABI

Scalar/Fr/Fp hex entry points now return `bool` and are nodiscard in supported
headers. Existing binary callers that ignored the former void return retain
the same symbol/calling arguments; source callers must handle the result.

### 10. Concurrency

All functions are stateless and reentrant.

### 11. Resource exhaustion and denial of service

Scanning is bounded by output capacity plus one observed terminator byte. No
input-controlled allocation occurs.

### 12. Portability and build behavior

MSVC compilation and execution pass. GCC and Clang strict translation-unit
checks pass. Native Apple execution remains pending the next CI run.
Unsupported byte-order toolchains fail at the shared `std::endian` assertion.

### 13. Performance-sensitive structure

Binary conversion loops and the little-endian in-place fast path are retained.
Strict parsing adds validation only to text APIs, not field arithmetic kernels.

### 14. Tests and evidence

Tests cover valid lower/upper case, prefixes, empty/prefix-only, malformed
suffixes, sign/space, nulls, over-width input, exact-width Fr/Fp values, and
unchanged output on failure. Variable-length byte conversion matches an oracle.

### Specialized additions

Rust FFI declarations return `bool`; a focused Rust compile remains part of
the next retained-Rust packet. No assembly generation is owned here.

## Findings

`AUD-062` and `AUD-063` - fixed; Apple and Rust evidence gaps are recorded.

## Validation

- MSVC single-worker Debug target build: passed.
- Focused public/C ABI tests: 4/4 parsing/conversion cases passed.
- Strict malformed/null/width matrix preserves output on every failure.
- Supported callers were migrated and a repository search found no ignored
  production result.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Phase 12 scalar/portable databases | No remaining mapped diagnostic | Parsing policy and inactive-union UB required hand review |
| Clang Static Analyzer / clang-tidy | Retained profiles | No additional mapped diagnostic | Apple macro branch pending CI |
| Clippy | Declaration updated | `AUD-063` pending focused Rust compile | Generated FFI declaration has no safe Rust wrapper caller |
| ASan/UBSan | Prior scalar baselines | No mapped report | Does not prove unchanged-output API policy |
| Fuzzing | Planned parser/serialization expansion | N/A in this packet | Future fuzzing should compose public decoders and round trips |

## Residual risk and decision

The file audit is complete, while findings remain `fixed` until Apple CI and a
retained Rust compile confirm the two platform/language edges. The hexadecimal
APIs are intentionally strict and variable-time; callers must not use their
timing as a secret-dependent protocol signal.
