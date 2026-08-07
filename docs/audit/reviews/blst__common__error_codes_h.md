# File Review: `blst/common/error_codes.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/common/error_codes.h` |
| Roadmap phase / queue | Phase 18 public C/C++/Rust surfaces, Wave 1 representation packet |
| Reviewed commit | `aa4c7ed187f5a63496f0e94f9284108cbd46c7f1` |
| Reviewer | Codex audit team |
| Date | 2026-08-06 |
| Disposition | `in_scope_production`, complete |
| Subsystem | shared public/native error-code ABI |
| Profiles | C99-or-later, native C++, scalar, AVX, Rust FFI |
| Generation/provenance | Handwritten neutral owner extracted from upstream-compatible public error values |
| Upstream references | `54e6e55674722fc2797ebb4bbb71b26d881eb4b8:bindings/blst.h` error enum; local `BLST_BAD_SCALAR` extension |
| Related files | `bindings/blst_types.h`, `blst/common/errors.h`, `bindings/blst_scalar.h`, `bindings/blst_avx.h`, point decoders, pairing aggregation, `bindings/rust/src/bindings.rs`, `test/BindingHeader_tests.cpp` |

This record reviews the exact implementation commit above.

## Purpose and trust boundary

This header is the single source of truth for error results crossing the native
core, public C ABI, native C++ namespace bridge, AVX adapters, and Rust FFI. The
codes describe malformed attacker-controlled encodings, invalid curve/group
state, pairing-context mismatches, failed verification equations, forbidden
infinity, and malformed scalar arguments. The values contain no secret and may
be branched on or logged by callers.

Its trust boundary is representation-critical: every language must agree on
the four-byte calling convention and numeric assignments, and an older FFI
consumer must be able to preserve an error appended by a newer compatible
native library without creating an invalid language value.

## Contracts and invariants

- `BLST_SUCCESS` through `BLST_BAD_SCALAR` have explicit ABI values zero
  through seven.
- Values zero through six retain the upstream BLST assignments.
  `BLST_BAD_SCALAR == 7` is CatBLST's append-only extension.
- Existing numeric assignments never change. Future compatible values may be
  appended but may not be inserted or renumbered.
- `sizeof(BLST_ERROR)` is exactly four. A C99-compatible negative-array
  declaration rejects a native compiler configuration that chooses another
  enum width.
- `BLST_BAD_ENCODING` denotes malformed or noncanonical bytes;
  `BLST_POINT_NOT_ON_CURVE` denotes failure of the curve equation;
  `BLST_POINT_NOT_IN_GROUP` denotes an on-curve point outside a required
  prime-order subgroup.
- `BLST_AGGR_TYPE_MISMATCH` reports incompatible pairing orientations;
  `BLST_VERIFY_FAIL` reports a false verification equation;
  `BLST_PK_IS_INFINITY` reports an infinite public key where finite is
  required; `BLST_BAD_SCALAR` reports a scalar pointer, width, or value that
  violates that entry point's documented contract.
- Native functions return only a documented value. Failure-output mutation and
  validation work are operation-specific and are not implied by the enum.
- Rust uses a four-byte `#[repr(transparent)] BLST_ERROR(pub u32)` with
  associated named constants. Every 32-bit result is therefore a valid Rust
  value, including a future unknown code.

The authoritative C declaration is this global type. `blst/common/errors.h`
imports the exact type and constants into namespace `blst`; it does not define
a second C++ enum. Raw Rust declarations use the transparent wrapper. Safe Rust
wrappers return or compare the same value and can preserve unknown results.

## Dependency and call-path map

`bindings/blst_types.h` imports this header for every split C binding.
`blst/common/errors.h` imports it for native point decoding and aggregation.
Scalar and AVX point decoders produce encoding/curve/group errors; pairing
aggregation and verification produce aggregation, infinity, scalar, and verify
errors. AVX C adapters cast only between the exact native/public shared type.
Rust's checked-in raw declarations receive the four-byte result and safe
wrappers propagate or compare it.

Middle-product and FFT C APIs use Boolean/status-specific contracts rather than
`BLST_ERROR`; the `AUD-189` remediation removed the former inaccurate
middle-product claim while documenting the real result domain.
Experimental CUDA batch errors are a different enum outside this production
contract.

## Upstream/specification comparison

Upstream's public error enum supplies values zero through six. CatBLST first
centralized the public and native declarations under `AUD-064`, then appended
`BLST_BAD_SCALAR` for weighted-scalar argument failures. Explicit values retain
source and binary compatibility with upstream callers.

`AUD-189` adds a compile-time four-byte ABI requirement, defines append-only
evolution, and changes generated Rust from a closed enum to a transparent
integer wrapper. Associated constants preserve ordinary
`BLST_ERROR::BLST_SUCCESS` spelling and generate no conversion or dispatch.

## Checklist evidence

### 1. Identity and purpose

The neutral owner, exact commit, supported surfaces, producers, consumers, and
upstream lineage are recorded above. The header is handwritten and generated
bindings consume it.

### 2. Contracts and invariants

Numeric assignments, physical width, meanings, append-only evolution, native
producer obligation, and unknown-value FFI behavior are explicit. Individual
functions remain responsible for documenting which subset they return and
whether failure leaves output unchanged.

### 3. Upstream and specification lineage

Values zero through six were compared with the pinned upstream public header.
The value-seven extension, neutral ownership, four-byte guard, and forward-safe
Rust representation are intentional CatBLST divergences.

### 4. Correctness

All eight current constants have explicit consecutive values. The internal C++
bridge imports the exact global type. Scalar and AVX producers use the shared
constants, and Rust associated constants reproduce the same values. The
corrected documentation no longer claims a middle-product consumer that does
not exist.

### 5. Memory and lifetime safety

The header owns no storage beyond compile-time constants and defines no pointer
operation. The Rust transparent wrapper accepts every `u32`, preventing an
unknown native result from becoming an invalid Rust enum discriminant.

### 6. Undefined and implementation-defined behavior

C and C++ normally choose a four-byte enum here, but the underlying enum width
is otherwise implementation-defined. The typedef assertion makes the required
ABI a compile-time condition. Rust's physical representation is explicitly
four bytes and no longer relies on the closed valid-discriminant set of a Rust
enum.

### 7. Side channels and cryptographic hygiene

Error values are public outcomes. The header performs no computation, access,
clearing, or formatting and makes no constant-time claim for the operation that
produced an error.

### 8. Randomness and secret lifecycle

N/A: no randomness, secret value, storage, or cleanup. An error can reveal
which public validation stage failed; APIs that require uniform failure
behavior must translate it at their own boundary.

### 9. API and ABI

The numeric values and four-byte return representation are public ABI. C keeps
the upstream `BLST_ERROR` type and enumerator spellings. C++ imports that exact
type. Rust preserves the common associated-constant spelling while accepting
future values. The compile-time guard has no symbol or runtime cost.

### 10. Concurrency

N/A: immutable compile-time declarations only.

### 11. Resource exhaustion and denial of service

N/A: no work, allocation, loop, or attacker-sized input. Formatting or retry
policy belongs to consumers.

### 12. Portability and build behavior

The header is valid C99 and C++. Strict MSVC C11 and C++20 standalone public-
header probes pass with the four-byte assertion active. Non-four-byte enum
configurations are intentionally unsupported and fail during compilation.

### 13. Performance-sensitive structure

Constants and a transparent integer result have the ordinary four-byte ABI.
There is no allocation, copy, table, branch, exception, runtime translation, or
dynamic abstraction. Rust named constants compile to direct integer values.

### 14. Tests and evidence

Public-header tests assert `sizeof(BLST_ERROR) == 4` and exercise the current
endpoint value. Targeted Rust layout tests assert size four, alignment four,
and preservation of `0xffff_ffff`, proving the wrapper's unknown-value domain.
The full checked-in bindings compile after the representation change.

### Specialized additions

N/A for cryptographic arithmetic, assembly, FFT, fuzz parsing, or AVX code
generation. AVX relevance is only that its public adapters share these exact
values.

## Findings

- `AUD-064` - validated: public and internal declarations use one exact type.
- `AUD-189` - validated: the C/C++ ABI rejects a non-four-byte error enum;
  meanings/evolution are documented; and Rust can preserve unknown future
  results without invalid-enum undefined behavior.

## Validation

- MSVC strict standalone C probe, normalized command:
  `cl /std:c11 /W4 /WX /Zs <source-and-generated-include-options> header_c.c`:
  passed.
- MSVC strict standalone C++ probe, normalized command:
  `cl /std:c++20 /permissive- /W4 /WX /Zs <source-and-generated-include-options> header_cpp.cpp`:
  passed.
- `blst_tests --gtest_filter=BindingHeader.*`: 17/17 passed, including the
  four-byte C++ assertion and endpoint-value check.
- `cargo test --manifest-path bindings/rust/Cargo.toml bindgen_test_layout_blst_error`:
  passed.
- `cargo test --manifest-path bindings/rust/Cargo.toml bindgen_test_layout_blst_uniq_result`:
  passed; the complete checked-in bindings compiled in the targeted Cargo runs.
- No negative `-fshort-enums` probe, broad build, benchmark, sanitizer, CUDA,
  BTx/BTE audit, Apple CI, or portable-AVX Release build was run.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Existing Phase 12 scalar databases | `AUD-064` historical declaration-drift context; no unresolved finding here | Cannot establish FFI enum validity or future-version policy |
| Clang Static Analyzer / clang-tidy | Existing retained-profile baselines | No unresolved diagnostic | Enum width is configuration-dependent and covered by the compile guard |
| Clippy | Targeted Rust tests and full bindings compile | No new diagnostic | Native calling-convention equivalence still requires C/C++ evidence |
| ASan/UBSan | Not rerun for this declaration packet | N/A | No executable memory access in this header |
| Fuzzing | Not applicable to an integer result vocabulary | N/A | Decoder inputs are fuzzed/reviewed with their implementations |

## Residual risk and decision

Complete for the shared error representation at the reviewed commit. A
negative GCC/Clang `-fshort-enums` compile was not run; the source-level typedef
is nevertheless designed to reject that representation and the negative probe
should be added to the broader configuration matrix. Previously published Rust
artifacts that modeled `BLST_ERROR` as a closed enum cannot be retroactively
made forward-safe and should not be paired with a newer independently built
native library. New supported Rust builds use the transparent wrapper.
