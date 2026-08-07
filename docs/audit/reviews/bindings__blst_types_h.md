# File Review: `bindings/blst_types.h`

## Record

| Field | Value |
| --- | --- |
| Path | `bindings/blst_types.h` |
| Roadmap phase / queue | Phase 18 public C/C++/Rust surfaces, Wave 1 representation packet |
| Reviewed commit | `aa4c7ed187f5a63496f0e94f9284108cbd46c7f1` |
| Reviewer | Codex audit team |
| Date | 2026-08-06 |
| Disposition | `in_scope_production`, complete |
| Subsystem | public scalar, field, point, and common binding representations |
| Profiles | C99-or-later, native C++20, scalar, AVX declarations, Rust FFI |
| Generation/provenance | Handwritten split of the upstream public binding header; native C++ representations remain owned by their core headers |
| Upstream references | `54e6e55674722fc2797ebb4bbb71b26d881eb4b8:bindings/blst.h`; local split and representation hardening history |
| Related files | `bindings/blst_scalar.h`, `bindings/blst_avx.h`, `blst/common/error_codes.h`, scalar field/point owner headers, `bindings/rust/src/bindings.rs`, `test/audit_scalar_header.c`, `test/BindingHeader_tests.cpp` |

This record reviews the exact implementation commit above.

## Purpose and trust boundary

This header owns the fixed public C storage vocabulary shared by every split
binding header. It defines byte and limb names, raw scalar storage, scalar field
and base/extension field storage, and G1/G2 projective and affine storage. In C
mode these are data-only structures. In C++ mode the public names are exact
aliases of the native scalar implementation types, avoiding wrapper objects,
copies, conversions, allocations, or runtime dispatch at the C boundary.

The types may carry secrets (`blst_scalar` and field values) or attacker-
controlled decoded representations (field and point storage). The storage type
alone does not validate reduction, nonzero secret-key status, curve membership,
or subgroup membership. Those properties belong to the operation that creates
or validates the value and to the precondition of every later consumer.

## Contracts and invariants

- CatBLST requires C99 or later. Pre-C23 C obtains `bool` from `<stdbool.h>`;
  C23 and C++ use their language keyword.
- `blst_byte` is exactly `uint8_t`; `blst_limb_t` is exactly `uint64_t`.
- `blst_scalar` is exactly 32 bytes with byte alignment and may contain any
  256-bit little-endian integer. `blst_scalar_fr_check` establishes `[0,r)`;
  `blst_sk_check` establishes the secret-key range `[1,r)`.
- `blst_scalar` and `blst_fr` are representation-distinct even though both are
  32 bytes. The former is an ordinary little-endian integer; the latter is a
  reduced Montgomery field representation.
- C field/point structures contain their documented nested components in
  declaration order. Their expected sizes are 32, 48, 96, 288, 576, 144, 96,
  288, and 192 bytes for Fr, Fp, Fp2, Fp6, Fp12, P1, P1 affine, P2, and P2
  affine respectively, with limb alignment.
- C++ public names are exact aliases of the native `blst::scalar` owners. The
  header asserts type identity, standard layout, trivial copying, size,
  alignment, and leading-member offsets.
- Raw field storage does not itself establish a reduced Montgomery value. Raw
  point storage does not itself establish on-curve or prime-order-subgroup
  status. Infinity and other representation semantics are documented with the
  owning operation declarations and native types.
- Callers own all storage and lifetimes. This header defines no object with
  dynamic ownership, destructor work, initialization side effect, or mutable
  global state.

The authoritative C operation declarations are in the split binding headers,
principally `bindings/blst_scalar.h`. The authoritative native C++ contracts are
in the aliased scalar field and point headers. Checked-in Rust `repr(C)` storage
and safe-wrapper contracts are generated from and layered over this C surface;
this file implements no independent safe Rust operation.

## Dependency and call-path map

`blst_scalar.h`, `blst_avx.h`, `blst_middle_product.h`, and `blst_fft.h` include
this header. It imports the shared error representation and, in C++ mode, the
native Fr/Fp/Fp2/Fp6/Fp12 and point owners. Scalar exports, retained C adapters,
FFT, middle product, AVX conversion boundaries, tests, and Rust FFI all consume
the resulting names.

Opaque pairing and uniqueness contexts are deliberately declared by
`blst_scalar.h`, not here. `AUD-187` gives those unrelated objects distinct
incomplete tags so a C/C++ compiler rejects accidental cross-use.

The ordinary public C ABI always uses scalar-layout field and point storage.
Selected or explicit AVX operations either consume those scalar boundary types
or use separately aligned opaque AVX shells. `BLST_AUX_ALIGNAS` supports those
shell declarations but creates no object or runtime instruction here.

## Upstream/specification comparison

Upstream BLST keeps these layouts in `bindings/blst.h`. CatBLST split them into
`blst_types.h`, replaced generic global `byte` and `limb_t` names with
library-prefixed typedefs, requires C99, and makes the supported native C++
types the exact C++ aliases. The physical scalar, field, and point layouts
remain upstream-compatible.

CatBLST intentionally clarifies that upstream's `blst_scalar` storage is not a
type-level proof of reduction or secret-key validity. `AUD-188` corrected the
former “canonical scalar” wording without changing any representation or
function. The common error definition moved to a neutral internal/public
header under `AUD-064`; `AUD-189` completes its width and Rust
evolution contract separately.

## Checklist evidence

### 1. Identity and purpose

The file, exact commit, Phase 18 position, handwritten provenance, upstream
source, supported profiles, and all direct public representation families are
recorded above.

### 2. Contracts and invariants

Byte order, representation state, range-establishing validators, Montgomery
distinction, C/C++ type ownership, sizes, alignments, offsets, lifetime, and
validation-state limitations are explicit. Operation-specific aliasing and
error behavior remain with their declarations rather than being generalized
incorrectly at the storage layer.

### 3. Upstream and specification lineage

The layouts and original error values were compared with the pinned upstream
`bindings/blst.h`. Intentional divergences are the split header, prefixed public
typedefs, C99 minimum, exact native C++ aliases, compile-time layout assertions,
and corrected scalar-state vocabulary.

### 4. Correctness

All fixed-width array extents resolve to the expected byte counts with a
64-bit limb. Nested structures reproduce the Fp2/Fp6/Fp12 and G1/G2 coordinate
composition. `AUD-188` now agrees with the actual fixed-width import helpers,
which preserve arbitrary 256-bit values rather than silently reducing them.

### 5. Memory and lifetime safety

Every type is caller-owned fixed storage. The C++ arithmetic/point aliases are
asserted trivially copyable and standard-layout. The header performs no access,
construction in raw storage, allocation, deallocation, cast, or pointer
arithmetic. Passing malformed object contents to a narrow arithmetic operation
remains a caller contract violation, not type-level validation.

### 6. Undefined and implementation-defined behavior

Fixed-width standard integer types avoid host `long` width differences. C++
type identity removes the prior need to reinterpret one live C++ field/point
type as another. The C and C++ declarations necessarily rely on the supported
platform C ABI for cross-language layout, backed by exact compile-time sizes,
alignments, and offsets. Unsupported AVX alignment-macro fallback behavior is
deferred to the explicit AVX public-header review.

### 7. Side channels and cryptographic hygiene

N/A for executable behavior: typedefs, aliases, and assertions generate no
runtime path. The header neither clears nor duplicates secret storage and makes
no blanket constant-time promise for downstream operations.

### 8. Randomness and secret lifecycle

N/A for randomness. `blst_scalar` can carry secret material, but ownership,
generation, copying, validation, and erasure belong to key-generation and
language-wrapper contracts. This header creates no secret and owns no cleanup.

### 9. API and ABI

The sizes and alignment are public ABI. C++ names are exact native aliases and
therefore add no boundary wrapper. Rust mirrors the C structures with
`repr(C)` and layout tests. All public identifiers are library-prefixed except
the deliberately retained upstream `BLST_ERROR` family. The include guard is
non-reserved and project-prefixed.

### 10. Concurrency

N/A: there is no state, initialization, caching, thread-local storage, or
synchronization.

### 11. Resource exhaustion and denial of service

N/A: all storage is fixed-size and the header performs no work or allocation.
Variable-size operation and scratch calculations are reviewed with their
owners.

### 12. Portability and build behavior

Strict MSVC C11 and C++20 standalone public-header probes pass. The C++ probe
uses the generated configuration include directory required by the native type
owners. The header rejects C before C99. Current production support is 64-bit,
little-endian MSVC/GCC/Clang-family hosts; other platform evidence remains a
milestone concern.

### 13. Performance-sensitive structure

All abstractions are compile-time aliases or aggregates. There is no copy,
conversion, heap use, virtual dispatch, dynamic type, runtime backend dispatch,
or added hot-path instruction. The hardening commit changed documentation and
compile-time guards only for this file.

### 14. Tests and evidence

Strict standalone probes cover language legality and include dependencies.
`BindingHeader.*` covers sizes, declarations, common scalar conversions, and
linked public behavior. Checked-in Rust binding layout tests cover the FFI
storage. Broader mathematical correctness belongs to the field/curve reviews.

### Specialized additions

N/A for arithmetic proofs, generated assembly, fuzz parsing, and FFT
mathematics. AVX relevance is limited to the alignment helper and scalar/opaque
boundary layout; explicit AVX types retain their own review.

## Findings

- `AUD-188` - validated: raw 256-bit scalar storage is no longer described as
  inherently canonical or secret-key-valid.
- `AUD-187` - related and validated in `blst_scalar.h`: the two opaque context
  families now have distinct incomplete tags.
- `AUD-064`, `AUD-069`, and `AUD-189` are related representation
  findings whose canonical enum and C++/Rust boundaries are recorded in their
  owning reviews.

## Validation

- MSVC strict standalone C probe, normalized command:
  `cl /std:c11 /W4 /WX /Zs <source-and-generated-include-options> header_c.c`:
  passed.
- MSVC strict standalone C++ probe, normalized command:
  `cl /std:c++20 /permissive- /W4 /WX /Zs <source-and-generated-include-options> header_cpp.cpp`:
  passed.
- `blst_tests --gtest_filter=BindingHeader.*`: 17/17 passed.
- The full checked-in Rust bindings compiled; the targeted error-layout tests
  listed in the `error_codes.h` review passed.
- No broad build, benchmark, sanitizer, CUDA, BTx/BTE audit, Apple CI, or
  portable-AVX Release build was run for this packet.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Existing Phase 12 scalar databases | No unresolved finding mapped to these declarations | Does not prove cross-language ABI or documentation semantics |
| Clang Static Analyzer / clang-tidy | Existing retained-profile baselines | No unresolved declaration diagnostic | Exact C/C++ layout is covered by compiler assertions instead |
| Clippy | Full checked-in Rust bindings compile | No new diagnostic | Clippy does not prove native C layout equivalence |
| ASan/UBSan | Not rerun for this header-only packet | N/A | No executable behavior; downstream typed access retains separate coverage |
| Fuzzing | Not applicable to fixed declarations | N/A | Parsers and validators are reviewed with their implementations |

## Residual risk and decision

Complete for this header's exact storage and alias contracts at the reviewed
commit. Cross-language ABI evidence remains compiler/platform-specific and must
be refreshed for new supported toolchains. The unknown-compiler empty fallback
for `BLST_AUX_ALIGNAS` is intentionally left to the explicit AVX public-header
review. This completion does not validate arbitrary raw field/point contents or
the mathematical, timing, aliasing, and failure contracts of their consumers.
