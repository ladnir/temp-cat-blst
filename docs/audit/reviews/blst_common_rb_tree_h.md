# File Review: `blst/common/rb_tree.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/common/rb_tree.h` |
| Reviewed commit | `1e43e0e65c149ea58c999a27f243c7d676560823` |
| Reviewer | Codex |
| Date | 2026-07-30 |
| Disposition | `in_scope_production` |
| Subsystem | internal uniqueness ABI bridge |
| Profiles | scalar, AVX-512, C, Rust |
| Generation/provenance | Handwritten internal header introduced by `1e43e0e` |
| Upstream references | N/A: local typed bridge for the accepted `AUD-011` API |
| Related files | `blst/common/rb_tree.cpp`, `blst/c_api.cpp`, `bindings/blst_scalar.h`, `bindings/rust/src/bindings.rs` |

## Purpose and trust boundary

This internal C++ header gives `rb_tree.cpp` and `c_api.cpp` one typed contract
without exposing the tree representation. It owns only the internal result
enum and declarations. It must not define storage, allocate, or become a
second public API.

## Contracts and invariants

The five result values are explicit ABI integers 0 through 4. Storage crosses
the boundary as `void*` so alignment can be checked before conversion to the
private tree type. Size and alignment return `size_t`; initialization returns
success; insertion returns the exact result enum.

## Dependency and call-path map

Included only by `rb_tree.cpp` and `c_api.cpp`. `c_api.cpp` translates these
internal namespaced declarations into the public C declarations in
`bindings/blst_scalar.h`; bindgen produces the Rust declarations.

## Upstream/specification comparison

Upstream declares the three legacy functions directly around `src/rb_tree.c`.
This local header exists because the fork compiles the implementation as C++
inside namespace `blst` and now needs a checked four-function bridge. License
and SPDX provenance match the source implementation.

## Checklist evidence

### 1. Identity and purpose

Path, commit, provenance, consumers, and narrow responsibility are recorded.

### 2. Contracts and invariants

All declarations use fixed result codes, opaque storage, `size_t` lengths, and
const message data. Representation remains private.

### 3. Upstream and specification lineage

N/A as a standalone upstream file; it is a local declaration split required by
the accepted hardening of upstream `src/rb_tree.c`.

### 4. Correctness

Declarations exactly match definitions and the C bridge. Link and Rust FFI
tests exercise every function and every result except impractical depth 64.

### 5. Memory and lifetime safety

Using `void*` allows implementation-side alignment validation before forming a
private typed pointer. The header creates no object or ownership.

### 6. Undefined and implementation-defined behavior

No arithmetic, casts, shifts, or object access occurs in the header.

### 7. Side channels and cryptographic hygiene

N/A: declarations only; message-publicness is documented by the implementation
review.

### 8. Randomness and secret lifecycle

N/A: declarations only.

### 9. API and ABI

Enum values are explicit and agree with the public C header and generated Rust
enum. The header is C++-internal and is not installed as the supported C API.

### 10. Concurrency

N/A: no state. The implementation contract requires external synchronization
for one mutable tree.

### 11. Resource exhaustion and denial of service

The declaration surface exposes capacity, alignment, and explicit
capacity/depth errors; enforcement is reviewed in `rb_tree.cpp`.

### 12. Portability and build behavior

The include guard is project-specific, the only dependency is `<stddef.h>`,
and the header compiled under MSVC C++ and GCC C++ including `-m32`.

### 13. Performance-sensitive structure

The bridge is static typing only and adds no runtime layer or dispatch.

### 14. Tests and evidence

Scalar, portable, targeted AVX-512, Linux ASan, Rust FFI, symbol, and 32-bit
compile results are recorded in `PHASE_05.md`.

### Specialized additions

N/A: no assembly, AVX-512 arithmetic, FFT, or Rust unsafe implementation.

## Findings

`AUD-011` — validated.

## Validation

See `docs/audit/PHASE_05.md`.

## Residual risk and decision

Complete. The internal and public enums are intentionally separate because the
C++ boundary cannot include the public opaque C declarations without type
conflicts. Explicit numeric values and cross-language tests control drift.
