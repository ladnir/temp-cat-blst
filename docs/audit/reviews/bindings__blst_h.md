# File Review: `bindings/blst.h`

## Record

| Field | Value |
| --- | --- |
| Path | `bindings/blst.h` |
| Roadmap phase | Phase 18 public C/C++/Rust surfaces |
| Reviewed implementation commit | `98faac3fb001f59f72fd1b131c278c101953f1a0` |
| Reviewer | Codex audit team |
| Date | 2026-08-06 |
| Disposition | in-scope production |
| Profiles | C99-or-later, native C++, scalar and AVX declarations |
| Related files | `bindings/blst_aux.h`, all split binding headers, `test/BindingHeader_tests.cpp` |

This record reviews the exact implementation commit above.

## Purpose and compatibility contract

This is the canonical upstream-compatible include name. It intentionally owns
no types, declarations, definitions, or configuration logic; it delegates the
complete supported C API surface to `blst_aux.h`. Existing BLST source that
includes `blst.h` must continue to receive the upstream declarations, while
CatBLST additions exposed by the compatibility umbrella must agree with the
split headers.

The include must be repeatable, standalone, valid in C99 or later and C++, and
must not alter ABI layout or generate a runtime dependency. All mathematical,
validation, aliasing, timing, failure, and ownership contracts belong to the
split declaration that owns each operation.

## Complete review

1. All 21 original lines and the current one-line documentation addition were
   inspected.
2. The include guard `BLST_H` is project-prefixed and not reserved.
3. The file has no definitions, writable state, allocation, exceptions,
   dispatch, threading, arithmetic, pointer operations, or secret handling.
4. `blst.h -> blst_aux.h` is acyclic and every nested header is guarded.
5. C and C++ type selection remains owned by `blst_types.h`; the umbrella adds
   no alternate representation or ABI decision.
6. Frozen BTx/BTE headers continue to include this name. Exposing an additional
   declaration header changes no frozen implementation or symbol.
7. Header history was traced through the local split (`c71496b`), compatibility
   restoration (`d3cb4d5`), and non-reserved guard cleanup (`a55bc3a`).
8. `AUD-185` was the only defect: FFT was absent from the transitive surface
   even though the compatibility documentation promised all split APIs.

## Validation

- A strict MSVC C11 `/W4 /WX` consumer included only `bindings/blst.h`, linked
  the existing scalar library, and successfully executed an n=1
  `blst_fr_fft`.
- A strict MSVC C++20 `/permissive- /W4 /WX` consumer did the same while also
  checking the scalar standard-layout and AVX opaque-shell size contracts.
- `BindingHeader_tests.cpp` now includes no direct FFT header and fails at
  preprocessing if `BLST_FFT_H` is not supplied through `blst.h`; the complete
  translation unit parses under strict MSVC warnings.
- No broad build, benchmark, sanitizer, CUDA, BTx/BTE audit, Apple CI, or
  portable-AVX Release build ran.

## Residual risk

Runtime and semantic assurance remains with the individually reviewed split
headers and implementations. Phase 18 still requires packaging/install
consumers and supported non-MSVC compiler evidence before its global exit gate.
