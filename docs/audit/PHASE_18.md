# Phase 18: C, Native C++, and Rust Surfaces

## Scope and state

Phase 18 audits the supported C99-or-later ABI, the main native C++ API, and
the safe/generated Rust surfaces. It verifies representation, ownership,
pointer extent, alignment, aliasing, object lifetime, validation domains,
failure publication, exception containment, timing, language agreement, and
upstream C/source compatibility. Every retained C entry remains subject to the
program-wide rule that no exception may cross the ABI.

The phase is in progress. Declarations are closed only after their owning
implementation has been reviewed; implementation findings do not make the
entire language surface complete. BTx/BTE is frozen and CUDA is experimental.

## Wave 1: scalar public header and boundary

At commits `02f5e5ff078b251ca3b5f86f4f83ea64af664bc7` and
`aa4c7ed187f5a63496f0e94f9284108cbd46c7f1`, Wave 1 completed:

- the 253-declaration `bindings/blst_scalar.h` inventory across C, selected
  native C++, implementations, and raw/safe Rust consumers;
- the public representation and shared result types in `blst_types.h` and
  `error_codes.h`;
- scalar C-boundary implementation review in `c_api.cpp` and `exports.cpp`;
- selected field facades and focused pairing, uniqueness, MSM, and conversion
  lifetime boundaries.

Findings `AUD-186` through `AUD-197` cover caller scratch lifetime, distinct
opaque tags, scalar canonicality wording, C/Rust result ABI, wire-buffer
extent, native-Fr MSM contracts, restored upstream C++ defaults, nodiscard,
empty scratch sizing, experimental-CUDA header cleanup, a reserved guard, and
opaque context lifetime. `AUD-070` received a linked exact-alias follow-up.

## Focused evidence

- Scalar Release rebuilt serially; focused MSM tests pass 12/12.
- Binding-header tests pass 17/17; strict MSVC C11 and C++20 `/W4 /WX /Zs`
  standalone probes pass.
- Focused pairing/aggregate tests pass 9/9 and uniqueness tests pass 3/3.
- The standalone field-boundary Release executable passes.
- Targeted Rust layout/unknown-value tests for both result newtypes pass and
  the full bindings compile.
- Optimized Pippenger and opaque-context inspection finds no allocation,
  lifetime helper call, construction loop, larger frame, or tile-wrapper size
  change. No benchmark was warranted.

No full suite, sanitizer sweep, CUDA, frozen BTx/BTE audit, Apple CI, native
AVX-512 run, or portable-AVX Release build was performed. Negative
`-fshort-enums`, installed-SDK, cross-compiler, and full cross-language package
evidence remains for the Phase 18 release gate.
