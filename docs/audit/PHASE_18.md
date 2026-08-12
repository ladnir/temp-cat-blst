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

## Wave 2: expert AVX C and Rust surfaces

The expert AVX header, C boundary, and handwritten Rust AVX module received a
complete source/declaration/caller pass. `AUD-316` through `AUD-323` record the
result. The packet removes a C++ object-identity violation at final
exponentiation, separates generic Fp12 multiplication/squaring from the GT
facade, diagnoses the narrow GT-MSM contract in Debug, makes Rust allocation
fallible, repairs unavailable-backend equality and empty-MSM identity, and
brings the raw Rust AVX declarations to exact 77/77 parity with C.

The six G1 GLV scratch/window hooks remain supported: BTX-DKG production uses
the automatic split-scratch family and its tuning benchmarks use the explicit
window family. Local focused coverage now exercises all six without changing
or auditing frozen BTx/BTE or the external repository.

Portable-AVX Debug rebuilt the focused test executable and passed the two new
C++ regressions. Rust AVX tests pass 7/7, and the unavailable-backend equality
test passes separately without the `avx` feature. Strict MSVC C11/C++20 header
parses pass, the MatrixNoK dispatch gate passes, and all 21 accepted proof rows
and 61 mapped hashes remain valid. No benchmark, CUDA, portable-AVX Release,
full suite, sanitizer, native AVX-512, remote, Apple, or push was performed.
The implementation and exact review state are recorded at `e6163ce`.
