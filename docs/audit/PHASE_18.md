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

## Wave 3: Rust generation, packaging, and high-level contracts

The read-only audit of the Rust root, checked-in scalar FFI artifact, binding
generator, Cargo package, and native build bridge identified `AUD-324` through
`AUD-333`. The implementation landed in reviewable packets:

1. make the scalar raw binding exactly reproducible and keep split AVX, FFT,
   and middle-product declarations in their handwritten modules;
2. make Cargo target/profile/feature/prebuilt selection fail closed, with
   portable AVX explicitly documented as developer emulation; and
3. complete safe-Rust contracts and restore the nonzero `SecretKey` invariant.

This wave also adds lightweight CI policy gates. They must catch declaration
drift, unsafe generator options, invalid feature/target/profile combinations,
stale compiler-support claims, and configuration-warning suppression before a
large native build starts. Native Apple ARM64 remains supported; unsupported
cross builds must supply an explicit toolchain rather than silently compiling
for the host.

Implementation/evidence commits are:

- `f27d7026873ce67fd721b0a0ecb5c7236cd5e88b`: scalar-only bindgen input,
  pinned bindgen 0.65.1 policy, fail-closed transforms, exact 253/253 C/Rust
  declaration ownership, split FFT declarations, and generated pairing
  alignment;
- `835b97d802948fabb1e88383a975df07d98cf0ee`: exact end-to-end regeneration,
  pinned Rust 1.65 formatting and repository configuration, public-name
  allowlists, distinct opaque context tags, system-header exclusion, and a
  path-scoped byte-comparison workflow;
- `b2c93e9d3ec4b4478edca51c9cd67187f859eba4`,
  `7b2ab8caf254e80611d57a1958912b93e7d7bd54`, and
  `906fc2c6964fcff3f8d188b5d373fd25334aa828`: profile/target/feature/toolchain/
  prebuilt routing, checked `blst_std`, no-std Serde, MSRV, deterministic native
  artifacts, and four warning-clean policy tests;
- `778f2f60d1fd4ba19569ba0dd1d4e3f5f1a50b88`: complete reviewed high-level
  Rust contracts, aggregate-flow guidance, non-default/nonzero `SecretKey`,
  private FFI output construction, and focused invariant tests; and
- `6d66fdce136ecdf86d6ba3f3ae0c169f81cfa5b1` and
  `906fc2c6964fcff3f8d188b5d373fd25334aa828`:
  portable AVX removed from obvious Release presets and the ordinary
  CI/deployment matrix, with a separate path-scoped Debug developer workflow;
  explicit optimized emulation remains available.

Focused evidence passes: exact raw declaration parity, Bash/Perl transform
fixtures, 12 generated ABI-layout tests plus the normality test, 9 FFT tests,
13 pairing/uniqueness wrapper tests,
4/4 build-policy tests, explicit-prebuilt linking, no-std Serde with warnings
denied, full-binding rustdoc with warnings denied, and 2/2 secret-key invariant
tests covering all derivation families and malformed scalar decoding. Portable
AVX remains covered by its focused Debug CI lane without being presented as a
deployable backend.

`AUD-324` through `AUD-333` are validated. Bindgen 0.65.1 with Clang/libclang
18.1.3 and the pinned Rust 1.65 formatter regenerated the checked artifact
byte-for-byte from both repository and temporary output locations. The focused
CI lane repeats that comparison only when the scalar header, generator inputs,
artifact, or checks change. No broad native rebuild, benchmark, CUDA, frozen
BTx/BTE audit, native AVX run, or push was performed in closure.
