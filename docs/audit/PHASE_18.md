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

## Wave 4: low-level Rust primitives and XMD input domain

The complete hand review of `bindings/rust/src/primitives.rs` traced all 23
production unsafe blocks through the scalar C ABI. No uninitialized output,
Rust memory-unsafety, exception escape, hidden heap allocation, or compact
Miller-descriptor overread was found. `AUD-334` through `AUD-336` narrow the
Miller contract to subgroup points, document little-endian scalar and
general-Fp12/GT semantics, and make 14 inherent cryptographic results and
decisions must-use without changing ABI or runtime validation.

`AUD-337` closes the SHA-256 input-domain gap found while tracing safe scalar
hashing. Overflow-free checks now cover both the 17-byte oversize-DST prehash
and the full b0 transcript before reading input. The void C boundary clears a
valid output for invalid lengths; safe Rust panics before FFI so invalid input
cannot be confused with a legitimate reduction to zero. The public header was
regenerated exactly into raw Rust.

Implementation commits are `bd0272b` and `dfd0740`. Scalar Debug XMD tests pass
7/7 and scalar Release safe-Rust primitives pass 9/9. Full-binding rustdoc is
warning-clean, formatting passes, bindgen 0.65.1 plus Rust 1.65 rustfmt
regenerates byte-for-byte, and the declaration inventory remains exactly
253/253. Optimized MSVC inlines the public length guard before hashing and emits
no helper call from that boundary. No benchmark, full suite, sanitizer, CUDA,
BTx/BTE audit, portable-AVX Release, native AVX, Apple CI, push, or parallel
heavy work was performed. The scalar Debug target rebuilt serially; frozen
BTx/BTE compiled unchanged only as an existing target dependency.

## Wave 5: Rust support-file closure

The remaining non-benchmark Rust test, example, and publication support files
are complete at implementation commits `16e5b49`, `5f6c423`, and `51b23a5`.
`AUD-338` documents the owned affine carriers' subgroup, validation, mutation,
empty, and allocation responsibilities without adding validation or production
overhead.
`AUD-339` makes publication fail closed on checkout discovery and any tracked,
staged, or untracked source change before creating its required source link. It
requires a real symbolic link, verifies that link resolves to the repository
root, and cleans its isolated target on normal or signaled exit. MSYS directory-
copy emulation is explicitly rejected.

The deterministic Pippenger support macro is correctly classified as test-only.
It now uses direct heap buffers, forms raw pointers for the C exact-alias
contract, and covers non-byte-aligned masking, 256-bit scalars, short buffers,
and over-width rejection. The example states its raw Miller subgroup
precondition, and every reviewed Rust support source has Apache-2.0/SPDX
provenance. Compressed-infinity and conformance behavior required no change.

Focused scalar Release evidence passes G1/G2 multipoint 16/16, compressed
infinity 1/1, scalar conformance 2/2, and the full-binding smoke example. Pinned
bindgen/Rust 1.65 regeneration is exact, C/raw Rust declarations remain 253/253,
and rustdoc plus Rust API gates pass with warnings denied. Shell syntax, dirty
rejection, a native POSIX injected-no-op Cargo cleanup path, and MSYS safe
native-link refusal pass; no registry was contacted. No native rebuild,
benchmark, full suite, sanitizer, CUDA, BTx/BTE audit, portable-AVX Release,
native AVX, Apple CI, push, or live publication ran.
