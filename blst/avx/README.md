# AVX-512 Backend

This directory contains category-blst's production x86_64 AVX-512/IFMA
backend. The selected implementation uses fixed BLS12-381 MatrixNoK arithmetic
derived from [VROOM](https://github.com/SimonLangowski/VROOM) and adapted to
category-blst's APIs, contracts, tests, and proof interfaces.

## Build and deployment

`BLST_AVX512=ON` compiles the native backend and expert AVX API.
`BLST_DEFAULT_AVX=ON` also selects AVX for the library's default
implementation. Selection occurs at compile time; the library does not perform
runtime dispatch between scalar and AVX implementations.

Native binaries require an x86_64 CPU and operating system with AVX2,
AVX-512F, AVX-512DQ, and AVX-512IFMA support. The supported native compilers
are Clang 21 or newer and GCC 15 or newer. Running a native binary on an
incompatible system may fault.

The portable vector emulation exists only for developer Debug and differential
testing. It is not a supported deployment backend or performance evidence.

## Source layout

- [`cpu/vector/`](cpu/vector/) provides the selected vector abstraction.
- [`rns/`](rns/) provides fixed MatrixNoK conversion, reduction, and bounded
  carriers.
- [`fields/`](fields/) provides Fp2/Fp12 arithmetic, inversion, and final
  exponentiation.
- [`curve/`](curve/) provides G1/G2 points, GLV/GLS multiplication, maps, and
  multi-scalar multiplication.
- [`pairing/`](pairing/) provides the Miller-loop and pairing owners.
- [`gt/`](gt/) provides GT scalar decomposition and exponentiation.
- [`fft/`](fft/) provides AVX GT FFT integration.
- [`generated_matrixnok_constants.hpp`](generated_matrixnok_constants.hpp) and
  [`generated_matrixnok_curve_constants.hpp`](generated_matrixnok_curve_constants.hpp)
  contain the fixed generated constants.

Public contracts are defined by the C, C++, and Rust declarations. Internal
MatrixNoK templates have narrow representation and bound preconditions; they
do not form a separate supported API and do not add runtime validation to hot
arithmetic.

## Constants, provenance, and proof scope

Production uses checked-in compile-time constants. Their initialization
performs no runtime parameter search, arbitrary-precision computation,
allocation, or mutable global setup. The standard-library validators
independently reconstruct the [RNS and field constants](../../scripts/audit/validate_matrixnok_constants.py)
and [curve constants](../../scripts/audit/validate_matrixnok_curve_constants.py).

The retained VROOM source pins and adaptation decisions are recorded in the
[upstream manifest](../../docs/avx/VROOM_UPSTREAM_MANIFEST.md) and
[port map](../../docs/avx/MATRIXNOK_PORT_MAP.md). Exact machine-checked claims
are bound to source and constant hashes in the
[correspondence ledger](../../proof/avx/CORRESPONDENCE.csv), with the claim
boundaries described in the [proof plan](../../docs/avx/PROOF_PLAN.md).

Those records do not prove the entire AVX backend or public API. Compiler
lowering, ISA behavior, microarchitectural constant-time behavior, generic FFT
orchestration, C ABI state machines, scalar observation bridges, and the
concrete reciprocal provider used by Fp12 inversion remain outside the Lean
claims. Some subgroup-identification facts also rely on the external
mathematics cited by the proof plan. A change to a mapped constant, bound,
representation, or operation schedule reopens the corresponding proof record.
