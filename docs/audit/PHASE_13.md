# Phase 13 - Foundational Scalar Trust Chain

## Scope and opening state

Phase 13 begins the file-by-file hand audit with configuration and common
representation policy. The first review packet covers `blst/config.h.in` and
`blst/common/defines.h`. It also records dependent defects in the C ABI and FFT
boundary without claiming complete reviews of those larger files.

The review baseline is commit
`33914f23da23fced86ca51ce9917cf52b59bd96e`. The candidate review of
`defines.h` reconstructed its type, compiler, constant-time, sizing, ABI, and
performance contracts. `config.h.in` was inspected as the configuration input
that fixes those contracts.

## Opening findings

- `AUD-042` (high): `blst_miller_loop_n` can throw across C/Rust FFI because
  it allocates two vectors inside the C function.
- `AUD-043` (medium): empty FFT root preparation can shift a 64-bit value by
  64 after `log2ceil(0)`.
- `AUD-044` (low): `ceil_div` can overflow and divide by zero.
- `AUD-045` (medium evidence gap): the MSVC value barrier lacks optimized
  object-code evidence.
- `AUD-046` (low): compiler warning suppression leaks across translation
  units.
- `AUD-047` (low): unsupported 32-bit limb configuration remains partially
  exposed.

The missing license/provenance header in `defines.h` is tracked as review
remediation rather than as a separate finding.

## Remediation order

1. Inventory the retained C ABI for every possible exception escape without
   editing or building.
2. Remove allocation from `blst_miller_loop_n` while preserving its public
   signature and compact-pointer convention.
3. Harden `log2ceil`, empty FFT root preparation, and `ceil_div` with explicit
   Release results and Debug assertions.
4. Inspect optimized MSVC, Clang, and GCC code for a focused value-barrier
   probe. No barrier change is authorized without evidence and a separate
   review if the current output is unsafe.
5. Remove global warning suppression, restore provenance, and enforce the
   64-bit-only configuration at CMake configure time.
6. Reconcile the final exact-commit review records and manifest rows.

## C ABI exception policy

No C++ exception may cross a CatBLST C ABI boundary. Legacy `void` entry points
must be allocation-free and nonthrowing for inputs satisfying their documented
preconditions. Fallible entry points must translate failure to an explicit
status, `bool`, or null result and must not partially update caller output.
Adding `noexcept` is not a substitute for proving this contract because an
escaping exception would terminate the process.

## Validation and resource policy

Use the smallest directly affected build or test target for each remediation.
Portable-AVX Release, broad sanitizer or fuzz sweeps, full C/C++/Rust suites,
benchmarks, CUDA, and frozen BTx/BTE are excluded from this packet unless the
user explicitly approves them. Compiler probes are compile-only and are run
sequentially. Performance-sensitive changes require structural or generated-
code evidence; no heap allocation, runtime dispatch, or threading may be added.

## Current status

The findings and sequence are recorded. Implementation and final validation
remain pending. Neither `config.h.in` nor `defines.h` is complete until its
review record names the final exact commit and the manifest passes integrity
checks.
