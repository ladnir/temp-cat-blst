# Phase 14 - Scalar and Extension Fields

## Scope

Phase 14 reviews the scalar `Fr` and base `Fp` implementations before moving
upward into extension fields. The first wave reviewed every line of
`blst/scalar/fields/Fr.h`, `Fr.cpp`, `Fp.h`, and `Fp.cpp` at
`c3096482cdab88edc4f05e0c6583c7896fd75bab`, together with only the directly
corresponding C and Rust declarations, wrappers, and tests. CUDA code was
excluded. No build, benchmark, or source edit occurred during the review.

## First-wave findings

The review opened `AUD-067` through `AUD-074`. The approved straightforward
packet covers zero exponentiation, repeated-square boundary safety, public
square-root failure output, secret-boundary cleanup, restoration of retained
scalar declarations, bulk-size hardening, and file documentation/provenance.

`AUD-069` and `AUD-070` are deliberately excluded from that implementation
packet. The current C boundary relies on ABI-compatible but nonidentical C and
C++ function/object types, while several inherited C storage idioms do not
establish valid C++ object lifetime. Their repair requires an explicit
zero-copy representation design and optimized scalar/AVX code-generation
evidence; routine boundary copies are not acceptable.

## Validation policy

The straightforward packet uses only focused Fr/Fp/C-boundary tests,
standalone C/Rust declaration checks, assembly regeneration when a generator
changes, and small optimized code-generation probes. Full suites, sanitizer
sweeps, portable-AVX Release, benchmarks, CUDA, and CI pushes remain separate
permission-gated evidence.

## Straightforward remediation result

Commit `7c0e399` changes the x86 q/x repeated-square generators and all eight
affected ELF, COFF, Mach-O, and Win64 products to preserve full `size_t`
counters. Commit `b6a51bb` implements the remaining approved straightforward
field packet. `AUD-067`, `AUD-068`, `AUD-072`, `AUD-073`, and `AUD-074` are
validated. `AUD-071` is fixed with focused C evidence and remains pending a
retained Rust execution check. `AUD-069` and `AUD-070` remain open.

Validation remained focused. One existing Debug target was rebuilt with one
Ninja worker after the field headers invalidated dependent objects; only three
public C/header tests were executed and passed. A standalone scalar field probe
passed zero-exponent, zero-count, exact-alias, scalar conversion, square-root
failure, and oversized-count cases. The scalar conversion header compiled as
strict C11 under MSVC `/W4 /WX`; GCC 13 and Clang 18 strict syntax checks passed
for the probes and modified C++ boundaries. Assembly regeneration verified all
100 products. MSVC 19.50 Release code generation shows constant-count
`sqr_n_mul` remains one direct native call without the alias slow path, while
constant-zero exponent calls reduce to the required identity copy. No full
test sweep, sanitizer run, benchmark, CUDA job, portable-AVX Release build, CI
push, or Apple run was performed.
