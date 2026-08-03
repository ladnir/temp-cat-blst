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
