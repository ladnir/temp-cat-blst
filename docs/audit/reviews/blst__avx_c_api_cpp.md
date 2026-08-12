# File Review: `blst/avx_c_api.cpp`

## Record

| Field | Value |
| --- | --- |
| Roadmap phase | Phase 18 expert AVX C boundary |
| Reviewed commit | `e6163cecba9a9f66a8e6ec7aed58cbf2f0ba1672` |
| Date | 2026-08-12 |
| Disposition | in-scope production, complete |
| Profiles | native/portable AVX-512 C ABI; CUDA and frozen BTx/BTE excluded |
| Related findings | `AUD-091`–`AUD-095`, `AUD-101`, `AUD-160`, `AUD-161`, `AUD-229`, `AUD-316`, `AUD-317`, `AUD-322` |

## Purpose and trust boundary

This translation unit converts public opaque AVX storage into selected native
types, checks pointers, sizes, alignment, windows, and scratch capacity, and
publishes C-compatible results. It is an untrusted language boundary: no C++
exception may escape and no unrelated C++ object type may be materialized by a
reinterpret cast merely because layouts match.

## Review results

- Public carrier size/alignment and actual C++ alias assertions remain. The
  former cross-type `AvxGt`/`MillerResult` assertions and cast are removed.
- Final exponentiation invokes the established MatrixNoK implementation on the
  neutral Fp12 data member, preserving exact aliasing without a 1536-byte copy.
- Generic Fp12 multiply/square use the neutral field ring. Square selects
  `fp12_flat_sqr`; GLV is a scalar-decomposition technique and is not applicable
  to one generic Fp12 squaring operation.
- Every fallible entry rejects invalid pointer/shape/capacity inputs before
  publication. Exact alias is retained for same-type arithmetic; partial or
  cross-type overlap is not promised.
- MSM remains variable-time over public scalars. Ordinary and GLV paths require
  subgroup points. Scratch object lifetime/alignment is established at the
  boundary and no hidden allocation or global scheduler exists.
- Friendly batch membership owns challenge derivation and catches callback
  exceptions. Expert relation hooks retain their caller-controlled soundness.
- All six retained G1 GLV scratch/window functions have a local differential
  test and an external caller rationale.

## Evidence and residual risk

The existing portable-AVX Debug target rebuilt successfully. The two focused
new tests pass, including generic non-GT Miller storage, exact aliases, both GLV
executors, and scalar agreement. MatrixNoK structural validation and exact
proof correspondence pass without rebinding the proved GT owner. Native
AVX-512 code generation/runtime, sanitizer/fuzz, Apple, and broad tests remain
milestone evidence. No benchmark was run.
