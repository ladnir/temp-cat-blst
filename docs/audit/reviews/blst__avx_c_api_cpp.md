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

## Phase 16 MSM boundary addendum

Commit `5b08edc8f890ef1cbf4c6443162ec4da985879ee` adds fixed-count,
overflow-safe byte-range validation before every affected MSM output is
published. G2 aggregate scratch sizing now applies the same 8,192-point GLS
cutoff as execution. Focused portable-Debug execution covers the exact cutoff,
all overlap families, partial overlap, failure atomicity, and retained
one-point projective aliasing. The checks execute once at the C boundary; no
bucket loop or decomposition changed. `AUD-389` and `AUD-390` are validated.

## 2026-08-27 secret G1 MSM addendum

`blst_avx_p1s_mult_ct` is a distinct secret-scalar boundary rather than a mode
bit on variable-time Pippenger. It accepts resident AVX projective G1 points
and reduced native Montgomery-form `blst_fr` scalars. Release assumes subgroup
membership; Debug asserts it. The implementation converts scalar field values
to canonical bytes only inside the fixed-time GLV preparation and erases every
secret byte and split component after use.

For eight or more points, the boundary requires aligned caller scratch. It
checks pointer products for overflow, scratch capacity and alignment, and
pairwise disjoint live ranges before arithmetic or output publication. Invalid
input leaves output unchanged. The same scratch is reused sequentially across
32-, 16-, and 8-point chunks, followed by the existing constant-time
single-point tail. No allocation, exception, dynamic dispatch, or backend
conversion is present in the hot path.
