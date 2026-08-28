# File Review: `bindings/rust/src/avx.rs`

## Record

| Field | Value |
| --- | --- |
| Roadmap phase | Phase 18 handwritten Rust AVX surface |
| Reviewed commit | `e6163cecba9a9f66a8e6ec7aed58cbf2f0ba1672` |
| Date | 2026-08-12 |
| Disposition | in-scope production, complete |
| Profiles | Rust with and without the `avx` feature |
| Related findings | `AUD-040`, `AUD-161`, `AUD-318`–`AUD-323` |

## Purpose and trust boundary

This module declares the expert AVX C ABI and provides small safe wrappers for
randomized batch membership and GT MSM. The five aligned public byte arrays are
kept public for raw FFI/source compatibility. They are opaque backend storage,
not serialization, and arbitrary mutation can violate representation and
algebraic preconditions.

## Review results

- The raw declaration inventory exactly matches all 77 C functions and root
  reexports include the newly restored conversions and generic Fp12 operations.
- Carrier sizes/alignment are tested even when AVX execution is unavailable.
- `PartialEq` is stable backend-independent byte-storage equality. Algebraic
  equality remains available through explicit raw AVX functions. `Eq` is
  deliberately not promised for expert-mutable opaque carriers.
- Friendly G1/G2 membership allocation uses `try_reserve_exact`; failure returns
  false and the by-value seed is erased before return.
- `gt_msm` checks scalar shape before FFI, requires prevalidated GT elements,
  asserts that contract only in Debug, and adds no Release membership cost.
- All caller-size-dependent GT-MSM allocations are fallible and map to
  `AllocationFailed`. `MaybeUninit` scratch vectors retain length zero while C
  owns writes, so Rust never claims uninitialized elements are live.
- Empty GT MSM returns the canonical multiplicative identity rather than zero
  Fp12 storage. Nonempty output is published only after successful native
  execution and conversion.
- There is no thread pool, dynamic dispatch, callback allocation, or secret
  scalar API in this module. MSM scalars are explicitly public and variable
  time.

## Evidence and residual risk

`cargo fmt --check` passes. Focused `cargo test --features avx --lib avx::tests`
passes 7/7, including the expected caught Debug precondition panic. A separate
default-feature run passes unavailable-backend equality 1/1. Allocation-failure
injection, Miri, sanitizer, native AVX-512 runtime, Apple, Rustdoc, and the broad
feature matrix remain milestone evidence. No benchmark was run.

## 2026-08-27 secret G1 MSM addendum

`G1MsmWorkspace::secret_msm_into` accepts `&[G1Projective]` and
`&[Scalar]`, so subgroup and reduced-scalar invariants remain typed across the
Rust/C++ boundary. It checks equal lengths and workspace capacity before FFI,
passes the resident Montgomery scalar array without Rust-side byte staging,
and returns `Error` rather than panicking. `for_points` reserves the larger of
the public Pippenger and secret-MSM scratch requirements outside the hot path.
The method performs no allocation, point conversion, or dynamic dispatch.

The focused differential test compares the chunked kernel with independent
constant-time products at empty, tail, exact-chunk, and mixed boundaries up to
65 points. It also checks length mismatch and undersized-workspace rejection.
Peach Criterion measurements are recorded in the C++ owner review; moving the
large tables into this reusable workspace did not regress the selected
thirty-two-point kernel.
