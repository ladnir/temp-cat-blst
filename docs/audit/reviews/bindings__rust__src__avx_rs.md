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

## 2026-08-27 public native-scalar MSM addendum

`G1MsmWorkspace::msm_into` now accepts `&[Scalar]`, eliminating BTX's
Rust-side canonical byte vector. The FFI receives the resident Montgomery
scalar array. The method checks dimensions and reusable capacity before the
call and allocates nothing. The selected C++ implementation is variable-time
and therefore retains the explicit public-scalar contract. The former
canonical-byte Pippenger operation remains available as
`pippenger_msm_into` for controlled comparison and low-level use.

The differential test spans every fixed chunk boundary, the 192/193 algorithm
boundary, and the larger Pippenger path through 257 points. It also checks
length mismatch and undersized-workspace rejection. Native Peach benchmark
results and C malformed-input coverage are recorded in the C++ owner reviews.
The downstream BTX benchmark uses Lagrange coefficients for nonconsecutive
validator indices rather than sparse integer scalars; its eight-point AVX path
measures 60.032 us versus 404.37 us for the scalar backend on Peach.

## 2026-09-03 prepared-pairing range addendum

`PreparedG2Batch::miller_loop_product_indices` accepts equally sized native G1
and prepared-index slices. It reports shape and index errors before FFI and
allocates no memory. The immutable prepared batch is `Send + Sync`; the new
operation needs no mutable workspace.

`MillerLoopResult` is a distinct transparent owner for a pre-final Miller
product. It exposes only identity construction, multiplication, and final
exponentiation. The typed interface therefore does not represent a Miller
product as `Gt`, even though the C ABI uses compatible opaque Fp12 storage.

Peach tests compare direct, complete, partitioned, reordered, and concurrent
evaluation. Native AVX unit and conformance tests pass. The C adapter also
passes its ASan and UBSan corpus replay. The benchmark records 1.70% and 1.74%
serial partition overhead at 128 and 256 pairs for the provisional ten-pair
grain. Four-worker kernel scaling reaches 3.04 and 3.42 times the full prepared
serial path at those anchors. Detailed evidence is in
`docs/audit/BTX_PREPARED_PAIRING_RANGES_2026-09-03.md`.

## 2026-09-04 alpha-pairing range addendum

`PreparedG2BatchWorkspace::with_capacity` bounds native scratch by the largest
range assigned to one worker rather than by the complete prepared table.
`pairing_range_uninit_into` accepts contiguous prepared indices and
`MaybeUninit<Gt>` output. The wrapper validates the range and workspace
capacity before FFI, allocates no memory during evaluation, and treats the
output as initialized only after the C operation succeeds. The ordinary
initialized-output method delegates to this same boundary.

The immutable prepared owner may be shared between workers. Each concurrent
call requires a distinct mutable workspace and a disjoint output range. The
wrapper performs no scalar/AVX conversion: successful output is already the
native AVX GT representation. Full ownership, failure, and performance notes
are in `docs/audit/BTX_ALPHA_PAIRING_RANGES_2026-09-04.md`. The focused native
AVX compatibility test passes under Clang 21 ASan/UBSan on Peach.
