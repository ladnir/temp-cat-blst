# Phase 2 Rust Hardening Evidence

## Record

- Date: 2026-07-30
- Baseline commit: `0d79438`
- Implementation commits: `5951b81` and `ddaeaf9`
- State: validated, independently reviewed, and committed
- Scope: `AUD-001`, `AUD-002`, `AUD-003`, and review-discovered `AUD-014`
  and `AUD-015`
- Upstream decisions closed: `1172c18`, `96a46e6`, and `6271918`
- Benchmarks: disabled and not run
- CUDA: disabled

This record closes findings and upstream decisions, not whole-file reviews.
The affected production files remain pending in `FILE_MANIFEST.csv` until
their one-file-at-a-time checklist records are complete.

## Implemented changes

### AUD-001: pairing DST ownership

Rust `Pairing` now owns one aligned allocation containing the exact native
context followed immediately by its copied DST. Exact trailing placement
engages the native self-relative DST representation, so moving the Rust value
between threads does not invalidate the native pointer.

Reinitialization reuses adequate trailing capacity. Growth allocates and
initializes a replacement before swapping and freeing the old allocation.
`Drop` retains the exact allocation layout. The type remains non-`Clone` and
non-`Sync`; its `Send` rationale is documented. The public C header now states
the ordinary borrowed-DST contract and the special trailing-allocation rule.

### AUD-002: weighted-pairing scalar length

`mul_n_aggregate` computes `ceil(nbits / 8)` without addition overflow and
returns `BLST_BAD_SCALAR` before native dispatch when the scalar slice is too
short. Returning the existing error type is safer for this API than upstream's
panic and lets high-level callers propagate the failure normally.

### AUD-003: Pippenger planner

The threaded MSM implementation routes work directly whenever the machine is
oversubscribed relative to the point count. Planner inputs and every
CPU-derived window subtraction are clamped to a nonzero width. Tests inject a
reported 1,024 CPUs for both G1 and G2 and exercise large-bit planner inputs
that previously caused subtraction underflow.

Worker results use uniquely written `UnsafeCell<MaybeUninit<T: Copy>>` slots
aligned to 64 bytes. This removes the small-work channel allocation and avoids
false sharing while preserving direct, statically dispatched native calls.

### AUD-014: scalar layout and empty work

Independent review found that `(nbits + 7) / 8` and
`bytes_per_scalar * npoints` could overflow before native MSM calls. Both std
and no-std paths now use checked arithmetic and slice-derived pointers. Empty
point sets, empty addition, and zero-bit multiplication return group identity
without indexing empty slices.

### AUD-015: native output and scratch storage

Independent review also found initialized `Vec<T>` lengths being installed
over uninitialized affine, scratch, grid, and test storage. FFI output and
scratch allocations now remain `MaybeUninit` until their native initialization
contract is satisfied. Fully initialized affine allocations are reconstructed
with their original pointer, length, and capacity. Concurrent result slots
have one writer and are read only after acquire synchronization.

## Validation matrix

| Profile | Configuration | Result |
| --- | --- | --- |
| Rust native | `full-rust-bindings` | Pairing 7/7; planner 2/2; oversubscription 2/2; empty 2/2; zero-bit 2/2 |
| Rust portable | `full-rust-bindings portable` | Pairing 7/7; oversubscription 2/2; empty 2/2; zero-bit 2/2 |
| Rust no threads | `full-rust-bindings no-threads` | `cargo check` passed |
| Rust no-std portable | `BLST_TEST_NO_STD=1`, `full-rust-bindings portable` | Test build passed; layout 1/1; empty 2/2; zero-bit 2/2 |
| Static checks | all Rust sources | `cargo fmt --check`, configuration checks, and `git diff --check` passed |

The implementation agents and an independent review lane examined ownership,
alignment, initialization validity, synchronization, planner arithmetic,
std/no-std parity, panic behavior, and hot-path structure.

## Residual validation

- Native FFI prevents useful Miri coverage; add sanitizer coverage in the
  Linux Clang release-validation matrix.
- Injected worker-panic tests are not yet present.
- Existing affine-conversion and validation barrier paths could wait if a
  worker or native C call unwound contrary to their contract. Ordinary
  input-driven Rust indexing and sizing failures were removed.
- The full slow no-thread Rust suite exceeded the bounded interactive window;
  focused tests and native, portable, no-thread, and no-std compilation passed.
- Remaining planner arithmetic is bounded by production-derived windows,
  actual pool counts, and realizable slice lengths; a later full arithmetic
  review can convert it to checked helper functions for defense in depth.
