# Phase 8 Caller-Threaded Rust Runtime Evidence

## Record

- Date: 2026-07-31
- Implementation commit: `0d655ff`
- State: validated and committed
- Scope: Rust internal scheduling, upstream `07b109d`, and `AUD-021`
- Benchmarks: not run or audited
- CUDA: experimental and not exercised
- BTx/BTE: frozen and untouched

## Decision

The Rust crate is a bottom-of-stack cryptographic component and must not own
process-wide scheduling policy. Every Rust API call now executes entirely on
its invoking thread. The crate does not create threads, discover a CPU count,
own an executor, or queue background work.

This supersedes the earlier plan to translate upstream `07b109d`. Bounded
channels would improve allocations inside the wrong architecture; removing all
internal scheduling removes the channels and their deadlock/lifetime surface.
The upstream row is therefore rejected as obsolete under the stronger local
runtime contract.

## Implementation

- Removed the global `threadpool::ThreadPool`, leaked `static mut` pointer,
  `Once`/`Mutex` initialization, and `ThreadPoolExt`.
- Removed the unsafe conversion of borrowed jobs to `'static` closures.
- Removed all production channels, barriers, worker counters, validity flags,
  CPU-count planning, and Rust-side parallel Pippenger tiling.
- Removed the normal `threadpool` dependency and obsolete `no-threads`
  feature/build-script configuration.
- Made the existing sequential multi-Miller-loop, aggregate verification,
  weighted aggregate verification, affine conversion, MSM, addition, and
  validation implementations universal across host and freestanding builds.
- Consolidated the duplicate std/no-std Pippenger sources into
  `bindings/rust/src/pippenger.rs`.
- Preserved checked scalar-buffer arithmetic, `MaybeUninit` output/scratch
  handling, and the optimized native Pippenger kernel.
- Documented caller-owned scheduling and added a repository policy check that
  rejects internal schedulers and production thread creation.

Public method signatures are unchanged. Scheduling behavior deliberately
changes from implicit multi-core execution to synchronous caller-thread
execution. Validation errors are now deterministic and input-ordered; notably,
MSM validation reports `BLST_PK_IS_INFINITY` instead of collapsing all worker
failures to `BLST_POINT_NOT_IN_GROUP`.

## Caller-managed parallelism

- Pairing callers can use one owned `Pairing` context per shard and reduce with
  `Pairing::merge`.
- Multi-Miller-loop callers can process paired slices independently and
  multiply the partial Fp12 results before final exponentiation.
- MSM callers can split matching point/scalar ranges and add partial group
  results using their own runtime.

No executor trait, callback interface, Rayon dependency, or replacement local
pool was added. A scheduler-neutral shard/plan API should be considered only
if later single-process measurements demonstrate that caller-side composition
cannot meet an actual throughput requirement.

## Focused validation

Following `VALIDATION_POLICY.md`, validation was limited to the affected Rust
surface:

| Surface | Evidence | Result |
| --- | --- | --- |
| Default Rust crate | `cargo test --lib` | 38/38 passed |
| Full Rust binding | `cargo test --lib --features full-rust-bindings` | 74/74 passed |
| Portable scalar/full binding | `cargo test --lib --features "full-rust-bindings portable"` | 74/74 passed |
| Freestanding configuration | `BLST_TEST_NO_STD=1 cargo check --features "full-rust-bindings portable"` | Passed |
| Compile-only surfaces | Default and full-binding `cargo check` | Passed without warnings |
| Dependency graph | Normal Cargo dependency tree | No thread pool, CPU-count, Rayon, or Crossbeam dependency |
| Structural policy | `check_rust_threading_policy.ps1` | Passed; the only `std::thread` use is an explicitly allowed caller-thread test |
| Formatting/integrity | Cargo formatting and Git whitespace | Passed |

The deterministic default-surface multi-Miller-loop regression compares the
batch result with an ordered product and covers the empty identity. Existing
G1/G2 tests cover 2,000-point MSM, affine conversion, addition, empty/zero-bit
work, aggregate verification, and weighted verification. New G1/G2 validation
checks prove exact infinity error reporting.

No C/C++, AVX-512 implementation, assembly, CUDA, or benchmark execution was
warranted because the native arithmetic implementation did not change.

## Residual considerations

- A single large Rust operation no longer automatically consumes all logical
  CPUs. This is intentional resource-policy behavior, not a performance bug in
  the arithmetic kernel.
- Caller-managed parallel throughput has not been benchmarked because
  benchmarks are outside this audit packet. If evaluated later, run only one
  benchmark process at a time.
- At Phase 8 close, WASI execution remained in a separately approved CI packet.
  The later Phase 12 scope decision rejected that platform commitment; the
  forced no-std build remains evidence that this change does not depend on host
  threading support.
