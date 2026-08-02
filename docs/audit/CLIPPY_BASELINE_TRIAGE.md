# Clippy Baseline Triage

## Scope and configuration

The first report-only Clippy baseline ran at commit `b3a537359ffa` on Windows
x86-64 with Rust/Cargo/Clippy 1.94.1. It checked the Rust library and build
script in Debug profile with Clippy's `all`, `pedantic`, and `nursery` groups.
Tests, examples, benchmarks, frozen BTx/BTE, experimental CUDA, and unsupported
platform targets were not linted. The native library built by each Cargo
profile retained `BLST_TEST=OFF`, `BLST_BENCHMARK=OFF`, `BLST_CUDA=OFF`, and
`BLST_PARALLEL_MSM=OFF`.

The feature profiles were:

| Profile | Cargo features | Native configuration |
| --- | --- | --- |
| default | none | supported host scalar/assembly Debug |
| full | `full-rust-bindings` | supported host scalar/assembly Debug |
| portable | `portable` | portable scalar Debug |
| avx | `avx` | portable AVX Debug on Windows; not native AVX evidence |
| full-serde | `full-rust-bindings,serde-secret` | supported host scalar/assembly Debug |

Portable AVX Release was not configured or compiled. The reproducible runner
is `scripts/audit/run_clippy_local.ps1`. The complete local evidence run is
`b3a537359ffa-20260802T235710Z`; all five cached profiles completed in 37
seconds.

## Results

Diagnostics produced twice by the two signature-orientation macro expansions
were deduplicated by lint, source file, line, and column.

| Profiles | Raw warnings | Canonical warnings |
| --- | ---: | ---: |
| default, portable, avx | 63 each | 63 each |
| full, full-serde | 169 each | 133 each |

Only three canonical diagnostics belong to `clippy::all`:

| Lint | Disposition |
| --- | --- |
| `clippy::uninit_vec` | Proposed `AUD-040`. `gt_msm` marks uninitialized `u8` capacity as initialized before passing it to C++, violating `Vec::set_len`'s safety contract even though the C++ implementation zeroes the whole flag buffer before reading it. |
| `clippy::manual_div_ceil` | Included in `AUD-040`. The current add-before-divide expression can overflow. Do not adopt `usize::div_ceil` until the crate has a declared compatible MSRV; use an overflow-free quotient/remainder calculation instead. |
| `clippy::manual_is_multiple_of` | Declined. The remainder check is exact and obvious, while `is_multiple_of` would raise the crate's compiler requirement substantially without changing behavior. |

The remaining 60 default-surface and 130 full-surface diagnostics come from
the pedantic/nursery groups:

- missing `# Errors`/`# Panics`, `must_use`, wildcard imports, `Self`, format,
  const-function, and control-flow suggestions are documentation or style
  debt, not defects;
- public underscore-prefixed fields are generated ABI padding and must retain
  their names and layout;
- raw-pointer borrow/cast suggestions preserve the same FFI operation and add
  no safety proof;
- the mutable receiver on `Pairing::ctx` is intentional because returning a
  mutable context pointer from `&self` would weaken the aliasing contract; and
- eight canonical reference/slice transmutes rely on explicit transparent
  wrapper layout. They require a deliberate unsafe-boundary review in Phase
  18 rather than mechanical lint rewrites.

No automatic fix was applied. Documentation/style cleanup can be batched with
the owning file reviews so public documentation and MSRV decisions are made
coherently.

## Proposed finding

### AUD-040 - Rust AVX MSM scratch storage violates `Vec` initialization contract

`gt_msm` reserves a `Vec<u8>` and immediately uses unsafe `set_len` without
initializing the new elements. The foreign implementation currently fills the
whole `bucket_used` region before reading it, but that does not satisfy Rust's
precondition at the earlier `set_len` call. The adjacent scratch vector uses
`MaybeUninit` safely, although its `set_len` is unnecessary. Its manual ceiling
division can also overflow before allocation.

Proposed remediation, pending approval:

- keep both work buffers allocation-only and represent foreign-written bytes
  as `MaybeUninit` storage;
- remove both `set_len` calls and pass pointers into allocated capacity plus
  the exact required byte counts, so no redundant zeroing is introduced;
- compute the scratch element count with overflow-free quotient/remainder
  arithmetic compatible with the current compiler floor;
- retain the C++ full-buffer initialization and size/alignment checks;
- add a focused allocation-plan boundary test and run the existing AVX MSM
  error/identity tests under the `avx` feature; and
- inspect optimized Rust assembly if the storage rewrite does not trivially
  disappear, rejecting redundant initialization or additional hot-path work.

## Residual limits

Clippy does not validate native FFI implementations, allocation feasibility,
constant-time behavior, foreign object lifetime, or wrapper layout claims.
Phase 18 must manually review each unsafe block and transparent-layout
conversion. Native AVX-512 execution remains Phase 16 evidence; this Windows
`avx` profile uses the Debug emulation backend.
