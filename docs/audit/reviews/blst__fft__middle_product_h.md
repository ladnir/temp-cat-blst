# File Review: `blst/fft/middle_product.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/fft/middle_product.h` |
| Roadmap phase / queue | Phase 17 Packet 17J |
| Reviewed implementation commit | `e54bf0a4b4f84a5b5b6d50a1326414d5a1c85713` |
| Reviewer | Codex audit team |
| Date | 2026-08-28 |
| Disposition | `in_scope_production` |
| Subsystem | Pairing middle product, FFT and direct multipair algorithms |
| Profiles | scalar, selected native AVX-512, and portable-AVX Debug emulation |
| Generation/provenance | Locally composed C++ using audited BLST pairing/FFT and MatrixNoK AVX primitives; Apache-2.0 and Category Labs |
| Upstream references | `not_applicable: CatBLST-specific middle-product owner` |
| Related files | `blst/fft/ctFft.h`, `blst/fft/context.h`, scalar/AVX point and GT owners, `blst/middle_product_c_api.cpp`, `bindings/blst_middle_product.h`, Rust wrappers, middle-product tests |

## Purpose and mathematical contract

For `a[0..n)` in G1 and `b[0..2n)` in G2, the owner computes

`out[k] = product_{i=0}^{n-1} e(a[i], b[n+k-i])`, for `0 <= k < n`.

An independent Fr oracle evaluates the same diagonal scalar products and lifts
them through pairing bilinearity. The basic full-convolution implementation,
truncated FFT implementation, scalar multipair implementation, and AVX
multipair implementation agree on retained sizes and infinity distributions.
The `b[0]` slot is part of the established `2n` representation but does not
enter this selected middle window.

G1, G2, and GT denote prime-order subgroups. Inputs must already be reduced,
on-curve subgroup points; infinity is the identity. No Release membership check
is performed. Points and lengths are public, and every algorithm is variable-
time in those inputs.

## Ownership, publication, and concurrency

- FFT and multipair precomputes own all backend-specific immutable state.
  Their fields remain expert-visible, but mutation after construction violates
  the contract.
- One immutable precompute may be shared across threads. Each evaluation needs
  a distinct mutable workspace and output.
- Selected-backend entry points use the original G1/G2 spans directly.
  Explicit non-selected scalar/AVX entry points allocate conversion vectors at
  the backend boundary.
- Workspace growth, conversion, and sparse batch reservation complete before
  output publication. Allocation or validation failure leaves output unchanged.
- Exact output aliasing is supported only where explicitly documented for
  `combine`; other input/output overlap is outside the contract.
- No global thread pool or mutable process-wide state is used. AVX prepared-root
  selection is scoped through a nestable thread-local pointer.

## Correctness and safety review

### Sizing and indexing

The signed-index FFT bound is `2^30`, yielding `n <= 2^29` for optimized middle
products. Zero, one-past-bound, and hostile `size_t` lengths fail before span
construction or allocation. All `2n`, `3n-1`, admissible rounding, diagonal,
fold, and extraction indices were checked against that bound. On the supported
64-bit platform, the largest fold product is below `2^59`.

### FFT path

The immutable B precompute stores ordered forward/inverse roots, one reusable
forward plan, one middle-range inverse plan, the transformed folded B operand,
and the tail used by `update`. Evaluation fills and transforms A, pairs each
frequency point, applies the inverse GT transform, and publishes exactly the
middle `n` values. Root direction, scale, plan shapes, and precompute dimensions
are checked before workspace growth or writes.

### Multipair path

Each nonidentity G2 point owns 68 prepared Miller lines; identity slots are
marked and never read. Nonidentity G1 points are collected with parallel source
indices. Dense and sparse scalar/AVX diagonal indexing stays within `[1,2n)`.
Miller results are final-exponentiated before becoming GT output. AVX batches
retain all six final-exponentiation scratch buffers in reusable workspaces.

`AUD-311` found that scalar sparse line/scale vectors could grow during a later
diagonal after earlier outputs were published. Both now reserve the known
nonidentity-G1 maximum before the output loop, matching AVX ownership and
removing allocation from the scalar sparse inner loop.

### Selected-backend dispatch

`AUD-312` found unconditional point conversion vectors even when public and
backend point types were identical. `if constexpr` now supplies original spans
for the selected scalar or AVX context; conversion helpers exist only in the
compile-time non-selected branch. No runtime condition, callback, virtual call,
type erasure, or additional copy was introduced.

### Exceptions and API behavior

Malformed native shapes throw `std::invalid_argument`; unsupported lengths
throw `std::length_error`; allocation can throw `std::bad_alloc`. Valid
arithmetic is nonthrowing. Invalid `MiddleProductAlgorithm` enum values are a
narrow-contract violation, diagnosed in Debug and performing no Release work.
The C boundary translates every exception and stages evaluation/update outputs.

## Performance-sensitive structure

FFT plans, inverse roots, transformed B values, prepared Miller lines, and AVX
final-exponentiation scratch are reused. Auto dispatch remains backend-specific,
constexpr, and based on the existing Peach measurements across FFT transform
buckets. This packet changes no arithmetic, cutoff, runtime dispatch, or
benchmark result. The direct selected spans remove one full point copy and heap
allocation; sparse reservation moves potential growth before the loop and
becomes a capacity comparison after warm-up.

## Findings

`AUD-088`, `AUD-277`, `AUD-282` through `AUD-291`, `AUD-303`, `AUD-311`,
`AUD-312`, and `AUD-404`.

## Validation

- MSVC 19.50 scalar-default plus portable-AVX Debug: 22/22 focused sizing,
  precompute, scalar/AVX multipair, dispatch, and C-boundary tests passed.
- MSVC 19.50 scalar-only Release: 14/14 focused tests passed.
- MSVC 19.50 AVX-default portable Debug: two selected middle-product and
  update/combine tests passed.
- Sparse scalar workspace capacity and both C `combine` exact-alias directions
  have direct regressions.
- Structural scan places all four conversion calls in compile-time non-selected
  branches and both scalar sparse reserves before publication.
- No new benchmark ran because arithmetic and dispatch are unchanged and the
  changes strictly remove/move allocations and copies.
- No CUDA, portable-AVX Release, full suite, sanitizer, native AVX-512 hardware,
  remote, Apple, or BTx/BTE audit/change was performed.

## Residual risk and decision

Auto cutoffs remain measurements for the Peach deployment CPU rather than a
universal microarchitecture guarantee. Native AVX-512, sanitizer, fuzz, Apple,
and broader release evidence remain milestone gates. Public expert fields and
invalid enum casts retain narrow native contracts. The containing file is
complete for line-by-line hand review.

## Packet 17K delta

Commit `ec18911a01e097bcbdacd9f16634049b4807c998` replaces the AVX inverse-root
call `set_scale_from_size(roots.size())` with `prepare_inverse_scale()`. The
prepared owner now derives the identical public `1/n` value from its owned root
count, eliminating the possibility of caller divergence without changing FFT
arithmetic, dispatch, thresholds, allocation publication, or output behavior.
The Packet 17K four-test scalar/AVX differential and invariant set passes, and
the containing-file completion remains valid at this delta commit.

## Native AVX pointwise-pairing adoption addendum

Implementation commit `713a012583657102c29a4287783585f822baecc7`
changed the native AVX FFT middle product in two connected ways. It batch
normalizes the transformed G1 operand once per evaluation, and it batch
normalizes and prepares the transformed G2 operand once per immutable
precompute. Each pointwise product then evaluates one affine G1 point against
the corresponding prepared G2 line table. Scalar middle product retains the
direct projective pairing path.

Policy commit `e54bf0a4b4f84a5b5b6d50a1326414d5a1c85713` makes the
AVX choice an explicit compile-time production policy rather than a runtime
option. The generic implementation still contains both branches because its
scalar specialization selects the direct branch; the AVX context selects only
prepared batch-affine pointwise pairing. A compile-time assertion guards that
policy. There is no dynamic dispatch, per-call algorithm flag, or additional
hot-loop allocation.

### General-purpose A/B evidence

The A/B used exact source commit `0fb92c8f`, except that the control build
changed only `AvxMiddleProductContext::batch_affine_pointwise` from `true` to
`false`. Both builds used Clang 21.1.8 in Release mode on Peach's Ryzen 9 7950X
and were pinned to CPU 2. Criterion used ten samples after a 500 ms warm-up.
Three complete runs were serialized; the third reversed the variant order to
check order and thermal bias. The table reports the arithmetic mean of the
three Criterion point estimates.

| Input length | Prepared batch-affine | Direct projective | Prepared advantage |
| ---: | ---: | ---: | ---: |
| 16 | 6.213 ms | 6.809 ms | 8.76% |
| 32 | 14.172 ms | 15.373 ms | 7.81% |
| 64 | 31.954 ms | 34.357 ms | 6.99% |
| 128 | 71.269 ms | 76.042 ms | 6.28% |
| 256 | 157.423 ms | 166.890 ms | 5.67% |

The prepared route won every individual A/B run at every tested length. In the
reversed-order run it retained advantages of 8.94%, 7.96%, 7.14%, 6.40%, and
5.77% from length 16 through 256. This supports general native-AVX adoption,
not merely the original BTX batch sizes. The diminishing percentage is
consistent with normalization and line preparation removing a bounded share
of total FFT and pairing work as the transform grows.

The decision is therefore to retain prepared batch-affine pointwise pairing as
the sole production native AVX FFT path. The direct projective implementation
is not deleted because it remains the scalar specialization. Results are
deployment-CPU evidence, not a universal microarchitecture guarantee; a future
AVX implementation change must repeat the same serialized forced-path A/B
before changing this policy.

The exact policy header rebuilds with Peach Clang 21.1.8, and all 11 focused
Release Rust middle-product tests pass with native AVX enabled. The policy
commit changes no arithmetic, allocation, transform shape, or runtime branch.
