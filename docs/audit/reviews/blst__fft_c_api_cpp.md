# File Review: `blst/fft_c_api.cpp`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/fft_c_api.cpp` |
| Roadmap phase / queue | Phase 17 Packet 17I |
| Reviewed implementation commit | `f4063cc34c28f0375f04908bfe5c456a128562b6` |
| Reviewer | Codex audit team |
| Date | 2026-08-11 |
| Disposition | `in_scope_production` |
| Subsystem | Fr FFT, prepared domains, and interpolation C boundary |
| Profiles | scalar and selected AVX native C++; supported C and safe Rust adapters |
| Generation/provenance | Locally composed handwritten adapter; Apache-2.0 and Category Labs |
| Upstream references | `not_applicable: upstream BLST has no corresponding FFT C API` |
| Related files | `bindings/blst_fft.h`, `bindings/rust/src/fft.rs`, `blst/fft/domain.h`, `blst/fft/fr_interpolation.*`, `blst/fft/ctFft.h`, `test/BindingHeader_tests.cpp` |

## Purpose and trust boundary

This file is the C ABI boundary for seven Fr operations: ordinary FFT, prepared
domain length and initialization, prepared-domain FFT, arbitrary-point
interpolation, suffix scratch sizing, and allocation-free FFT-domain suffix
interpolation. Pointers and sizes may be untrusted. Coefficients may be
secret-derived; transform sizes, directions, points, and indices control public
work schedules.

The adapter aliases `blst_fr` directly to the selected native `blst::Fr` and
proves exact type identity, standard layout, and trivial copyability. It does
not reinterpret unrelated storage, allocate across the ABI, retain caller
pointers, or own mutable global state.

## Contracts and invariants

- All Fr operands are reduced Montgomery-form values.
- FFT sizes are nonzero powers of two no larger than `2^30`.
- Ordinary FFT prepares directional roots and may allocate before mutation.
- Prepared domains have exactly `n + 2` elements: forward roots, inverse scale,
  and an encoded size tag. Initialization publishes only on success.
- Prepared transforms require exact shape, a matching tag, and nonoverlapping
  mutable values and immutable domain storage.
- General interpolation requires equal nonempty shapes and distinct points.
- Suffix interpolation requires valid unique domain indices, `ell <= n <=
  fft_len`, exact scratch capacity, and the documented overlap rules. The C ABI
  retains non-null pointer requirements even for a logical zero-length suffix.
- No C entry point permits a C++ exception to cross the ABI. Allocation-bearing
  functions catch locally; interpolation uses native `noexcept` status paths.
- Validation failures precede output mutation. Ordinary and prepared in-place
  FFTs can be partially transformed only after an unexpected internal failure.

## Dependency and call-path map

`domain.h` owns supported sizing, root construction, and caller-owned domain
layout. `ctFft.h` owns the transform. `fr_interpolation.cpp` owns validated
arbitrary-point and suffix reconstruction and translates internal failures to
typed status. The safe Rust module calls these C functions with slice-derived
shapes. Frozen BTx/BTE does not use this public adapter and was not audited or
modified.

## Checklist evidence

### 1. Identity, API, and ABI

Every exported function and declaration was matched to C and safe Rust callers.
Signatures and storage layout remain compatible. Four fallible declarations now
carry `BLST_NODISCARD`; no operation, parameter, or result representation was
changed.

### 2. Correctness and failure publication

The adapter validates pointers and public shapes before span construction.
Prepared transforms authenticate the encoded size tag before mutation. General
and suffix interpolation delegate to previously audited status-returning
owners. Forward/inverse round trips, prepared-domain reuse, invalid-size
rejection, oversubscribed interpolation, aliased general inputs, and scratch
sizing pass focused tests.

### 3. Memory, lifetime, and undefined behavior

Exact C/native type identity makes span construction type-correct without a
copy or aliasing exception. Checked native sizing precedes indexed domain access
at the C boundary. No pointer is retained, and ownership stays with the caller.
Native narrow-contract overlap restrictions are explicit in the public header.

### 4. Exceptions and allocation

Ordinary FFT allocates its root vector under a catch-all boundary. Domain
initialization is allocation-free but retains the same defensive boundary.
Prepared FFT is allocation-free. Interpolation status APIs are `noexcept` and
contain their own allocation failures. No exception can escape to C or Rust.

### 5. Side channels and secret lifecycle

FFT schedules and addresses depend on public sizes and direction. General
interpolation is variable-time in public points; for fixed points its
reconstruction shape does not depend on values. Suffix work depends on public
sizes and indices. Ordinary heap and caller scratch can hold value-derived
intermediates and are intentionally not securely erased; that lifecycle is
documented rather than imposing generic hot-path clearing.

### 6. Concurrency and resource exhaustion

There is no global mutable state. Completed domains may be shared immutably;
mutable value, output, and scratch ranges must be disjoint between concurrent
calls. The ordinary FFT and general interpolation may allocate proportionally
to public input size; invalid sizes fail before allocation. Prepared and suffix
paths use caller-owned storage.

### 7. Portability and performance

The public header parses as strict MSVC C11 and C++20 at `/W4 /WX`. The boundary
adds no runtime dispatch, copying, or per-element validation. `BLST_NODISCARD`
and documentation have no generated-code cost. No benchmark was warranted
because arithmetic and hot control flow are unchanged.

## Findings

`AUD-025`, `AUD-031`, `AUD-032`, `AUD-034`, `AUD-043`, `AUD-046`, `AUD-297`,
`AUD-300`, `AUD-309`, and `AUD-310`.

## Validation

- MSVC 19.50 scalar Debug and Release focused boundary tests: 5/5 in each.
- Strict standalone `bindings/blst_fft.h` C11 and C++20 `/W4 /WX` parse: passed.
- Public documentation ledger: 26 audited/focused families, passed.
- Audit manifest integrity and `git diff --check`: passed at closure.
- No benchmark, CUDA, portable-AVX Release, full suite, sanitizer, native
  AVX-512, remote, Apple, or BTx/BTE audit/change was performed.

## Residual risk and decision

Milestone sanitizer, fuzz, native AVX-512, and Apple evidence remain broader
release gates. Arbitrary prepared-domain contents are not mathematically
authenticated beyond shape/tag; callers must initialize and preserve immutable
domain storage as documented. This file is complete for line-by-line hand
review.
