# Phase 17 - FFT and Middle Product

## State

Phase 17 is in progress. Packet 17A reduces and separates the FFT and
interpolation boundaries before the deeper file-by-file algorithm audit. It
does not mark any FFT or middle-product production file complete.

## Packet 17B - Middle-product safety and public contracts

The first middle-product file pass found that public length arithmetic was not
bounded before FFT-size growth, span construction, or allocation. Because the
shared FFT uses signed `int` indices, its largest supported power-of-two size is
`2^30`; the optimized middle product therefore accepts `1 <= n <= 2^29`.
`AUD-283` replaces unchecked growth with explicit `invalid_argument` and
`length_error` failures, overflow-free admissible-size rounding, and C-boundary
rejection before constructing the `2*n` input span. Hostile-size tests exercise
the exact supported boundary, one past it, and `SIZE_MAX` without allocating a
large buffer.

`AUD-286` records the mathematical and cryptographic contract across native
C++, C, and safe Rust. For `a[0..n)` and `b[0..2*n)`, output `k` is
`product_i e(a[i], b[n+k-i])`. Inputs named G1/G2 must already be reduced,
on-curve prime-order subgroup points; infinity is the identity. Outputs and
`combine` operands are GT values carried in Fp12 storage. No Release membership
validation is performed. The APIs now state variable-time and allocation
behavior, output publication on failure, immutable handle ownership, workspace
concurrency, and exact-versus-partial alias rules. Result-bearing C declarations
are marked `BLST_NODISCARD`.

The same pass records `AUD-284`, `AUD-285`, and `AUD-287` through `AUD-289`:
inverse-root construction, repeated truncated-plan allocation, stale backend-
independent Auto selection, avoidable AVX result serialization, and cleanup/
provenance. Packets 17C through 17E close all recorded findings from this pass.

Focused validation exposed two additional defects. `AUD-290` found that the
scalar-to-AVX projective bridge preserved scalar's all-zero infinity rather
than returning MatrixNoK's canonical identity. Addition accepted that
representative, but scalar multiplication did not, so an infinity supplied
through the scalar ABI could corrupt the AVX G1 FFT. The bridge now returns
canonical G1/G2 AVX infinity before coordinate conversion. `AUD-291` corrects
backend-explicit algorithm tests that had silently called default-backend
wrappers in their AVX branches. The strengthened C reconstruction test compares
finite, zero-constant, and update contributions independently against scalar
results, and direct bridge tests require both canonicalized identities to
remain infinity after AVX multiplication.

Focused validation at implementation commit `5e14afc` used MSVC 19.50 in a
scalar-default Debug build that also compiles portable AVX. Seven sizing,
boundary, scalar-differential, exact-alias, and bridge tests pass; the corrected
explicit AVX FFT comparison also passes. The safe Rust middle-product module
passes 8/8. The 25-family public-documentation ledger and group-contract ledger
pass. An older portable-Debug executable reproduces `AUD-290`, while an older
scalar executable passes the same reconstruction, preserving before/after and
backend-isolation evidence. No benchmark, remote process, CUDA, sanitizer,
portable-AVX Release, full test sweep, or BTx/BTE audit was performed.

## Packet 17C - Reusable precomputation and typed GT boundary

Commit `f39063f` removes the repeated offline and per-evaluation setup recorded
by `AUD-284` and `AUD-285`. A complete ordered root table now produces inverse
roots as `[roots[0], roots[n-1], ..., roots[1]]`; no per-root field inversion is
performed. One shape-only truncated-forward plan and one fixed middle-window
inverse plan are built with the immutable B precompute. Both B transforms and
every later A transform read the same forward plan, while every inverse reads
the same inverse plan. Concurrent evaluations remain read-only with respect to
the precompute and retain separate caller-owned workspaces.

The B0 transform now runs directly in its final `b_fft` storage, removing one
full temporary G2 vector without changing the butterfly or fold algorithms.
`AUD-288` replaces the 576-byte serialize/decode C-boundary conversion with
typed scalar-GT and AVX-GT overloads; the latter uses the already-audited direct
MatrixNoK-to-scalar Fp12 bridge. Exact ABI layout assertions remain in place.

`AUD-289` removes the unused backend-selection macro, empty timing scope,
duplicated G1 bridge expression, and tabbed formatting residue. The three
middle-product owners now record Apache-2.0, Category Labs, retained
Supranational C-header copyright, and BLST primitive lineage accurately.

MSVC 19.50 scalar-default plus portable-AVX Debug builds the affected library
and test target with two workers. Eight focused inverse-root, plan-reuse,
sizing, and C-boundary tests pass. Repeated scalar evaluations with separate
workspaces leave both plan-storage addresses unchanged and agree exactly.
Structural scans find no old selector macro, empty timer metric, or C-boundary
serialization/decode path. The build recompiled frozen BTx as a dependency and
reported its pre-existing nodiscard warning; no BTx/BTE file was changed or
audited. No benchmark, remote process, CUDA, sanitizer, portable-AVX Release,
or full test sweep ran.

## Packet 17D - Reusable AVX batch final-exponentiation scratch

Commit `30c7558` closes `AUD-282`. `avx::GT::final_exp_n` now requires an
explicit caller-owned `GT::FinalExpScratch` instead of constructing six fresh
vectors on every call. The scratch owns the two six-term decomposition arrays,
the two-term inverse input, the first-component terms, prefix products, and
inverses. Callers may reserve capacity, query it, and reuse the owner; the
individual arithmetic buffers remain private. Buffer sizes retain their
high-water mark, so a later smaller or equal batch performs no resize and no
batch-final-exponentiation allocation.

The reusable AVX multipair middle-product workspace owns this scratch and
passes it through every finalization. The one-shot update path uses an explicit
local scratch because it has no reusable workspace lifetime. Allocation growth
still occurs before output publication and may throw to the existing native/C
translation boundary; arithmetic, batch-inversion order, nonoverlap contract,
and delayed publication are unchanged. This packet makes no claim that the
entire middle-product operation is allocation-free.

MSVC 19.50 scalar-default plus portable-AVX Debug builds the affected test
target with two workers. Direct batch final exponentiation matches individual
results across full, smaller, repeated, and zero-length calls while retaining
capacity. Two evaluations through one multipair workspace agree and preserve
scratch capacity. The existing AVX multipair-versus-FFT comparison passes at
sizes 3, 5, 9, and 16. The ownership-only `GT.h` rebind validates all 21 AVX
correspondence rows and 61 exact mapped hashes. Frozen BTx rebuilt as a header
dependency and emitted its existing nodiscard warning; it was neither changed
nor audited. No benchmark, remote run, CUDA, sanitizer, portable-AVX Release,
or full suite ran. Only benchmark-gated `AUD-287` remains open in the current
middle-product packet.

## Packet 17E - Backend-specific Auto dispatch

`AUD-287` could not be resolved from the existing benchmark because that case
precomputed B outside the timed loop and recreated multipair workspace while
reusing FFT workspace. The production `Auto` overload is one-shot, so commits
`d9dc202` through `03457c7` add an isolated non-shipping executable that times
forced scalar/AVX FFT and multipair calls including their actual precomputation
and workspace ownership. Isolating the executable also avoids compiling or
linking frozen BTx/BTE benchmark sources.

One benchmark process at a time ran on Peach's AMD Ryzen 9 7950X, pinned to CPU
2, with Clang 21.1.8, `Release`, native AVX-512/IFMA, five repetitions, and a
0.10-second minimum case time. The focused boundary medians in milliseconds
were:

| n | Scalar FFT | Scalar multipair | AVX FFT | AVX multipair |
| ---: | ---: | ---: | ---: | ---: |
| 31 | 87.4 | 82.2 | 24.3 | 22.3 |
| 32 | 83.1 | 88.1 | 22.9 | 23.6 |
| 49 | 200 | 193 | 55.8 | 52.0 |
| 50 | 200 | 201 | 55.7 | 54.0 |
| 51 | 200 | 208 | 55.9 | 56.1 |
| 76 | 451 | 442 | 126 | 120 |
| 78 | 452 | 466 | 126 | 126 |

The nonmonotonic result is expected: the optimized FFT size is the next power
of two at least `2*n`, so FFT cost jumps after 16, 32, and 64 while multipair
cost grows smoothly. Commit `e3e1a45` therefore replaces the stale `n <= 24`
cutoff with allocation-free constexpr policies. Scalar selects multipair for
`n <= 31`, `33..49`, and `65..76`; AVX extends the second range through 50.
Every later size selects FFT. This introduces no runtime CPU dispatch, timing
probe, allocation, or data-dependent cryptographic branch; `n` is public.

MSVC 19.50 scalar-default plus portable-AVX Debug builds the affected target
with two workers. Exact range-policy tests and actual scalar and portable-AVX
Auto dispatch tests pass. Peach native Release rebuilds the isolated benchmark
against the production selector. That build exposed one inactive scalar GT
boundary overload in an AVX-default translation unit; commit `4dcb4fa` marks
only that overload `[[maybe_unused]]`, after which Clang rebuilds without the
warning. The benchmark is deployment-CPU evidence, not a universal threshold
guarantee for every microarchitecture. No concurrent benchmark, CUDA,
sanitizer, portable-AVX Release, full suite, or BTx/BTE audit ran.

## Packet 17A - Generic FFT boundary and isolated interpolation

The former `blst/fft_c_api.cpp` was an 866-line anonymous namespace combining
root generation, Fr transforms, prepared-domain storage, a general polynomial
library, product-tree interpolation, specialized suffix interpolation, scratch
partitioning, and C ABI adaptation.

The repository already had one generic FFT implementation in `ctFft.h`.
Packet 17A therefore does not add a second native `fr_fft` implementation.
Instead:

- `blst/fft/domain.h` provides compile-time generic root and prepared-domain
  utilities used by the existing FFT templates;
- `blst/fft/middle_product.h` reuses the same forward-root constructor;
- `blst/fft/fr_interpolation.h` exposes only the retained native Fr
  interpolation operations and typed statuses;
- product-tree nodes and polynomial helpers remain private to
  `blst/fft/fr_interpolation.cpp`; and
- `blst/fft_c_api.cpp` is a 169-line typed-layout, validation, span, exception,
  and C ABI adapter over the generic FFT and isolated interpolation code.

No runtime polymorphism, callback wrapper, global state, scheduler, or hidden
thread creation was introduced.

## Retained public surface

The C FFT surface contains only:

- one-shot forward/inverse Fr FFT;
- prepared-domain storage length and initialization; and
- allocation-free transform with a prepared domain.

The retained C interpolation surface contains only:

- general arbitrary-point interpolation;
- the FFT-domain suffix scratch-length query; and
- allocation-free FFT-domain suffix interpolation with caller-owned scratch.

The redundant forward-transpose alias, full FFT-domain convenience wrapper,
and allocating suffix wrapper were removed. Rust exposes the same reduced
operations under the explicit `fft` module instead of re-exporting them from
the crate root. General interpolation and FFT-domain suffix interpolation stay
separate because they use different algorithms and have different complexity.

Prepared domains retain exact `n + 2` storage: `n` forward roots, `1/n`, and an
encoded size tag. The C/Rust boundary validates the exact storage length and
tag; the generic internal helper assumes its caller already established those
invariants.

## Performance constraints

This is a structural change. It preserves the generic FFT butterflies, inverse
scaling, product-tree thresholds, general interpolation algorithm, specialized
suffix algorithm, and flat caller-owned suffix scratch layout. The hot paths
remain statically dispatched templates with no virtual calls or type erasure.

No benchmark was run. Product-tree allocation, truncated suffix FFT work,
direct-versus-FFT Fr middle-window crossover, and middle-product algorithm
selection remain later review packets and must be benchmarked one process at a
time.

## Rust package repair

The checked-in Cargo include list still named the removed upstream `src/`
layout. Repository builds worked only because `build.rs` found the surrounding
checkout; a published crate omitted the native CMake/source tree. Its broad
`**/*.rs` pattern could also follow the temporary publish link into local build
outputs.

The manifest now uses explicit crate-local Rust paths and packages the actual
root CMake project, native sources, generated assembly, public headers, license,
and frozen BTx/BTE configure dependency. Package verification uses a fresh
`CARGO_TARGET_DIR` so an existing checkout CMake cache cannot be reused for the
relocated package source tree. This repair does not audit or modify frozen
BTx/BTE code.

## Focused validation

- Dedicated scalar native build: 7/7 focused FFT, C boundary, suffix, and
  middle-product smoke tests passed.
- Portable AVX-selected build with `BLST_DEFAULT_AVX=ON`: the same 7/7 tests
  passed.
- Linux GCC ASan: the same 7/7 tests passed with leak detection and
  halt-on-error enabled.
- Rust Release FFT module: default 9/9 and portable 9/9 passed.
- Packaged Rust crate: 358 files, 4.9 MiB; all three new FFT/interpolation files
  are present, no `blst/out/**` entries are present, and archive verification
  passed from a fresh Cargo target directory.
- Audit inventory: 612 files and all 63 upstream commits are accounted for.
- `git diff --check`: passed apart from local line-ending notices.

No benchmark, full test sweep, native AVX-512 hardware run, CUDA build, or
BTx/BTE audit was performed.

## Remaining FFT decisions

1. Audit `ctFft.h` as the single transform engine, including size types, root
   shape contracts, permutations, normalization, truncated plans, and dead
   experimental code.
2. Audit the retained product-tree implementation and decide whether its
   allocation pattern should become a contiguous reusable workspace.
3. Evaluate the specialized suffix as an adaptive direct/Fr-middle-window
   operation without coupling it to the pairing-specific middle-product API.
4. Decide whether the C ABI should retain `bool` or adopt typed operation
   statuses while pre-production breakage remains acceptable.
5. Document the C/native Fr representation convention and evaluate a future
   unified storage definition.
6. Complete differential oracles, fuzzing, sanitizer, native AVX-512, and Apple
   Silicon evidence during the remaining phase.
