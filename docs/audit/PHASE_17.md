# Phase 17 - FFT and Middle Product

## State

Phase 17 is in progress. Packet 17A reduces and separates the FFT and
interpolation boundaries before the deeper file-by-file algorithm audit. It
does not mark any FFT or middle-product production file complete.

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
