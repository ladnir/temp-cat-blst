# AVX Backend

This directory contains the promoted AVX backend, following the structure
described in `docs/backend-refactor-plan.md`.

## Current State

- the promoted AVX core lives under `cpu/`, `rns/`, `fields/`, `curve/`, `pairing/`, `gt/`, and `fft/`
- build wiring exists via `blst/avx/CMakeLists.txt`
- copied AVX base-field headers are now consumable through
  `blst_avx_headers`
- the AVX-owned tests and benchmarks now live under
  `test/avx/` and `benchmark/avx/`

## Directory Boundary Rule

Copied upstream AVX core must stay in these paths:

- `blst/avx/cpu/vector/`
- `blst/avx/cpu/reduction/`
- `blst/avx/cpu/precompute/`
- `blst/avx/rns/`

`blst`-owned AVX integration code must stay in these paths:

- `blst/avx/adapters/`
- `test/avx/test_data/`

## Import Rule

During the literal import phase, changes to copied upstream files are limited to:

- include-path adjustments
- build glue
- symbol collision fixes
- mechanical instruction-binding changes needed to route vector operations
  through the `blst/backend/avx512/` wrapper layer

No algorithm, scheduling, reduction-order, or representation changes belong in
this directory until parity is established.
