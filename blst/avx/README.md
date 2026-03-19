# VROOM Landing Zone

This directory is the Phase 2 landing zone for the literal in-tree VROOM import
described in `docs/vroom-backend-reset-plan.md`.

## Current State

- the literal upstream `cpu/`, `src/`, and `test_data/` trees are present
- build wiring exists via `blst/avx/CMakeLists.txt`
- copied VROOM base-field headers are now consumable through
  `blst_vroom_headers`
- the first `blst`-owned base-field tests and benchmarks live under
  `blst/avx/tests/` and `blst/avx/bench/`

## Directory Boundary Rule

Copied upstream VROOM core must stay in these paths:

- `blst/avx/cpu/vector/`
- `blst/avx/cpu/reduction/`
- `blst/avx/cpu/precompute/`
- `blst/avx/src/`
- `blst/avx/test_data/`

`blst`-owned code must stay in these paths:

- `blst/avx/adapters/`
- `blst/avx/tests/`
- `blst/avx/bench/`

## Import Rule

During the literal import phase, changes to copied upstream files are limited to:

- include-path adjustments
- build glue
- symbol collision fixes
- mechanical instruction-binding changes needed to route vector operations
  through the `blst/backend/avx512/` wrapper layer

No algorithm, scheduling, reduction-order, or representation changes belong in
this directory until parity is established.
