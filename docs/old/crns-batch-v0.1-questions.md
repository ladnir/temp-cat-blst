# CRNS Batch v0.1 Questions

## Definitions

- `Literal VROOM port`: the vendored VROOM arithmetic code compiled inside `blst` with only thin `blst` boundary adapters around it.
- `Direct benchmark`: a benchmark that calls the vendored VROOM arithmetic object directly, without scalar `blst` encode or decode around the hot path.
- `Wrapped benchmark`: a benchmark that starts from `blst` types, converts once at the boundary, and then calls the vendored VROOM arithmetic.

## Decisions Taken

1. The hand-ported `Fp12Batch.*` path was deleted from the active build.
2. `crns_batch` now means the literal vendored VROOM path only.
3. The current milestone stops at `Fp12 mul`.
4. The hot path stays source-faithful to VROOM instead of being rewritten into `blst` style.
5. Tests and benchmarks are attached to the vendored path, not to a parallel hand port.

## Current Active Surface

- Base field:
  - `blst/crns_batch/VroomLiteral.h`
  - `VroomFpPack8`
  - `mul`
- Fp12:
  - `blst/crns_batch/VroomLiteral.h`
  - `fp12_from_blst`
  - `fp12_to_blst`
  - `fp12_mul`
- Tests:
  - `test/Vroom_vendor_tests.cpp`
- Benchmarks:
  - `benchmark/Vroom_vendor_bench.cpp`

## Current Server Results

- Host: `18.188.231.170`
- CPU: Intel Xeon Platinum 8488C
- Benchmark mode: pinned single core with `taskset -c 0`

### Base Field

- In-tree vendored direct benchmark:
  - `BM_VroomVendorBatchModMul8Direct = 193 ns`
- Upstream VROOM comparison:
  - `BM_BatchModMul_8_52 = 192 ns`

### Fp12

- In-tree vendored direct benchmark:
  - `BM_VroomVendorFP12MulDirect = 1813 ns`
- In-tree vendored wrapped benchmark:
  - `BM_VroomVendorFP12MulWrapped = 1814 ns`
- Upstream VROOM minimal benchmark built on the same server with the same compiler family:
  - `BM_FP12_Mul_Original = 1917 ns`

## What This Means

- Base-field literal parity is achieved.
- `Fp12 mul` literal parity is achieved for this server toolchain.
- The active backend is now suitable as the reference implementation for any later throughput-specific transpose or pack work.

## Questions To Revisit Later

1. Should the next step stay literal and bring over `Fp12 sqr` and `Fp12 cyclotomic sqr` before any throughput-specific layout work?
2. Should the future throughput path be built as a second layer on top of this vendored arithmetic, or as a second vendored subtree with a different physical layout?
3. When throughput work starts, should the benchmark target be packed independent `Fp12` operations first, or a pairing-specific inner kernel first?
