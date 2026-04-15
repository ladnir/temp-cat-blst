# PFE vs BTX Cost Deep Dive

This note compares the current `BTE/pfe` and `BTE/btx` implementations at the
component level, using the prepared-path cost model:

- CRS is already fixed.
- Batch-size-specific prepared caches are already available.
- "Offline" means batch precomputation / FFT-heavy work.
- "Online" means finalize once offline precompute is done.

## Executive Summary

The corrected measurements now line up much better with the operation counts.

At `n = 512` in `x64-Release-Scalar`:

| Phase | PFE | BTX | Ratio |
|---|---:|---:|---:|
| Offline prepared | `2.428 s` | `1.489 s` | `1.63x` slower |
| Online prepared | `1.237 s` | `0.315 s` | `3.93x` slower |
| End-to-end prepared | `3.889 s` | `1.927 s` | `2.02x` slower |

Those results are consistent with the implementation-level operation counts:

- `BTX` online does about `n` pairings.
- `PFE` online does about `4n` pairings.
- `BTX` offline does fewer FFTs than `PFE`.

The main mismatch in earlier numbers was benchmark contamination: the first
timed `BTX` prepared finalize call was lazily building the alpha cache, so the
"online" benchmark was accidentally charging one-time precomputation to the
online phase.

## Component Cost Model

Notation:

- `n` = batch size
- `m = next_pow2(2n)`
- for power-of-two `n`, `m = 2n`

### BTX Prepared Path

#### Offline / beta precompute

Current code path:

- `BTE/btx/core.cpp::compute_cross_terms_prepared_fft`
- `blst/fft/middle_product.h::middle_product(...)`

Costs per batch:

- FFTs:
  - `1` forward `G1` FFT of size `m`
  - `1` truncated inverse `GT` FFT of size `m`
- Pairings:
  - `m` pointwise pairings
  - for power-of-two `n`, this is `2n`
- MSMs:
  - `0`

#### Online / finalize

Current code path:

- `BTE/btx/core.cpp::FinalizeBatchDecryption(prepared, ...)`
- `blst::update(...)`

Costs per batch:

- FFTs:
  - `0`
- Pairings:
  - `n`
- MSMs:
  - `0`

#### One-time prepared cache cost per batch size

Current code path:

- `ensure_precomputed_beta`
- `ensure_precomputed_alpha`
- `precompute_middle_product_b(...)`

Costs per cached batch size:

- beta cache:
  - `2` forward `G2` FFTs of size `m`
- alpha cache:
  - `2` forward `G2` FFTs of size `m`
- pairings:
  - `0`

So once caches are warmed, the steady-state prepared batch cost is:

- offline: `2` FFTs + `m` pairings
- online: `n` pairings

### PFE Prepared Path

#### Offline / batch precompute

Current code path:

- `BTE/pfe/core.cpp::PrecomputeBatchDecryption(prepared, ...)`

Costs per batch:

- FFTs:
  - `1` forward `G1` FFT of size `m`
  - `1` truncated inverse `G1` FFT of size `m`
  - `1` forward `GT` FFT of size `m`
  - `1` truncated inverse `GT` FFT of size `m`
- Pairings:
  - `n`
  - these are the `tm_self[j] = e(c1[j], dk[j])` terms
- MSMs:
  - `0`

#### Online / finalize

Current code path:

- `BTE/pfe/core.cpp::FinalizeBatchDecryption(...)`

Costs per slot:

- `e(l1_input, dk[j])`
- `e(c1[j], weighted_dk_excluding_slot[j])`
- `e(c2[j], g2_gen)`
- `e(c1[j], sum_dk)`

So per batch:

- FFTs:
  - `0`
- Pairings:
  - `4n`
- MSMs:
  - `0`

#### One-time setup preprocessing

Current code path:

- `BTE/pfe/core.cpp::Setup`

Costs once for fixed batch size:

- `1` scalar-kernel FFT of size `m`
- `1` forward `G2` FFT of size `m`
- `1` truncated inverse `G2` FFT of size `m`
- pairings:
  - `0`

So once setup is fixed, the steady-state prepared batch cost is:

- offline: `4` FFTs + `n` pairings
- online: `4n` pairings

### Threshold Side

For both schemes:

- create partial share:
  - `1` MSM of size `n`
- combine threshold shares:
  - `1` MSM of size `t`

No FFTs or pairings occur in those two steps.

## Expected Ratios From a Simple Heuristic

Using the working heuristic:

- `2n` pairings `≈` `1` FFT of size `2n`

and assuming `n` is a power of two so `m = 2n`:

### Offline

- `BTX`: `2` FFTs + `2n` pairings `≈ 3.0` FFT-equivalents
- `PFE`: `4` FFTs + `n` pairings `≈ 4.5` FFT-equivalents

Expected ratio:

- `PFE / BTX ≈ 4.5 / 3.0 = 1.5x`

### Online

- `BTX`: `n` pairings
- `PFE`: `4n` pairings

Expected ratio:

- `PFE / BTX ≈ 4x`

### End-to-End

- `BTX`: `3.5` FFT-equivalents
- `PFE`: `6.5` FFT-equivalents

Expected ratio:

- `PFE / BTX ≈ 6.5 / 3.5 ≈ 1.86x`

## Where the Earlier Measurements Went Wrong

The earlier head-to-head "online prepared" benchmark did not prewarm the `BTX`
alpha cache.

In `BTE/btx/core.cpp`, the prepared finalize path calls:

- `ensure_precomputed_alpha(...)`
- which lazily builds the alpha middle-product cache on first use

That cache build performs one-time FFT work. If the benchmark fixture is
constructed and the first timed iteration immediately calls prepared finalize,
that one-time FFT cost gets charged to the "online" phase.

This especially distorted large-batch runs where each repetition performed only
one timed iteration.

The benchmark harness has been updated to prewarm that cache before timing:

- `benchmark/PfeBtx_head_to_head_bench.cpp`
- `benchmark/Btx_bench.cpp`

## Corrected Measurements

### Corrected Online Prepared

`x64-Release-Scalar`, 3 repetitions, aggregates only:

| Batch | PFE | BTX | Ratio |
|---|---:|---:|---:|
| `32` | `79.47 ms` | `18.45 ms` | `4.31x` |
| `64` | `154.42 ms` | `35.97 ms` | `4.29x` |
| `128` | `310.76 ms` | `105.97 ms` | `2.93x` |
| `512` | `1237.18 ms` | `314.83 ms` | `3.93x` |

Takeaway:

- After cache prewarm, the online phase is now very close to the expected `4x`
  gap.
- The `128` run is noisier than the others, but it still strongly favors `BTX`.

### Corrected End-to-End Prepared

`x64-Release-Scalar`, 3 repetitions, aggregates only:

| Batch | PFE | BTX | Ratio |
|---|---:|---:|---:|
| `32` | `183.05 ms` | `124.58 ms` | `1.47x` |
| `64` | `393.22 ms` | `304.87 ms` | `1.29x` |
| `128` | `837.27 ms` | `503.75 ms` | `1.66x` |
| `512` | `3889.57 ms` | `1926.90 ms` | `2.02x` |

Takeaway:

- Once the online benchmark is cleaned up, the end-to-end ratios move much
  closer to the expected `~1.86x`.
- The `512` result is slightly worse for `PFE` than the simple heuristic
  predicts, but it is in the same neighborhood.

### Offline Prepared

`x64-Release-Scalar`, 3 repetitions, aggregates only:

| Batch | PFE | BTX | Ratio |
|---|---:|---:|---:|
| `512` | `2428.25 ms` | `1489.18 ms` | `1.63x` |

Takeaway:

- The offline gap is slightly worse than the `1.5x` FFT-equivalent model, but
  still reasonably close.

## Interpretation

The corrected data supports the original intuition:

- `BTX` should be materially faster online because it uses about `1/4` as many
  pairings as `PFE` in prepared finalize.
- `BTX` should also be faster overall because it uses fewer FFTs offline.

The implementation now matches that story:

- `BTX` is clearly faster online.
- `BTX` is faster offline.
- `BTX` is faster end-to-end.

`PFE` is still a useful comparison baseline, but with the current implementation
and cost structure it is not the faster scheme.

## Commands Used

Scalar release benchmark build:

```powershell
cmake --preset x64-Release-Scalar
cmake --build out/build/x64-Release-Scalar --target bte_pfe_benchmarks --parallel 1
```

Corrected online prepared:

```powershell
out/build/x64-Release-Scalar/benchmark/bte_pfe_benchmarks.exe --benchmark_filter="BM_PfeVsBtx_Online_.*/(32|64|128|512)/real_time" --benchmark_min_time=0.12s --benchmark_repetitions=3 --benchmark_report_aggregates_only=true
```

Corrected end-to-end prepared:

```powershell
out/build/x64-Release-Scalar/benchmark/bte_pfe_benchmarks.exe --benchmark_filter="BM_PfeVsBtx_EndToEnd_.*/(32|64|128|512)/real_time" --benchmark_min_time=0.12s --benchmark_repetitions=3 --benchmark_report_aggregates_only=true
```

Offline prepared at `512`:

```powershell
out/build/x64-Release-Scalar/benchmark/bte_pfe_benchmarks.exe --benchmark_filter="BM_PfeVsBtx_(Offline|Online|EndToEnd)_.*/512/real_time" --benchmark_min_time=0.15s --benchmark_repetitions=3 --benchmark_report_aggregates_only=true
```
