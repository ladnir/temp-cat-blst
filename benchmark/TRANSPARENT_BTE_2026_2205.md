# Transparent BTE performance experiment

This report compares the online algebra in [*DKG Is All You Need* (ePrint 2026/2205)](https://eprint.iacr.org/2026/2205.pdf) with the existing BTX implementation in this repository. DKG is outside the comparison. The paper's uninstantiated simulation-extractable NIZK is excluded from both timing and security claims.

## Method

Measurements ran sequentially on Peach, pinned to CPU 2 of an AMD Ryzen 9 7950X. The build used Clang 21.1.8, Release, `-O3 -march=native`, CatBLST assembly and ADX, native AVX-512 types, and the AVX backend as the default. No two benchmark processes ran concurrently.

Unless stated otherwise, each reported value is the median of five repetitions; the `B=512` measurements use three repetitions. Google Benchmark measured wall time. Ciphertext randomness, sharing coefficients, and Cauchy weights use full-width `Fr` values rather than favorable small integers. `Paper best` uses the constant one-slot case at `B=1`, direct multiplication at `B=2`, and one native projective G1 MSM plus one native projective G2 MSM per Cauchy row from `B=3` through the largest measured size, `B=512`. This selection is benchmark-local because its crossover depends on the machine and backend.

`Total decrypt` includes ciphertext-only precomputation, post-share finalization, GT serialization, and a 32-byte HKDF output per slot. It excludes ciphertext proof verification, share collection, and share verification for both schemes.

## End-to-end decryption

| Batch | Paper best | Existing BTX | Paper / BTX |
|---:|---:|---:|---:|
| 1 | 0.143 ms | 0.140 ms | 1.03x |
| 2 | 0.769 ms | 0.726 ms | 1.06x |
| 4 | 1.749 ms | 1.576 ms | 1.11x |
| 8 | 4.315 ms | 3.505 ms | 1.23x |
| 16 | 12.033 ms | 7.841 ms | 1.53x |
| 32 | 38.099 ms | 17.450 ms | 2.18x |
| 50 | 105.642 ms | 36.261 ms | 2.91x |
| 64 | 159.689 ms | 38.560 ms | 4.14x |
| 128 | 526.845 ms | 84.985 ms | 6.20x |
| 256 | 1,820.196 ms | 186.200 ms | 9.78x |
| 512 | 5,969.554 ms | 403.300 ms | 14.80x |

At `B=256`, native row MSM cuts the paper implementation's total time from about 2.80 seconds to 1.82 seconds. BTX remains 9.78 times faster overall. The ciphertext-only Cauchy transform remains the cause: 1,824.850 ms for native row MSM versus 157.631 ms for BTX's prepared middle product, an 11.6-fold difference. Post-share finalization is much closer.

## Cauchy-transform implementations

| Batch | Repeated multiply | Native row MSM | Lemma 1 | BTX middle product |
|---:|---:|---:|---:|---:|
| 8 | 4.694 ms | 3.200 ms | 8.002 ms | 2.631 ms |
| 32 | 63.229 ms | 33.705 ms | 140.167 ms | 13.992 ms |
| 50 | 150.708 ms | 100.948 ms | 260.000 ms | 30.894 ms |
| 64 | 244.599 ms | 153.197 ms | 403.529 ms | 31.669 ms |
| 128 | 961.735 ms | 520.086 ms | 1,075.371 ms | 71.112 ms |
| 256 | 3,815.926 ms | 1,824.850 ms | 2,766.648 ms | 157.631 ms |
| 512 | — | 5,888.205 ms | 6,970.342 ms | — |

The native row-MSM path is the best measured implementation of the paper's transform. It reduces the direct path by 46.7% at `B=32`, 45.9% at `B=128`, and 52.2% at `B=256`. The selected CatBLST kernels use shared-doubling schedules for small public MSMs and Pippenger for larger rows. Points and scalars remain in native AVX and Montgomery representations; one aligned scratch buffer is reused across all rows.

The Lemma 1 algorithm is asymptotically better, but its group-polynomial transforms and arbitrary-point remainder trees have large practical constants. It lost to native row MSM at every measured size through `B=512`. The actual crossover is therefore above 512 on this machine; it was not extrapolated into a claimed number.

The prepared-pairing quadratic experiment was also rejected. Although it eliminates quadratic G2 scalar multiplication and reuses prepared G2 lines, its wide Miller loops cost 66.861 ms at `B=32` and 264.390 ms at `B=64`, slower than both native row MSM and repeated multiplication.

### Cost of full-width Cauchy weights

The native row-MSM path was instrumented to distinguish construction of the
Cauchy weights from their use by the G1 and G2 MSM kernels. A control run used
the same ciphertexts, allocations, row shape, affine conversion, and pairings.
It first constructed the real Cauchy matrix, then replaced each non-diagonal
weight by one immediately before the MSMs. The control output is not a valid
protocol result; it isolates the effect of scalar magnitude and density.

| Batch | Weights | G1 MSM | G2 MSM | Pairings | Total |
|---:|:---|---:|---:|---:|---:|
| 16 | Cauchy | 2.149 ms | 4.723 ms | 2.874 ms | 9.875 ms |
| 16 | Unit control | 0.382 ms | 0.635 ms | 2.874 ms | 4.019 ms |
| 32 | Cauchy | 8.307 ms | 19.348 ms | 5.753 ms | 33.767 ms |
| 32 | Unit control | 1.528 ms | 2.547 ms | 5.750 ms | 10.182 ms |
| 50 | Cauchy | 21.640 ms | 67.850 ms | 8.987 ms | 99.143 ms |
| 50 | Unit control | 3.738 ms | 19.995 ms | 8.986 ms | 33.392 ms |

Constructing the weights itself costs only 0.017 ms, 0.067 ms, and 0.218 ms
at these three sizes and is unchanged in the control. The full-width weights
increase combined G1 and G2 MSM time by 5.86 ms at `B=16`, 23.58 ms at `B=32`,
and 65.76 ms at `B=50`. They increase the complete precomputation by 2.46x,
3.32x, and 2.97x, respectively. At `B=50`, the unit-weight control is within
8.1% of BTX's 30.894 ms middle product, whereas the actual paper transform is
3.21x slower. The G2 kernel changes from the shared-doubling implementation to
Pippenger above 32 points, which explains the discontinuity in the unit-weight
G2 control between `B=32` and `B=50`.

An uninstrumented repeat measured 9.954 ms, 33.965 ms, and 99.109 ms for the
real MSM path, confirming that the phase timers do not materially change the
result.

## Other online phases at B=256

| Operation | Paper | Existing BTX | Observation |
|:---|---:|---:|:---|
| Encrypt one ciphertext | 0.163 ms | 0.0738 ms | Paper is 2.21x slower |
| Produce one party's batch share | 0.0630 ms | 6.630 ms | Paper is 105x faster |
| Verify one paper share | 0.276 ms | — | Two-pair equation |
| Verify threshold shares / existing aggregate | about 2.21 ms for 8 paper shares | 13.612 ms | Different verification workflows |
| Combine 8 shares | 0.108 ms | 0.275 ms | Paper is 2.55x faster |
| Finalize and derive 256 keys | 34.350 ms | 28.133 ms | Paper is 1.22x slower |

The paper's strongest performance result is its share path. A party publishes one G1 point after one scalar multiplication of the sum of all first ciphertext components, whereas the current BTX share path scales its cryptographic work with the batch. This advantage is useful but does not compensate for the paper's arbitrary-point cross-term transform at the measured batch sizes.

## Interpretation

The paper construction primarily buys transparent, batch-size-independent public parameters and arbitrary batch sizes. Existing BTX buys substantially lower online decryption latency by using structured public parameters and a structured middle product.

At `B=256`, the optimized paper implementation is still not a performance replacement for BTX. It is, however, substantially closer than the first implementation, and at very small batches (`B=2` and `B=4`) the measured decryption totals are within about 6–11%.

This conclusion is limited to the implemented algebraic core. A complete comparison still needs a concrete, reviewed NIZK for the paper construction and the corresponding ciphertext-admission cost for both systems.

## Validation

- The direct, native-MSM, prepared-pairing, and Lemma 1 Cauchy transforms agree at batch sizes 1, 2, 3, 8, 17, and 33.
- Threshold round trips and malformed-share rejection pass on native AVX and scalar builds.
- The complete native test binary ran 363 tests: 358 passed, 4 vector tests were skipped because their external vector file was not configured, and one prepared-G2 C-API test failed. The same test fails at the untouched base commit `05d5a4f6` under the identical native build configuration.
