# AVX batch subgroup checks: current design and earlier measurements

## Current design

The friendly G1 check keeps 81 ternary rows above its 81-point direct-check
cutoff. The G1 cofactor has a factor 3, so the three coefficients are already
the largest consecutive symmetric alphabet with a one-cancellation-per-row
argument on every nonzero cofactor component.

The friendly G2 check uses direct relations through 51 points and 22 rows
with coefficients `0, ±1, ..., ±6` above that. Its randomized path uses a
caller-owned 12-point bucket workspace. Both groups also use the direct path
when an input is semantically infinity but its physical RNS `Z` is not exactly
zero. This avoids a false rejection by the AVX complete-add facade. A fresh
32-byte seed is required even for a nonempty direct check. The G2 challenge
stream has its own SHA-256 domain and uses rejection sampling for uniform
coefficients.

The G2 cofactor is coprime to every integer from 1 through 12. Condition on
all challenges except one input whose quotient component is nonzero. For
either alphabet, at most one coefficient can cancel that component in a row.
The current randomized G2 path has an ideal false-accept probability at most
`13^-22 < 2^-80`; the SHA-256 PRF distinguishing advantage is additional.
For ordinary points, the direct path uses the already-proved equivalence
between prevalidation plus the native AVX subgroup relation and full
individual validation. It handles semantic infinity explicitly.

`CurveBatchValidation.lean` derives the selected G2 quotient's cardinality
from the ambient curve-order identity and the previously proved exact order
of the selected generator. It proves coefficient injectivity on that quotient,
the conditioned 22-row uniqueness statement, and `13^22 > 2^80`. It also
retains the earlier 35- and 41-row statements. The ambient curve-order
identity, PRF model, and physical AVX addition and doubling correspondence
remain explicit proof boundaries.

## Peach performance measurements

On 2026-11-02, the restored `master` at `05d5a4f6` and the initial direct-check
candidate (with a 22-point cutoff) were built with Clang 21.1.8, native
AVX-512/IFMA, and Release settings. Google Benchmark ran on CPU 15 with seven
repetitions at each count. Baseline and candidate ran one after the other under
`/tmp/hypercat-global-benchmark.lock`; each entry is the median elapsed time
for an accepted batch of consecutive generator multiples.

| G2 points | Baseline batch | Direct-check candidate | Change |
| ---: | ---: | ---: | ---: |
| 1 | 576.569 µs | 18.475 µs | 96.8% faster |
| 8 | 632.836 µs | 147.742 µs | 76.7% faster |
| 16 | 695.412 µs | 295.360 µs | 57.5% faster |
| 22 | 742.328 µs | 406.082 µs | 45.3% faster |
| 23 | 750.313 µs | 745.064 µs | similar |
| 256 | 2571.929 µs | 2543.920 µs | similar |

The direct path has no statistical false-accept term once on-curve validation
and the native subgroup relation agree with their field interpretation.

### Direct-path crossover refinement

The Peach scratch build used the same native AVX candidate with its direct
cutoff temporarily raised to 80 for the coarse sweep and 52 for the exact
49–52 sweep. The matching restored-master build retained the randomized path.
Both used the Release settings and CPU 15 above. The two benchmark processes
ran serially under the shared lock, with seven repetitions
per size. Each number is the median elapsed time for an accepted generator
batch; neither time includes point construction.

| G2 points | Direct relations | 22-row buckets | Direct advantage |
| ---: | ---: | ---: | ---: |
| 40 | 742 µs | 880 µs | 15.7% |
| 44 | 817 µs | 911 µs | 10.3% |
| 48 | 891 µs | 941 µs | 5.3% |
| 49 | 910 µs | 947 µs | 3.9% |
| 50 | 929 µs | 955 µs | 2.7% |
| 51 | 948 µs | 962 µs | 1.5% |
| 52 | 966 µs | 970 µs | 0.4% |
| 56 | 1039 µs | 1004 µs | -3.5% |
| 60 | 1113 µs | 1034 µs | -7.6% |
| 64 | 1188 µs | 1065 µs | -11.5% |

The production cutoff is 51, leaving the near tie at 52 on the bucket path.
This is a Peach-specific tuning choice; the security argument for either
branch does not depend on the cutoff.

## Earlier Peach measurements

`bench.cpp` in the former Peach scratch checkout measured the median of seven
samples, each with five calls, on accepted batches of consecutive generator
multiples. The old ternary and experimental 35/41-row binaries ran serially
under `/tmp/hypercat-global-benchmark.lock` on 2026-10-30.

| G2 points | Old ternary batch | Experimental 35/41-row batch | Reduction |
| ---: | ---: | ---: | ---: |
| 256 | 8.476 ms | 4.963 ms | 41.4% |
| 8192 | 164.797 ms | 125.900 ms | 23.6% |

These numbers do not compare against the current 22-row bucket implementation
and are not evidence that the earlier variant should replace it. The G1
256-point control was 5.276 ms before the challenge refactor and 5.241 ms
afterward, within the noise expected for unchanged ternary work.

Protected comparisons measured 58.591 ms versus 60.950 ms for the
thirteen- and nine-value variants at 4096 points, and 92.638 ms versus
92.219 ms at 6144 points. At 8192 points they measured 132.958 ms versus
128.330 ms. The nine-value variant uses a three-multiple table instead of a
five-multiple table. Around the 50-point cutoff and the 4096-point width
switch, timings are hardware dependent; these constants are Peach-tuned.

## Paper-derived G1 idea

[Divide-and-Pair](https://eprint.iacr.org/2026/749) studies division and
pairing trade-offs for cofactor membership tests. BLS12-381 G1 has a rational
point of order 3 at `(0,2)`, suggesting a cubic-character filter based on
`y-2`. A scalar-field exponentiation by `(p-1)/3` measured 16.943 µs under
the shared benchmark lock. Applying one filter per point would add about
4.3 ms at 256 points before any wider batch rows. This route did not justify
an implementation. The paper does not benchmark BLS12-381 G1 or G2.
