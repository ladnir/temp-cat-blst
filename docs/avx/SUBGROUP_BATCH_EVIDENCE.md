# AVX batch subgroup checks: current design and earlier measurements

## Current design

The friendly G1 check keeps 81 ternary rows above its 81-point direct-check
cutoff. The G1 cofactor has a factor 3, so the three coefficients are already
the largest consecutive symmetric alphabet with a one-cancellation-per-row
argument on every nonzero cofactor component.

The friendly G2 check uses direct relations through 22 points and 22 rows
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
The direct path uses the already-proved equivalence between prevalidation plus
the scalar relation and full individual validation.

`CurveBatchValidation.lean` proves coefficient injectivity on a finite
additive quotient group of cardinality equal to the concrete G2 cofactor,
the conditioned 22-row uniqueness statement, and `13^22 > 2^80`. It also
retains the earlier 35- and 41-row statements. The curve-order identity, PRF
model, and physical AVX addition and doubling correspondence remain explicit
proof boundaries.

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
