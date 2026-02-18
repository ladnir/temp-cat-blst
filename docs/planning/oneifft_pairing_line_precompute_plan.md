# One-iFFT Pairing Line Precompute Plan (Fr-NTT Correctness)

## 1) Setup and Objective

We target the online middle-product pipeline where:

- `A` has length `n` over `G1`
- `B` has length `2n` over `G2` and is fixed across many queries
- FFT length is `m = pow2_ceil(2n)`
- output is the middle window of length `n`

Current one-iFFT online path (precomputed `B`) is:

1. Truncated forward FFT of `A` in `G1`
2. `m` pointwise pairings
3. Inverse FFT in exponent domain
4. Middle-window extraction

Goal: reduce online pairing cost by exploiting fixed `D[j]` (the folded `B` spectrum), while preserving algebraic correctness for Fr twiddles.


## 2) Algebra Model (explicit)

This project uses FFT in the exponent with twiddles in `Fr` (order `r`).

1. **Target domain:** `GT` is the prime-order subgroup of `Fp12*` with order `r`.
2. **Scalar ring:** twiddles are in `Fr`; inverse FFT is an `Fr`-linear transform on `GT^m`.
3. **IFFT operations:** group op in `GT`, inverse in `GT`, and scalar multiplication by `Fr` (implemented as exponentiation).

Invariant:

> Every intermediate element fed to Fr-NTT/iFFT must lie in `GT`; scalar action is modulo `r`.


## 3) Optimization Definition (Phase 1 only)

For each frequency bin `j`, `D[j] in G2` is fixed offline while `A_hat[j] in G1` changes online.  
So we precompute pairing line data (prepared-`Q`) for each fixed `D[j]` once, and reuse it across queries.

Offline:

- build folded `D[j]` (already implemented)
- build `PreparedG2(D[j])` for each `j`

Online:

- compute `A_hat[j]`
- run Miller using prepared `D[j]`
- apply final exponentiation **per bin**
- run iFFT over `GT`

This keeps one-iFFT structure unchanged and only accelerates pointwise pairing.


## 4) Deferred Final Exponentiation: Not Allowed Here

Miller outputs lie in `Fp12*` (order `p^12 - 1`), not guaranteed in subgroup `GT` of order `r`.  
Fr-NTT requires an `Fr`-module action, which is defined for subgroup elements of order `r`.

Therefore:

- final exponentiation must occur per bin before iFFT input is formed
- defer-FE is invalid for current Fr-NTT exponent FFT pipeline

If a deferred-FE flag exists, enforce:

- if `fft_kind == ExponentFftFr` and `defer_final_exp == true`: return error/assert

Current status in this repo at plan time:

- no `defer_final_exp` mode exists yet; this guard is a required rule if such a mode is added


## 5) API and Data Model Changes

Add one-iFFT pairing mode:

- `OneIfftPairingMode::Direct`
- `OneIfftPairingMode::PrecomputedLines`

Extend one-iFFT precompute payload:

- `oneifft_d_points[j]` (optional; debug/fallback)
- `oneifft_d_prepared[j]` (prepared G2 line data)
- optional `d_is_zero[j]` mask for skip-fast path

Hot-path constraints:

- no dynamic allocations in online loop
- prepared data immutable after precompute
- workspace owns online temporaries


## 6) Correctness Argument

### Lemma 1: prepared-lines pairing equivalence

For fixed `Q in G2`, prepared Miller and direct Miller compute the same function:

- `MillerPrepared(P, Q) = MillerDirect(P, Q)` for all `P in G1`.

After final exponentiation, resulting `GT` bin values are identical in direct and prepared modes.

### Lemma 2: Fr-iFFT well-definedness

Because every bin entering iFFT is in `GT` (order `r`), twiddle scalar actions by `Fr` are valid.  
Thus iFFT computes the same transform in both modes, so extracted middle outputs match exactly.

No claim is made about moving FE past iFFT.


## 7) Implementation Checklist

1. Add backend wrapper type for prepared G2 (opaque per pairing backend).
2. Extend one-iFFT precompute object with prepared vectors.
3. Implement `PrecomputedLines` online path:
   - Miller with prepared `D[j]`
   - immediate final exponentiation
   - unchanged iFFT + extract
4. Keep `Direct` path untouched for fallback and A/B benchmarking.
5. If any future deferred-FE option appears, hard-reject for Fr exponent FFT.


## 8) Test Plan

Add/extend tests in `test/Fft_tests.cpp`:

1. `Fft.OneIfftLines_EqualsDirect_SmallDeterministic`
2. `Fft.OneIfftLines_EqualsDirect_Randomized`
3. `Fft.OneIfftLines_PrecomputeReuse_ManyQueries`
4. `Fft.OneIfftLines_BoundaryN` with `n in {1,2,3,63,64,65,127,128,129,511,512}`
5. `Fft.OneIfftLines_WorkspaceNoAllocHotPath`
6. `Fft.OneIfftRejectsDeferredFinalExp` (only if flag/mode exists)
7. subgroup sanity (debug): after FE, each bin satisfies subgroup check if API exists

Also keep oracle cross-check (`O(n^2)` middle-product helper) for small sizes.


## 9) Benchmark and Metrics Plan

Primary benchmark target:

- `bench_middle_product_one_ifft_detailed`

Counters:

- `pairing_direct_ms` / `pairing_precomp_ms`
- `precompute_lines_once_ms`
- `precompute_lines_bytes` (if measurable)
- existing `fft_g1_ms`, `ifft_gt_ms`, `fill_ms`, `extract_ms`, `online_total_ms`

Evaluate:

1. online latency only
2. amortized latency over query count `Q in {1,2,4,8,16,32,64,128}`
3. break-even `Q*` for line precompute


## 10) Risks and Mitigations

Risk: prepared line tables increase memory footprint.  
Mitigation: mode is opt-in and reports memory counters.

Risk: backend pairing API lacks stable prepared-`Q` representation.  
Mitigation: adapter layer + direct fallback.

Risk: regressions on small `n`.  
Mitigation: no auto-select yet; benchmark-guided selection.


## 11) Definition of Done

Done when all are true:

- prepared-lines mode matches direct mode on deterministic + randomized tests
- online hot path remains allocation-free with workspace
- benchmark report includes pairing/precompute breakdown and amortized crossover
- no invalid deferred-FE path is exposed for Fr exponent FFT
