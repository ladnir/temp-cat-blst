# AVX Constant-Time Add And GT Fix Plan

This document captures the agreed implementation plan for the AVX curve-add
correctness work and the separate `GT::doubling` performance fix.

## Status

- Patch 1, Patch 2, and Patch 3 are implemented in the current working tree.
- Focused AVX regression tests and the full `x64-Debug-Avx` test suite pass.
- A serialized GT benchmark run is still recommended before merging so we can
  record the `GT::doubling` before/after ratio on the target machine.

## Goals

- Make AVX G1/G2 `add` and `dadd` entry points constant-time with respect to
  point values where the API is expected to be constant-time.
- Preserve scalar BLST's identity conventions:
  - projective identity is `Z == 0`
  - affine identity sentinel is `(X == 0, Y == 0)`
- Keep Pippenger MSM variable-time. Its scalar-window bucket access pattern is
  part of the existing vartime MSM contract.
- Make AVX affine `cneg` mask-based instead of branch-based.
- Fix the independent `GT::doubling` performance anomaly caused by converting a
  cyclotomic square to standard form through a full multiply by one.

## Non-Goals

- Do not make Pippenger MSM constant-time.
- Do not change external C API signatures.
- Do not change `GT::data`'s internal representation as part of the first GT
  performance patch.
- Do not replace optimized fixed-structure AVX arithmetic with generic helper
  layers that obscure code generation in hot paths.

## Decisions

1. `a4 != nullptr` is a public predicate.
   Branching on it is allowed. The selected path must still avoid branches on
   secret point data.

2. Mixed add should handle affine identity inside the low-level formula helper.
   Put the `(Q.x == 0 && Q.y == 0)` select in both `PointMixedAddUnexpanded`
   and `PointMixedAdd`, so all callers inherit the scalar convention.

3. Verify `from_affine(0,0)` behavior.
   The AVX conversion path should agree with scalar BLST that all-zero affine
   means identity. Add tests if coverage is missing.

4. Use mask/select style for point-data conditionals.
   `vec_is_zero` or equivalent zero checks may produce a mask, but the mask must
   feed a constant-time select, not a data-dependent branch.

5. Make `G1AffinePoint::cneg` and `G2AffinePoint::cneg` constant-time.
   This is useful even though MSM remains vartime, because `cneg` is a general
   conditional primitive.

6. Treat `GT::doubling` as a separate performance patch.
   The GT fix needs bounds/narrowing validation and benchmarks, so it should
   land after the add/dadd correctness patch.

## Patch 1: Constant-Time Curve Add Surface

### Scope

- `blst/avx/curve/ec.h`
- `blst/avx/curve/point.h`
- focused AVX tests under `test/avx/`

### Implementation Outline

1. Add or reuse a projective constant-time select helper.
   `ProjectivePoint::ctime_select` already exists. Prefer using it directly
   instead of adding a wrapper unless readability suffers.

2. Make `dadd(J+J)` branchless on point data.
   For `a4 == nullptr`, compute the complete `PointAdd(P, Q, ring)` result
   unconditionally, then select:
   - return `Q` when `P.z == 0`
   - return `P` when `Q.z == 0`

   The complete RCB formula already handles `P == Q` and `P == -Q`, so this
   removes the secret-dependent `ProjectiveEqual` branch and avoids computing
   an unnecessary double in the common `a4 == nullptr` path.

3. Keep the public `a4 != nullptr` branch.
   For the `a4` path, compute `PointDaddA(P, Q, *a4, ring)` and apply the same
   identity selects at the end. Before implementation, confirm whether
   `PointDaddA` needs additional equality handling beyond identity. If it does,
   compute the required alternative result and select without branching on
   point data.

4. Fix mixed add affine identity.
   In both mixed-add helpers, compute the existing formula result
   unconditionally, then select `P` when the affine input is the scalar BLST
   identity sentinel `(X == 0, Y == 0)`.

5. Make affine `cneg` constant-time.
   Replace the current `if (c)` branches in `G1AffinePoint::cneg` and
   `G2AffinePoint::cneg` with unconditional negation plus a mask/select back to
   the original value when `c == 0`.

### Tests

Add focused tests for both G1 and G2:

- `add(J+J)`:
  - generic points
  - `P == Q`
  - `P == -Q`
  - left identity
  - right identity
  - both identity
  - non-trivial `Z`
  - aliasing patterns

- `dadd(J+J)`:
  - same behavioral matrix as `add(J+J)`
  - include `a4 == nullptr`
  - include existing `a4 != nullptr` wrapper coverage if practical

- `add(J+A)` and `dadd(J+A)`:
  - generic points
  - projective left identity
  - doubling via affine input
  - inverse pair
  - affine identity sentinel `(0,0)`
  - non-trivial `Z`
  - aliasing patterns

- `from_affine(0,0)`:
  - all-zero affine converts to projective identity for both G1 and G2

- affine `cneg`:
  - `c == 0` leaves the point unchanged
  - `c == 1` negates `Y`
  - all-zero affine remains all-zero

### Validation

Portable AVX correctness:

```bash
ctest --test-dir out/build/x64-Debug-Avx -R "(Avx|CApi\\.AvxExpert)" --output-on-failure
```

If a fresh portable build is needed:

```bash
cmake -S . -B out/build/avx-portable -DCMAKE_BUILD_TYPE=Debug -DBLST_AVX512_PORTABLE=ON
cmake --build out/build/avx-portable --target blst_tests -j
ctest --test-dir out/build/avx-portable -R "(Avx|CApi\\.AvxExpert)" --output-on-failure
```

Native AVX validation should be run on an AVX-512 host before merging.

## Patch 2: MSM Sentinel Regression Coverage

### Scope

- `test/avx/Avx_msm_tests.cpp`
- C API tests if direct expert API coverage is useful

### Implementation Outline

No MSM algorithm change is planned. The curve mixed-add fix should make MSM
inherit scalar BLST's affine identity convention automatically.

Add regression tests that include at least one all-zero affine point with a
nonzero scalar and compare AVX Pippenger to scalar BLST or to a naive AVX sum
that also honors `from_affine(0,0)`.

### Validation

```bash
ctest --test-dir out/build/x64-Debug-Avx -R "(AvxMsm|CApi\\.AvxExpertP[12]Pippenger)" --output-on-failure
```

## Patch 3: GT Standard Narrowing Performance Fix

### Scope

- `blst/avx/fields/fp12.h`
- `blst/avx/curve/GT.h`
- GT and pairing/final-exponentiation tests
- AVX group-op benchmarks

### Implementation Outline

1. Add an elementwise narrowing helper to `FP12Ring`.
   It should convert `CommonElement<N>` to `StandardElement` by applying
   `Ring::prep_standard` to each of the 12 coordinates. This performs the
   Montgomery reduction work currently obtained as a side effect of
   `fp12_flat_mul(one, value, ring)`, without paying for a full Fp12 multiply.

2. Replace `to_standard(...)` call sites where the input already has the right
   Fp12 value and only needs bounds narrowing:
   - `GT::doubling`
   - `GT::mult_by_z`
   - `GT::mult_by_z_div_by_2`

3. Keep `GT::data` as `FP12::StandardElement`.
   Do not move to a wider internal representation in this patch.

### Tests

- Existing AVX GT wrapper tests.
- `test/avx/Avx_gt_glv4_tests.cpp`.
- Pairing/final-exponentiation tests that exercise `mult_by_z` and
  `mult_by_z_div_by_2`.
- Add a direct test that `GT::doubling` still matches scalar BLST.

### Benchmarks

Run benchmarks one at a time. Never run two benchmarks concurrently.

Primary check:

```bash
out/build/<preset>/benchmark/blst_benchmarks --benchmark_filter="Avx.*GT.*(Add|Dadd|Double)"
```

Expected result: `GT::doubling` should become materially cheaper than
`GT::add`, rather than about 1.5x slower.

Secondary check:

```bash
out/build/<preset>/benchmark/blst_benchmarks --benchmark_filter="Avx.*(Pairing|Final|GT)"
```

Confirm `GT::add`, `GT::mult`, and pairing paths are not regressed.

## Open Questions

- Which native AVX host and preset should be treated as the authoritative
  benchmark environment for the GT patch?

## Recommended Next Turn

Run a single serialized GT group-ops benchmark on the intended performance
machine, then review the generated diff for code style and benchmark impact.
