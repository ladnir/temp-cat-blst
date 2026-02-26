# AVXBLS Additive Phased Implementation Plan (Separated by Type)

## 1) Objective

Port AVXBLS into `blst` in additive phases only.

Additive-only means:

- add new AVX files and tests,
- do not change existing scalar/CUDA behavior,
- do not re-route existing public APIs to AVX paths,
- do not remove existing code.

Primary development support:

- `blst/backend/avx512/intrin_compat.h` (already added), used with `BLST_AVX512_PORTABLE` for non-AVX machines.

---

## 2) Required File Separation

All new AVX functions are separated by algebraic type.

Field files:

- `blst/fields/Fp_avx.h`
- `blst/fields/Fp_avx.cpp`
- `blst/fields/Fp2_avx.h`
- `blst/fields/Fp2_avx.cpp`
- `blst/fields/Fp4_avx.h`
- `blst/fields/Fp4_avx.cpp`
- `blst/fields/Fp6_avx.h`
- `blst/fields/Fp6_avx.cpp`
- `blst/fields/Fp12_avx.h`
- `blst/fields/Fp12_avx.cpp`

Curve files:

- `blst/curve/G2_avx.h`
- `blst/curve/G2_avx.cpp`

Common shared AVX header (minimal primitives/constants only):

- `blst/fields/AvxCommon.h`

Rules:

1. `Fp*` arithmetic functions live only in `blst/fields/Fp_avx.*`.
2. `Fp2*` arithmetic functions live only in `blst/fields/Fp2_avx.*`.
3. `Fp4*` arithmetic functions live only in `blst/fields/Fp4_avx.*`.
4. `Fp6*` arithmetic functions live only in `blst/fields/Fp6_avx.*`.
5. `Fp12*` arithmetic functions live only in `blst/fields/Fp12_avx.*`.
6. `G2` AVX logic lives in `blst/curve/G2_avx.*`.
7. AVX data layout types are owned by their associated header:
   - `fp_*` and `fpx2_*` types in `Fp_avx.h`
   - `fp2_*` and `fp2x2_*` types in `Fp2_avx.h`
   - `fp4_*` and `fp4x2_*` types in `Fp4_avx.h`
   - Fp6 packed aliases used by Fp6 functions in `Fp6_avx.h`
   - Fp12 packed aliases used by Fp12 functions in `Fp12_avx.h`
8. `AvxCommon.h` is limited to shared SIMD wrappers, compile-time constants, and cross-file utility macros.
9. Higher-level headers include lower-level headers (`Fp2_avx.h` includes `Fp_avx.h`, `Fp4_avx.h` includes `Fp2_avx.h`, etc.).

---

## 3) Build Modes and Gating

Mode A:

- `BLST_AVX512_PORTABLE=ON` (required for all phases)

Mode B:

- `BLST_AVX512=ON` (run when AVX512 hardware exists)

Phase gate:

1. mode A build + tests pass,
2. mode B build + tests pass (or compile-only in CI when hardware missing),
3. only then move to next phase.

---

## 4) Test Strategy (Applies to Every Phase)

Every implemented function must have direct or traced coverage.

Function coverage requirement:

1. each public AVX function has direct unit tests,
2. each internal helper has either direct tests or direct validation through a deterministic harness that isolates its effect,
3. if helper cannot be isolated through public path, expose a test-only wrapper under `#ifdef BLST_TESTING`.

Reference oracle order:

1. existing non-AVX scalar implementation in `blst` (preferred),
2. simple scalar reference routine in test file (for layout/permute/blend or widened intermediates),
3. dual-check when available (both scalar API and independent math helper).

Comparison representation contract (mandatory):

1. scalar oracle values are in existing `blst` 64-bit-limb Montgomery representation,
2. AVX paths that use internal radix-48 packing must be converted back to canonical 64-bit Montgomery form before equality checks,
3. widened intermediates are compared in widened form only when function output is intentionally unreduced,
4. reduced kernels are always compared after canonical Montgomery normalization.

Test categories required for each function:

1. deterministic vectors (fixed edge patterns),
2. randomized differential tests (fixed seed),
3. algebraic identities/invariants when applicable,
4. cross-layout equivalence when multiple packed layouts exist.

Schedule semantics policy (critical for Fp6/Fp12 and sparse line kernels):

1. do not assume every output lane corresponds to an independent scalar multiplication.
2. many AVXBLS kernels are schedule kernels that intentionally mix lanes with `perm_*`, `blend_*`, and `shuf_*` before and after arithmetic.
3. only use direct lane-by-lane scalar assertions when the lane contract is explicitly documented and proven in tests.
4. when lane mapping is not yet fully documented, use schedule-agnostic differential tests:
   - build broadcast/degenerate packs (same scalar value in all schedule slots),
   - run AVX kernel and scalar oracle,
   - assert scalar oracle results are present in reduced outputs.
5. once lane contracts are derived, add strict slot-mapping tests in addition to schedule-agnostic tests.

Per-function test contract (mandatory):

1. for every function listed in each phase, add at least one dedicated test block named after that function,
2. each function test block must include:
   - deterministic vector compare,
   - randomized compare,
   - edge-case compare,
3. comparison target must be scalar/non-AVX reference when available,
4. if no scalar API exists, compare against a local scalar oracle in test code implementing the same math/permute semantics.

Edge-case catalog to reuse across phases:

- zero, one, max limb, alternating limb masks,
- near-modulus values (`p-2`, `p-1`, `p`, `p+1` where representable pre-reduction),
- carry/borrow stress patterns,
- single-bit set in each limb/lane,
- mixed-sign behavior for arithmetic shift helpers.

Fast-but-robust default volume:

- deterministic: at least 32 vectors per function,
- random differential: at least 200 vectors per function,
- stress subset: at least 1,000 vectors per phase across core kernels.

Regression policy:

1. every bug found during implementation adds a permanent regression vector,
2. regression vectors are tagged with bug id/commit reference,
3. when helper semantics are clarified/changed, add targeted new tests immediately.

---

## 5) Test File Layout

Approved test files:

- `test/avx512_phase0_layout_tests.cpp`
- `test/avx512_fp_tests.cpp`
- `test/avx512_fp2_tests.cpp`
- `test/avx512_fp4_tests.cpp`
- `test/avx512_fp6_tests.cpp`
- `test/avx512_fp12_tests.cpp`
- `test/avx512_g2_line_tests.cpp`
---

## 6) Phase 0: Layout/Infrastructure (Before Fp Arithmetic)

Implementation files:

- `blst/fields/AvxCommon.h`
- `blst/fields/Fp_avx.h` (initial Fp layout type declarations)
- test-only wrappers as needed in corresponding `*_avx.h/.cpp`

Functions to add now (classified under Fp support, implemented in `Fp_avx.*`):

- `conv_64to48_fp`
- `conv_64to48_fp_8x1w`
- `conv_64to48_fp_4x2w`
- `conv_64to48_mpi`
- `conv_48to64_fp`
- `conv_48to64_fp_8x1w`
- `conv_48to64_fp_4x2w`
- `conv_48to64_mpi`
- `carryp_mpi48`
- `carryp_fp_4x2w`
- `get_channel_8x1w`
- `get_channel_4x2w`

Required tests (`avx512_phase0_layout_tests.cpp`):

1. `64->48->64` round-trip for scalar Fp words.
2. `64->48->64` round-trip for packed `8x1w` and `4x2w`.
3. lane order checks for `VSET`-equivalent packing.
4. channel extraction correctness for all channels.
5. carry normalization checks (`< 2^48`) after `carryp_*`.
6. negative tests for malformed input assumptions in debug mode (if assertions exist).

Gate:

- all phase-0 tests pass in modes A/B.

---

## 7) Phase 1: Fp

Implementation files:

- `blst/fields/Fp_avx.h`
- `blst/fields/Fp_avx.cpp`

Types:

- `fp_8x1w`
- `fpx2_8x1w`
- `fp_4x2w`
- `fpx2_4x2w`

Functions to add:

- `add_fp_8x1w`
- `sub_fp_8x1w`
- `asx4_fp_8x1w`
- `asax6_fp_8x1w`
- `assx6_fp_8x1w`
- `add_fp_4x2w`
- `sub_fp_4x2w`
- `asx2_fp_4x2w`
- `add_fpx2_8x1w`
- `sub_fpx2_8x1w`
- `aaasx2_fpx2_8x1w`
- `aassx2_fpx2_8x1w`
- `asx4_fpx2_8x1w`
- `ssaax2_fpx2_8x1w`
- `ssasx2_fpx2_8x1w`
- `add_fpx2_4x2w`
- `sub_fpx2_4x2w`
- `asx2_fpx2_4x2w`
- `aass_fpx2_4x2w`
- `mul_fpx2_8x1w` (implemented from AVXBLS `mul_fpx2_8x1w_v2`)
- `mul_fpx2_8x1w_v2`
- `mul_fpx2_4x2w`
- `mul_fpx2_4x2w_v1`
- `redc_fpx2_8x1w`
- `redc_fpx2_4x2w`

Required tests (`avx512_fp_tests.cpp`):

1. add/sub tests against scalar Fp addition/subtraction lane-by-lane.
2. linear-combination helper tests (`as*`, `aa*`, `ss*`) against direct scalar formulas.
3. widened multiply tests for each `mul_fpx2_*` kernel vs scalar widened multiply oracle.
4. reduction tests (`redc_*`) vs scalar Montgomery reduction.
5. layout equivalence tests: `8x1w` and `4x2w` produce identical canonical scalar outputs.
6. alias equivalence tests: `mul_fpx2_8x1w == mul_fpx2_8x1w_v2` and `mul_fpx2_4x2w == mul_fpx2_4x2w_v1`.
7. identity/property tests: commutativity (mul/add), subtraction inverse checks.
8. edge-case modulus boundary tests plus randomized differential tests.

Gate:

- all Fp tests pass in modes A/B.

---

## 8) Phase 2: Fp2

Implementation files:

- `blst/fields/Fp2_avx.h`
- `blst/fields/Fp2_avx.cpp`

Functions to add:

- `add_fp2_8x1x1w`
- `add_fp2_4x2x1w`
- `sub_fp2_4x2x1w`
- `aass_fp2_4x2x1w`
- `assa_fp2_4x2x1w`
- `sass_fp2_4x2x1w`
- `ssaa_fp2_4x2x1w`
- `mul_by_u_plus_1_fp2_4x2x1w`
- `mul_fp2_2x4x1w`
- `add_fp2_2x2x2w`
- `sub_fp2_2x2x2w`
- `as_fp2_2x2x2w`
- `sa_fp2_2x2x2w`
- `mul_by_u_plus_1_fp2_2x2x2w`
- `mul_fp2_1x4x2w`
- `add_fp2x2_8x1x1w`
- `sub_fp2x2_8x1x1w`
- `asx4_fp2x2_8x1x1w`
- `mul_by_u_plus_1_fp2x2_8x1x1w`
- `add_fp2x2_4x2x1w`
- `sub_fp2x2_4x2x1w`
- `sax2_fp2x2_4x2x1w`
- `mul_by_u_plus_1_fp2x2_4x2x1w`
- `redc_fp2x2_4x2x1w`
- `redc_fp2x2_2x2x2w`
- `add_fp2x2_2x2x2w`
- `sub_fp2x2_2x2x2w`
- `mul_by_u_plus_1_fp2x2_2x2x2w`
- `mul_fp2_4x2x1w`
- `sqr_fp2_4x2x1w`
- `mul_fp2_2x2x2w`
- `sqr_fp2_2x2x2w`
- `mul_fp2x2_8x1x1w`
- `mul_fp2x2_4x2x1w`
- `sqr_fp2x2_4x2x1w`
- `mul_fp2x2_2x4x1w`
- `mul_fp2x2_2x2x2w`
- `sqr_fp2x2_2x2x2w`
- `mul_fp2x2_1x4x2w`

Required tests (`avx512_fp2_tests.cpp`):

1. each add/sub/combination helper validated against scalar `Fp2` formulas.
2. each `mul_by_u_plus_1_*` validated against scalar multiplication by `(u+1)`.
3. each `mul_fp2_*` and `sqr_fp2_*` validated against scalar `Fp2::mul/sqr`.
4. each widened `mul_fp2x2_*` validated against widened scalar reference and reduction path checks.
5. each reduction helper validated against scalar reduced result.
6. cross-layout equivalence (`8x1x1w`, `4x2x1w`, `2x2x2w`, `1x4x2w`) for same input batch.
7. conjugation-compatible sanity checks where applicable.
8. randomized and edge-pattern vectors across both components (`c0`, `c1`).

Gate:

- all Fp2 tests pass in modes A/B.

---

## 9) Phase 3: Fp4

Implementation files:

- `blst/fields/Fp4_avx.h`
- `blst/fields/Fp4_avx.cpp`

Functions to add:

- `Fp4Pack2x2x2x1::sqr_v1`
- `Fp4Pack1x2x2x2::sqr_v1`
- `Fp4Pack2x2x2x1::sqr_v2`
- `Fp4Pack2x2x2x1::sqr_v3`
- `Fp4Pack1x2x2x2::sqr_v2`

Required tests (`avx512_fp4_tests.cpp`):

1. each square variant vs scalar Fp4 square oracle.
2. variant consistency checks (`v1/v2/v3`) for same layout/input.
3. cross-layout consistency (`2x2x2x1w` vs `1x2x2x2w`) after unpack.
4. edge tests with sparse and dense coefficients.
5. randomized differential tests.

Gate:

- all Fp4 tests pass in modes A/B.

---

## 10) Phase 4: Fp6

Implementation files:

- `blst/fields/Fp6_avx.h`
- `blst/fields/Fp6_avx.cpp`

Functions to add:

- `mul_fp6x2_4x2x1x1w`
- `mul_fp6x2_2x2x2x1w`
- `mul_fp6x2_1x2x2x2w`
- `mul_by_xy00z0_fp6x2_2x2x2x1w`
- `mul_by_xy00z0_fp6x2_2x4x1x1w`

Required tests (`avx512_fp6_tests.cpp`):

1. each multiply function vs scalar Fp6 multiply oracle.
2. sparse `xy00z0` functions vs scalar sparse multiply oracle.
3. cross-layout consistency checks between supported layouts.
4. edge/random vectors with independent stress on each Fp2 limb.
5. widened-intermediate sanity (non-reduced ranges before reduction points).
6. schedule-agnostic broadcast tests for schedule kernels (`2x2x2x1w`, `4x2x1x1w`, sparse variants), to avoid false negatives from undocumented lane mixing.
7. strict lane-position checks only for layouts with fully documented slot contracts.

Gate:

- all Fp6 tests pass in modes A/B.

---

## 11) Phase 5: Fp12

Implementation files:

- `blst/fields/Fp12_avx.h`
- `blst/fields/Fp12_avx.cpp`

Functions to add:

- `cyclotomic_sqr_fp12_vec_v1`
- `cyclotomic_sqr_fp12_vec_v2`
- `compressed_cyclotomic_sqr_fp12_vec_v1`
- `mul_fp12_vec_v1`
- `mul_fp12_vec_v2`
- `mul_fp12_vec_v3`
- `mul_fp12_vec_v4`
- `mul_by_xy00z0_fp12_vec_v1`
- `mul_by_xy00z0_fp12_vec_v2`
- `sqr_fp12_vec_v1`
- `mul_fp12_vector`
- `mul_by_xy00z0_fp12_vector`
- `sqr_fp12_vector`
- `cyclotomic_sqr_fp12_vector`

Required tests (`avx512_fp12_tests.cpp`):

1. each multiply/square variant vs scalar `Fp12` oracle.
2. variant consistency (`mul v1..v4`) for identical inputs.
3. cyclotomic square vs scalar cyclotomic square for cyclotomic-valid inputs.
4. compressed cyclotomic square validated via decompress/expand consistency check.
5. sparse `xy00z0` multiply vs scalar sparse path.
6. algebraic checks (`(a*b)*c == a*(b*c)` sample-based, `a*1=a`, `a^2=a*a`).
7. edge/random vectors across all 12 Fp coefficients.
8. wrapper-vs-scalar checks for:
   - `mul_fp12_vector` vs `Fp12::mul`,
   - `mul_by_xy00z0_fp12_vector` vs `Fp12::mul_by_xy00z0`,
   - `sqr_fp12_vector` vs `Fp12::sqr`,
   - `cyclotomic_sqr_fp12_vector` vs `Fp12::cyclotomic_sqr` on cyclotomic inputs.

Gate:

- all Fp12 tests pass in modes A/B.

---

## 12) Phase 6: G2 / Pairing Line Kernels

Implementation files:

- `blst/curve/G2_avx.h`
- `blst/curve/G2_avx.cpp`

Functions to add:

- `line_dbl_vec_v1`
- `line_dbl_vec_v2`
- `line_add_vec_v1`
- `line_by_Px2_4x2x1w`
- `line_by_Px2_2x2x2w`
- `line_dbl_vector`
- `line_add_vector`

Required tests (`avx512_g2_line_tests.cpp`):

1. `line_dbl_vec_v1` and `line_dbl_vec_v2` vs scalar line-function oracle.
2. `line_dbl` variant consistency (`v1` vs `v2`).
3. `line_by_Px2_*` coefficient checks against scalar expansion.
4. deterministic Miller-step traces comparing AVX and scalar intermediate states.
5. randomized projective input sets (valid points only) with strict equality checks.
6. `line_add_vec_v1` tests must use non-exceptional mixed-add tuples only:
   - `R=(X1,Y1,Z1)` and affine `Q=(X2,Y2)` are on-curve,
   - `H = X2*Z1^2 - X1 != 0`,
   - avoid synthetic random-limb tuples as scalar-oracle inputs for this kernel.
7. wrapper-vs-scalar checks:
   - `line_dbl_vector` vs scalar `line_dbl`,
   - `line_add_vector` vs scalar `line_add`.

Phase-6 note from implementation:

- `line_add_vec_v1` follows the same non-exceptional `madd-2007-bl` contract as avxbls.
- Random Montgomery limb tuples are not a valid oracle domain for this kernel; use points derived from
  real group elements (for example, scalar multiples of `G2Point::generator()`), then compare AVX
  and scalar outputs directly.
- Important packing detail: for this kernel, `X1Y1` is packed as `X1 | Y1` (not `Y1 | X1`).

Gate:

- all G2 line tests pass in modes A/B.

---

## 13) Internal Helper Testability Policy

To satisfy "each function tested":

1. expose internal helpers via `#ifdef BLST_TESTING` declarations in corresponding `*_avx.h`,
2. keep production visibility unchanged outside test builds,
3. unit tests must call helper directly at least once with deterministic and random vectors.

For permutation/blend/shuffle helpers:

- reference oracle is a small scalar lane-rewrite function in test code,
- each helper gets full-lane expected-output checks.

---

## 14) Completion Criteria for This Additive Track

1. all separated AVX files exist by type (`Fp`, `Fp2`, `Fp4`, `Fp6`, `Fp12`, `G2`),
2. every implemented function has explicit test coverage,
3. scalar/non-AVX oracle comparisons are present for all arithmetic kernels,
4. edge-case and randomized differential coverage is present in each phase test file,
5. AVXBLS wrapper APIs are present and tested:
   - `mul_fp12_vector`, `mul_by_xy00z0_fp12_vector`, `sqr_fp12_vector`, `cyclotomic_sqr_fp12_vector`,
   - `line_dbl_vector`, `line_add_vector`,
6. no integration routing into existing public API paths has been introduced.

---

## 15) Test Expansion Guidance During Implementation

The baseline tests above are required minimums.

During implementation, add more tests whenever:

1. a bug is found,
2. a layout invariant is clarified,
3. a function precondition becomes stricter,
4. a new helper is introduced,
5. behavior differs between `BLST_AVX512_PORTABLE` and native AVX512.

This plan intentionally encourages growing tests as implementation findings emerge.

---

## 16) AVXBLS Cycle Benchmark Mode

Purpose:

- provide a direct in-tree cycle benchmark with avxbls-style sections and labels,
- measure the ported AVX kernels (and matching scalar baselines where available),
- compare numbers against avxbls runs without changing the main benchmark flow.

Runtime entry:

- `blst_benchmarks --user --avxbls`

Optional tuning flags:

- `--iters N`
- `--warmup N`

Implemented sections:

1. `TIMING - PAIRING`
2. `TIMING - FP`
3. `TIMING - FP2`
4. `TIMING - FP4`
5. `TIMING - FP6`
6. `TIMING - FP12`
7. `TIMING - LINE`

FP-kernel version note:

- this port intentionally keeps `mul_fpx2_8x1w_v2` and `mul_fpx2_4x2w_v1` only (no v1/v3/v4 FP multiply variants).

Success criteria:

1. `--user --avxbls` runs and prints all section headers above.
2. each line emits `#cycle = ...` for the mapped operation.
3. output labels stay stable so results can be diffed across runs/revisions.
