# File Review: Selected MatrixNoK Curve Facades

## Record

| Field | Value |
| --- | --- |
| Paths | `blst/avx/curve/matrixnok_curve.h`; `matrixnok_g1.h`; `matrixnok_g2.h` |
| Roadmap phase / queue | Phase 16 selected AVX curve closure |
| Reviewed commit | `4bb0c45a929ae1721a3fd9c942a758796e822b78` |
| Reviewer | Codex audit team |
| Date | 2026-08-13 |
| Disposition | `in_scope_production` |
| Profiles | native AVX-512/IFMA and portable developer Debug |
| Provenance | VROOM `d6ecd534`; Category Labs integration and proof binding |

## Purpose and trust boundary

These files bind the selected MatrixNoK Fp and Fp2 policies to G1 and G2.
They own complete addition, mixed addition, doubling, negation, identity, and
the GLV/GLS policy aliases. They do not validate points.

The arithmetic accepts valid on-curve homogeneous points. Subgroup membership
is unnecessary at this layer. Affine infinity is the exact `(0,0)` sentinel.
Projective identity selection observes an exact zero RNS `z`; callers must use
the canonical carriers produced by the selected point boundary.

## Contracts and invariants

- `add_point`, `add_mixed_point`, and `double_point` preserve the selected
  curve and subgroup membership when their inputs have it.
- A null `a4` selects the a=0 complete-add schedule. A non-null `a4` is the
  coefficient of the curve containing both `dadd` inputs.
- The fixed schedules allocate nothing, throw nothing, and own no shared
  mutable state.
- Point-data choices use mask selection. The public choice of the `a4`
  schedule may branch before arithmetic begins.
- Return-by-value construction permits exact outer output/input aliasing.

## Correctness and performance evidence

`AVX-P4-003` proves the selected reduced-field group law, doubling, GLV/GLS
decomposition, recoding, and dispatch. The physical RNS producer bounds,
compiler lowering, and IFMA semantics remain separate evidence boundaries.
The nonzero-a4 schedule remains differential-only.

Commit `4bb0c45a` changes names, contracts, and provenance only. The focused
portable-Debug executable passes all 43 selected curve, GT, and pairing tests.
The dispatch and 61-hash correspondence validators pass after rebinding.
Peach Clang 21.1.8 then builds the exact `4bb0c45a` owners in native Release.
All 79 AVX tests pass, including complete-add, nonzero-a4, GLV/GLS, point,
MSM, GT, and pairing differentials. The optimized current-hash probe contains
no selected-vector conversion or wrapper symbol. GCC 15.2 also passes all 79
tests. Its optimized probe is byte-identical before and after the warning-only
`AUD-399` compiler-policy repair. No benchmark ran.

## Findings and residual risk

`AUD-215`, `AUD-218`, `AUD-241`, `AUD-392`, `AUD-393`, and `AUD-394`.

Compiler and IFMA correctness remain the common P0 boundary. Sanitizer/fuzz
coverage and statistical leakage remain milestone evidence. These gaps do not
reopen the completed source review.
