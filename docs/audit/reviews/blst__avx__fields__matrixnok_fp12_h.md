# `blst/avx/fields/matrixnok_fp12.h` Review

## Scope

This record covers the complete selected MatrixNoK Fp12 constants and type
binding through commit `ce22bc2`.

## Contracts

- `Fp12Constants` is a stateless compile-time policy for the exact selected RNS
  basis. It exposes generated offsets, fractions, and Frobenius constants.
- Every returned constant has the physical and represented-multiple bound
  required by its consumer in `FP12Ring`.
- `Fp12Ring` is the sole selected alias and is exactly
  `FP12Ring<FpRingCandidate, Fp12Constants>`.
- Construction has no dynamic initialization, allocation, runtime selection,
  or borrowed lifetime.

## Findings

- `AUD-268`: the alias retained a version selector for deleted schedules.
- `AUD-269`: duplicate candidate/selected aliases and inherited policy state
  obscured the one production owner.
- `AUD-273`: selection and generated-constant ownership were not documented at
  their final owner.

## Validation

- The selected production and focused test target compile after direct alias
  migration.
- Constant regeneration/self-test and the structural dispatch validator bind
  the policy to the generated tables and reject obsolete aliases or selectors.
- `AVX-P3-003`, `AVX-P3-004`, and `AVX-P4-002` are rebound to the exact source
  hash. The complete 21-row/61-hash correspondence ledger validates.
- The native optimized compile-only probe has no runtime policy initialization
  or selection symbol.

## Residual risk

Generated-table correctness and field semantics are owned by the referenced P1
through P3 proof rows. Compiler lowering and native IFMA semantics remain the
shared P0 boundary. No benchmark, remote execution, broad suite, CUDA,
BTx/BTE audit, or portable-AVX Release build ran for this packet.
