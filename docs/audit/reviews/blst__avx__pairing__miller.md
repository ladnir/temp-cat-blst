# File Review: Selected MatrixNoK Miller Engine

## Record

| Field | Value |
| --- | --- |
| Paths | `blst/avx/pairing/miller.h`; `matrixnok_pairing.h` |
| Roadmap phase / queue | Phase 16 selected AVX pairing closure |
| Reviewed commit | `4bb0c45a929ae1721a3fd9c942a758796e822b78` |
| Reviewer | Codex audit team |
| Date | 2026-08-13 |
| Disposition | `in_scope_production` |
| Profiles | native AVX-512/IFMA and portable developer Debug |
| Provenance | VROOM `d6ecd534`; Category Labs selected integration |

## Purpose and trust boundary

The Miller template is the finite internal pairing core. Its constructor
requires validated finite G1 and G2 subgroup affine points. The public GT
wrapper intercepts either infinity and returns the pairing identity before it
constructs this core.

`miller_loop` returns an Fp12 value in the Miller-result domain. The result is
not GT until final exponentiation completes. `PairingInverter` is
variable-time and accepts only nonzero public or intermediate Miller values.

## Contracts and invariants

- The direct path follows the fixed 68-operation BLS-x line schedule.
- The core performs no Release validation, allocation, or throwing operation.
- All state is instance-owned. There is no global cache, lock, or shared
  mutation.
- The formulas support finite subgroup inputs. Wider on-curve behavior is not
  part of the selected API contract.
- Final exponentiation owns its separate inversion, alias, and nonzero-input
  contracts.

## Correctness and performance evidence

`AVX-P4-001` proves the double and add geometry, sparse Fp12 application,
physical bounds, fixed schedule, and direct/prepared ledgers. Wrapper theorems
cover observed infinity and the prepared identity table. The concrete base
reciprocal and compiler/ISA lowering remain differential evidence boundaries.

The dead generic `pairing` helper was deleted. Tests and the historical
comparison benchmark now compose `MillerEngine` and final exponentiation
directly. The benchmark target compiles but was not executed. All 43 focused
portable-Debug tests and both AVX validators pass. No arithmetic body changed.

## Findings and residual risk

`AUD-220`, `AUD-233`, `AUD-241`, `AUD-392`, `AUD-394`, and `AUD-397`.

The historical 8,128-byte native Miller frame remains the current evidence.
Exact-hash native reexecution, sanitizer/fuzz coverage, and compiler/IFMA
lowering remain milestone evidence rather than incomplete source review.
