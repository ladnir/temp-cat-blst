# MatrixNoK Port Map

## Porting rule

The MatrixNoK update is a semantic port from VROOM commit `d6ecd534`, not a
source merge. Every selected operation is mapped onto CatBLST's existing
owners, contracts, vector abstraction, fixed constants, and tests. The current
AVX backend remains available as the comparison baseline until the candidate
passes every correctness and performance gate.

## Minimal dependency graph

| Order | Capability | VROOM source | CatBLST owner | Decision |
| ---: | --- | --- | --- | --- |
| 1 | Selected vector operations | `cpu/vector/avx512ifma.hpp` | `blst/backend/avx512/intrin_compat.*`, `blst/avx/cpu/vector/avx512ifma.hpp` | Translate only the low-128 horizontal/broadcast operations required by BLS12-381 MatrixNoK |
| 2 | Configurable input radix | `cpu/vector/{conversion,multiplication}.hpp` | Same CatBLST paths | Parameterize the integer-to-RNS word width; preserve 52 bits as the default |
| 3 | Conditional reduction | `src/elementwise_reduce.hpp` | `blst/avx/rns/elementwise_reduce.h` | Generalize the fixed conditional-subtraction schedule and state the actual coverage invariant |
| 4 | Fixed BLS12-381 parameters | `scripts/bls12_381_intrns4_params.hpp`, `src/bls12_381_*params.hpp` | Generated CatBLST constant source plus validation support | Import exact 8x50-bit constants with provenance; independently regenerate and validate them |
| 5 | Permutation change of base | `src/{rns_reduce_perm,bls12_381_perm,change_base}.hpp` | New narrow helper under `blst/avx/rns` | Translate the BLS12-381 8x8 path only; no BN254 or generic runtime parameter surface |
| 6 | MatrixNoK batch reduction | `src/{matrix_nok_batch,bls12_381_perm_batch}.hpp` | `blst/avx/rns/bounded_ring.h` plus a narrow helper | Keep one implementation, not the duplicate generic and BLS wrappers in VROOM |
| 7 | Conversion/inversion | `src/conversion_inversion.hpp` | `blst/avx/fields/conversion_inversion.h` | Port only the 50-bit packing and selected-ring changes |
| 8 | Field tower | `src/{fp2,fp12}.hpp` | `blst/avx/fields/{fp2,fp12}.h` | Port only changes required by the selected ring; preserve audited CatBLST contracts |
| 9 | Curve and pairing | `src/ec.hpp` and pairing callers | `blst/avx/{curve,pairing,gt}` | Review formula deltas independently; do not import generic API or domain changes |

## Upstream changes already present locally

- Variable permutation uses intrinsic order `(index, value)`.
- Native production operations use `_vv512_*` and map directly to `_mm512_*`;
  portable emulation remains isolated in `_cc512_*`.
- Fixed runtime constants replace mutable GMP-based production setup.
- Local field, curve, GT, pairing, validation, and C/Rust contracts take
  precedence over the research repository.

## Rejected or rewritten upstream material

- BN254 branches, parameters, and limb heuristics.
- GPU/CUDA paths.
- Runtime GMP construction, parameter search, and mutable initialization.
- Duplicate `batch_change_base_perm_qr`, `bls12_381_perm_batch`, and generic
  MatrixNoK entry points where one narrow BLS12-381 owner is sufficient.
- Broad scalar-operation and curve-template additions not reached by the
  BLS12-381 MatrixNoK call graph.
- Mask-zero IFMA and widening/narrowing vector helpers added by the same
  upstream commit but not reached by BLS12-381 MatrixNoK.
- Removal of CatBLST hardening, documentation, API contracts, or differential
  tests merely to make a source diff smaller.

## Reduction note

VROOM's generalized `min_sub_peel` comment says the schedule follows the set
bits of `upper_bound - stop`. The implementation instead tries every descending
power of two no larger than that gap. That schedule is arithmetically valid:
each step conditionally subtracts its multiple only when the lane is at least
that large, and the attempted multiples sum to at least the required gap. The
CatBLST version will document and test this coverage invariant directly rather
than repeat the misleading explanation.

## First implementation packet

1. Add the missing `_cc`/`_vv` low-128 extract/add/shift/broadcast operations
   with native-versus-emulated differential tests.
2. Remove the native `AVXVector::extract0` full-register memory round trip.
3. Parameterize integer-to-RNS conversion word width while retaining the
   current 52-bit default and behavior.
4. Generalize conditional subtraction with compile-time coverage assertions
   and focused boundary tests.

This packet does not select MatrixNoK, change field representation, import
constants, or alter a public API. It establishes the small, testable primitive
base needed by the candidate.

## Fixed-parameter packet

The BLS12-381 tables from pinned VROOM commit `d6ecd534` are imported into
`generated_matrixnok_constants.hpp` in a separate namespace. No production
owner references them yet, so this packet cannot change field representation
or dispatch.

`validate_matrixnok_constants.py` parses the C++ header and independently
regenerates every base-change matrix, correction vector, and fixed-point
quotient estimate. It also checks the CRT products, pairwise coprimality,
Montgomery inverses, QR premultiplication, and combined-ring rotation. Its
hostile self-test confirms representative one-bit corruptions are rejected.

The 16 RNS moduli are deliberately not required to be prime. The actual
algorithmic requirements are that they are distinct, pairwise coprime, odd,
and below the 50-bit bound.

## Gates before candidate arithmetic

- Combined native-plus-portable tests agree for every new operation.
- Portable-AVX Debug passes the directly affected vector and reduction tests.
- Native optimized code contains direct intrinsic lowering without conversion
  objects, full-register spills for lane extraction, or out-of-line `_vv`
  calls.
- The 52-bit incumbent produces unchanged focused results.
- Exact MatrixNoK constants and bound claims have independent executable
  validation before any candidate multiplication is enabled.
