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

## Change-of-base packet

`matrixnok_change_base.h` owns one BLS12-381-only 8x8 r1 permutation kernel.
Its fixed eight-diagonal and fixed-batch loops are forced inline and unrolled,
using literal selected-vector rotations so Clang lowers them to immediate
`valignq`/`vshufi64x2` operations. It intentionally omits the conventional
quotient correction; the candidate ring may call it only after its
redundant-range preconditions are proved and tested. The generated row table
is explicitly 64-byte aligned and that layout is validator- and test-enforced.

Portable Debug tests compare both IFMA accumulators lane by lane against an
independent scalar 52-bit product split. Native Peach execution and optimized
code-generation evidence are recorded in `CHANGE_BASE_EVIDENCE.md`. The
kernel remains unreachable from production arithmetic in this packet.

## Gates before candidate arithmetic

- Combined native-plus-portable tests agree for every new operation.
- Portable-AVX Debug passes the directly affected vector and reduction tests.
- Native optimized code contains direct intrinsic lowering without conversion
  objects, full-register spills for lane extraction, or out-of-line `_vv`
  calls.
- The 52-bit incumbent produces unchanged focused results.
- Exact MatrixNoK constants and bound claims have independent executable
  validation before any candidate multiplication is enabled.

## Base-field candidate packet

`matrixnok_fp.h` now owns a side-by-side, non-dispatched BLS12-381
multiplication candidate. It uses register-only single-item conversion and
permutation paths while preserving the existing fixed-batch interfaces for
real batching. Its returning operation is fully inlined, while the
output-reference wrapper supports exact aliasing and delegates without a
second implementation. The candidate remains unreachable from production
field, curve, and pairing operations.

Portable Debug and Peach native AVX-512IFMA compare decoded results with the
portable scalar field implementation. Native code-generation and remaining
proof obligations are recorded in `BASE_FIELD_EVIDENCE.md`.

The optional focused performance target can compile the actual pinned VROOM
headers into a separate translation unit of the same benchmark executable.
It rejects every VROOM commit except the reviewed head, adds no production
dependency, and permits randomly interleaved same-process Clang comparisons.

## Fp2 candidate packet

`matrixnok_fp2.h` owns a side-by-side, non-dispatched Fp2 candidate. It uses
the ordinary four-product multiplication formula and a three-product square,
but reduces and changes base once per output rather than once per base-field
product. This preserves a small, explicit operation surface instead of making
the candidate imitate the incumbent delayed-expression template API.

Portable Debug and Peach native AVX-512IFMA compare every retained result with
scalar Fp2, including boundary values, exact aliases, and dependent chains.
Same-process Clang 21 measurements beat both the prior CatBLST 50-bit Fp2 and
pinned VROOM. Exact evidence is in `FP2_EVIDENCE.md`.

The candidate remains unreachable from production. The `AUD-214` selection
blocker was resolved at `19f37a1`: both ring owners are empty types using
compile-time-derived immutable fixed operations, with no hidden initialization
or allocation. Direct generic-helper differential tests, native code
generation, and serial remeasurement validate the redesign. Formal proof and
the remaining field-tower work still gate production selection.

## Fp12 and GT candidate packet

`matrixnok_fp12.h` supplies a stateless constants policy to the existing Fp12
arithmetic template. Fixed offsets and Frobenius values are bound to the exact
MatrixNoK parameter set rather than selected by limb width. This prevents the
incumbent and candidate 50-bit rings from silently sharing incompatible RNS
constants (`AUD-216`). The generic template keeps its existing default policy,
so production dispatch and source compatibility are unchanged.

Focused native tests compare multiplication, squaring, Frobenius maps,
conjugation, cyclotomic squaring, and CT/VT GT exponentiation with scalar
CatBLST. The core multiplication matches pinned VROOM's IFMA counts and code
size after removing harmful forced unrolling from the fixed outer R2 loop
(`AUD-217`). Same-process measurements put candidate Fp12 multiplication at
303 ns versus 362 ns for the incumbent and 304 ns for VROOM. Exact evidence,
rejected alternatives, residual cyclotomic/FFT gaps, and source hashes are in
`FP12_EVIDENCE.md`.

Follow-up commit `3e1f99a` validates the existing sparse line-multiplication
kernel, including exact aliasing, over MatrixNoK. It also generalizes the hard
final-exponentiation helpers over the compile-time Fp12 policy and matches
scalar CatBLST on three independent inputs. The production default policy
still compiles unchanged. A MatrixNoK inversion boundary is deliberately not
invented in this packet: the test imports an independently computed inverse.

The candidate remains isolated. MatrixNoK inversion, G2 performance, Miller
and pairing integration, formal proof, and release evidence still gate
selection.

## G2 candidate packet

`matrixnok_g2.h` applies the existing GLS executor to a production-shaped,
non-dispatched MatrixNoK G2 curve policy. G1 and G2 share the exact complete
homogeneous add schedule in `matrixnok_curve.h`; Fp2 delegates signed
small-basis expansion component-wise to the base ring. This centralization
prevents the `AUD-215` repair from silently diverging across curve groups.

Portable Debug compilation and Peach native AVX-512IFMA execution pass the
focused complete-operation and CT/VT GLS scalar differentials. The complete
add is fully inlined with no call or vector conversion boundary. Optimized CT
selection retains fixed control flow and a full table scan.

One serial same-process gate places candidate add, double, CT GLS, and VT GLS
28.9%, 15.1%, 21.6%, and 23.9% ahead of the incumbent. The candidate is faster
than pinned VROOM except for a 0.9% double delta. Performance selection passes.
`AUD-218` records the inherited greater-than-102-KiB peak CT GLS stack use for
a separate direct-table-construction redesign. Exact evidence is in
`G2_EVIDENCE.md`.

## FFT preview packet

A benchmark-only adapter compiles the complete pinned VROOM Fp12 stack in the
existing focused benchmark executable while preserving CatBLST's roots and
GLV4 schedule. A 256-element GT forward/inverse round trip passes. Native
Clang 21 inverse GT FFT latency is 15.8% to 16.9% below the incumbent from 256
through 16,384 elements, justifying the Fp12/GT port. Exact construction and
measurements are in `FFT_PREVIEW_EVIDENCE.md`.

The corrected G1 preview imports CatBLST's canonical generator rather than a
copied literal and independently verifies the homogeneous curve equation. It
passes the two-element butterfly and a 256-element forward/inverse round trip.
Native Clang 21 forward G1 FFT latency is 9.3% to 11.0% below the incumbent
from 256 through 16,384 elements. VROOM therefore establishes a useful G1
performance target, while the port-map decision remains unchanged: adapt the
MatrixNoK arithmetic to CatBLST's audited curve owners and contracts rather
than import `ec.hpp` wholesale.
