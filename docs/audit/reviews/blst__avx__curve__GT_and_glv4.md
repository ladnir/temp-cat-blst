# File Review: AVX GT and GLV4

## Record

| Field | Value |
| --- | --- |
| Paths | `blst/avx/curve/GT.h`; `blst/avx/gt/glv4.h` |
| Roadmap phase / queue | Phase 16 selected AVX GT closure |
| Reviewed commit | `e092f8d10c0be542ba95cd7bcef5d8b610a241b9` |
| Reviewer | Codex audit team |
| Date | 2026-08-28 |
| Disposition | `in_scope_production` |
| Profiles | native AVX-512/IFMA and portable developer Debug |
| Provenance | BLST target-group semantics; VROOM `d6ecd534`; Category Labs integration |

## Purpose and trust boundary

`GT.h` owns the selected native target-group API, the distinct Miller result,
prepared Miller lines, final exponentiation, and reusable C++ batch scratch.
`glv4.h` owns public and constant-time 256-bit exponent preparation and the
four-dimensional GT executors.

Ordinary GT arithmetic requires an existing GT element and performs no Release
membership validation. `deserialize_unchecked` is the deliberate exception:
it returns canonical general Fp12 storage. Callers may combine that candidate
only through the Fp12 API until `in_group` succeeds.

The selected implementation now owns the canonical codec and membership
predicates directly. Serialization is still the scalar CatBLST wire format.
Only the internal representation changed.

## Contracts and invariants

- `mult` and `mult_vartime` require `nbits <= 256`. Scalar bytes are
  little-endian and may be null only for zero bits.
- Constant-time multiplication uses a fixed binary or 64-round GLV4 schedule.
  Its permitted timing distinction is whether noncanonical reduction was
  required. Variable-time paths accept public scalars only.
- Secret-capable scalar temporaries are cleared. Public variable-time storage
  is not wiped.
- The 16-entry GLV4 table, table scan, split bounds, schedule length, and
  unsigned and balanced decompositions remain unchanged. Every table entry is
  assigned before the first possible read; aggregate zero-initialization is
  therefore unnecessary.
- Serialization maps all twelve native coefficients to the documented scalar
  coefficient order and emits canonical big-endian field elements. Decoding
  rejects every coefficient greater than or equal to the base-field modulus.
  Checked decoding additionally requires nonzero prime-order GT membership and
  is failure atomic.
- Public C GT boundaries are allocation-free and nonthrowing.
  `FinalExpScratch` is an explicit reusable C++ owner; growth can throw and one
  instance must not be used concurrently.
- Exact aliasing is supported only where the declaration states it. Batch
  final-exponentiation arrays must not overlap.

## Correctness and performance evidence

`AVX-P3-006`, `AVX-P3-007`, `AVX-P4-001`, and the GT owner in `AVX-P4-002`
are rebound to the exact current sources. The audit compared every intervening
change with the declarations named by those rows. The GLV4 change removes only
redundant initialization after all sixteen assignments were accounted for.
The GT changes add native codecs, equality and membership predicates, and cold
scratch preparation; they do not alter the mapped exponentiation, Miller-loop,
or final-exponentiation arithmetic. `GtExponentiation`, `GtCompiledExecutor`,
`GtCompiledPhysical`, and `Fp12CyclotomicPredicate` compile directly against
the current proof tree. Compiler lowering and IFMA semantics remain external.

The AVX regression proves that unchecked decoding accepts a canonical
non-GT Fp12 value, membership rejects it, and checked decoding is failure
atomic. On Peach, Clang 21.1.8 passes eight focused `AvxGt` tests, the selected
Fp12/GT exponentiation differential, and all five `GtDeserializeAudit` tests
against the exact current headers. The Rust native canonical GT round-trip also
passes when its test thread has a 16 MiB stack. The default Rust test-thread
stack is too small for the unoptimized inlined Miller loop; `AUD-403` records
that separate debug-build limitation. Existing Clang and GCC native Release
evidence, including the 32,184-byte optimized GT CT GLV4 frame, remains valid
for the unchanged mapped arithmetic.

## Findings and residual risk

`AUD-082`, `AUD-085`, `AUD-088`, `AUD-092`, `AUD-099`, `AUD-220`, `AUD-229`,
`AUD-241`, `AUD-277` through `AUD-282`, `AUD-392`, `AUD-395`, `AUD-396`,
`AUD-402`, and `AUD-403`.

The unchecked decoder remains intentionally easy to misuse because API
compatibility precludes a distinct type. Precise documentation and negative
tests are the selected control. Native compiler/ISA, sanitizer/fuzz, and
statistical leakage evidence remain milestone gates. Native exact-hash
execution is complete for the rebound GT and GLV4 owners. The current
repository-wide correspondence gate remains blocked only by the unaudited
`point.h` and `scalar_mult.h` mappings.
