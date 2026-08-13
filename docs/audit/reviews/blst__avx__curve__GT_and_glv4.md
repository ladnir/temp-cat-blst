# File Review: AVX GT and GLV4

## Record

| Field | Value |
| --- | --- |
| Paths | `blst/avx/curve/GT.h`; `blst/avx/gt/glv4.h` |
| Roadmap phase / queue | Phase 16 selected AVX GT closure |
| Reviewed commit | `4bb0c45a929ae1721a3fd9c942a758796e822b78` |
| Reviewer | Codex audit team |
| Date | 2026-08-13 |
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

## Contracts and invariants

- `mult` and `mult_vartime` require `nbits <= 256`. Scalar bytes are
  little-endian and may be null only for zero bits.
- Constant-time multiplication uses a fixed binary or 64-round GLV4 schedule.
  Its permitted timing distinction is whether noncanonical reduction was
  required. Variable-time paths accept public scalars only.
- Secret-capable scalar temporaries are cleared. Public variable-time storage
  is not wiped.
- The 16-entry GLV4 table, table scan, split bounds, schedule length, and
  unsigned and balanced decompositions remain unchanged.
- Public C GT boundaries are allocation-free and nonthrowing.
  `FinalExpScratch` is an explicit reusable C++ owner; growth can throw and one
  instance must not be used concurrently.
- Exact aliasing is supported only where the declaration states it. Batch
  final-exponentiation arrays must not overlap.

## Correctness and performance evidence

`AVX-P3-006`, `AVX-P3-007`, and `AVX-P4-001` accept the selected exponent and
Miller-related schedules. Compiler lowering, IFMA semantics, and scalar
observation bridges remain external.

The new AVX regression proves that unchecked decoding accepts a canonical
non-GT Fp12 value, membership rejects it, and checked decoding is failure
atomic. All 43 focused portable-Debug tests pass. The selected dispatch and
correspondence validators pass. Commit `4bb0c45a` changes no exponentiation,
Miller, or table arithmetic, so no benchmark ran.

## Findings and residual risk

`AUD-082`, `AUD-085`, `AUD-088`, `AUD-092`, `AUD-099`, `AUD-220`, `AUD-229`,
`AUD-241`, `AUD-277` through `AUD-282`, and `AUD-392`, `AUD-395`, `AUD-396`.

The unchecked decoder remains intentionally easy to misuse because API
compatibility precludes a distinct type. Precise documentation and negative
tests are the selected control. Native compiler/ISA, sanitizer/fuzz, and
statistical leakage evidence remain milestone gates.
