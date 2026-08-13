# MatrixNoK G1 Historical Selection Evidence

## Packet boundary

- Historical candidate implementation: `2fc0ae98f962f5ab2259ed2d368f6328713b7310`
- Original G1 selection implementation: `7e731d999780be0723c2335a4fc44ee54afd18aa`
- Upstream algorithm commit: `d6ecd534c54a41702ea80ae75b52250b6d51b372`
- Compared VROOM head: `26bd9aee2168f167583cf096add355547a7b56d4`
- C++ curve owner: `blst/avx/curve/matrixnok_g1.h`
- RNS/field owner: `blst/avx/rns/matrixnok_fp.h`

This document records the isolated candidate packet that preceded production
selection. MatrixNoK now owns the production AVX G1 path. The historical packet
did not alter public types, symbols, or dispatch. It evaluated CatBLST's audited
complete homogeneous formulas and GLV executors over the fixed, stateless
MatrixNoK base-field owner.

## Correctness design

The candidate converts only at the test boundary. Its independent oracle is
scalar CatBLST, not another invocation of the MatrixNoK implementation. The
focused native suite covers:

- zero, one, and fixed MatrixNoK constants;
- signed delayed field products and normalized cross terms;
- generator conversion and homogeneous-to-Jacobian conversion without an
  inversion;
- complete add, double, equal inputs, opposite inputs, both infinity operand
  positions, and a 32-step deterministic addition chain;
- constant- and variable-time GLV multiplication on the generator and doubled
  generator for zero, one, two, `r-1`, fixed dense values, and 16 deterministic
  dense canonical scalars; and
- a 16-element forward/inverse G1 FFT round trip against the original scalar
  points.

All 11 `AvxMatrixNoKPrimitives.*` tests pass under Clang 21.1.8 on native
AVX-512IFMA. The independent fixed-constant validator also passes, including
the newly retained RNS encodings of one, the field modulus, and signed-range
offsets.

The first delayed-expansion add translation exposed a correctness hazard:
expanding a negative small-m2 cross term without first lifting its represented
integer lost the RNS multiple. `FpRingCandidate::prep_expand` now adds the
smallest sufficient power-of-two multiple of `p` for a statically negative
RNS bound, then performs the fixed reduction and expansion. This changes no
field value and preserves the exact integer needed for base extension.

## Native code generation

Clang 21.1.8 compiled the G1 probes with `-O3 -DNDEBUG -march=native`.
The original selected complete-add probe contained no call, selected-vector
wrapper, conversion object, or lane-memory round trip. It contained 447 IFMA
instructions and a 600-byte logical stack adjustment with compiler spills.
After the G2 packet centralized the signed-expansion-safe schedule, the same
G1 probe remains fully inlined with no call; its current body is `0x188a`
bytes and contains 363 IFMA instructions. This later source state retains the
native scalar differential pass but will be remeasured with the G2 gate.

The generic MatrixNoK point-double body contains 347 IFMA instructions, no
internal calls, and no stack adjustment. Its explicit audit wrapper calls that
one template body and uses a 512-byte return adapter. A candidate-specific
delayed-expansion rewrite was correct but measured identically at 138 ns, so it
was rejected rather than retaining duplicate source.

The complete-add spill count and primitive delta remain optimization targets;
they are not hidden by the end-to-end selection result.

## Native performance gate

One Google Benchmark process ran at a time on Peach's AMD Ryzen 9 7950X. The
same Clang 21.1.8 Release binary contained CatBLST's incumbent, this candidate,
and actual VROOM commit `26bd9aee`. Each G1 FFT case ran five repetitions.
Google Benchmark reported CPU scaling enabled, so these are selection results,
not final release claims.

| Elements | Incumbent | CatBLST MatrixNoK | VROOM | vs incumbent | vs VROOM |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 256 | 20.2 ms | 18.3 ms | 18.6 ms | -9.4% | -1.6% |
| 1,024 | 111 ms | 100 ms | 99.9 ms | -9.9% | +0.1% |
| 4,096 | 562 ms | 508 ms | 503 ms | -9.6% | +1.0% |
| 16,384 | 2,724 ms | 2,463 ms | 2,421 ms | -9.6% | +1.7% |

The candidate is faster than the incumbent at every size and remains inside
the project's 2% primary-case tolerance relative to VROOM. Primitive probes
remain slower: complete add is 152 ns versus VROOM's 135 ns (+12.6%), and
double is 138 ns versus 124 ns (+11.3%). The FFT result is the selection gate;
the primitive figures identify where later compiler/layout work may still pay.

## Rejected alternatives

- direct Booth-table variable-time GLV was materially slower than wNAF;
- wNAF widths 4 and 6 lost to the retained width 5;
- preconverting FFT roots did not materially change the result;
- forcing complete add, double, or the fixed batch reducer out of line worsened
  latency; and
- a register-extraction rewrite of the fixed batch matrix loop did not improve
  the target arities.

No losing benchmark-only arithmetic variant remains in the source.

## Exact correspondence inputs

- `matrixnok_curve.h` SHA-256:
  `9fec0e6abbdf727011a7b4f859f0ae551235ea9aee9d9f58c49ff61e12d13556`
- `matrixnok_g1.h` SHA-256:
  `cb0651d44a8390a624d846f78e51fc24f94afa98b349b6a6599f54bc18c5dfe6`
- `matrixnok_fp.h` SHA-256:
  `36143f479eee5460595bf632b5a9e803388a6549c505d9939000e7c853bfeadf`
- `matrixnok_fixed.h` SHA-256:
  `73ba31d8a347208ac84173790a61929a34f443a0d47783cffa21559661ac026f`
- `generated_matrixnok_constants.hpp` SHA-256 at the measured selection commit:
  `d1bd0045e8d7f223cb332842bd4c9023daefc715079c0649fbce953cabdfea43`

The current repaired constant hash is recorded in
`proof/avx/CORRESPONDENCE.csv`; the eight ConvertTo estimate changes do not
alter G1 arithmetic or this historical performance measurement.
- `scalar_mult.h` SHA-256:
  `b1a3a30196a5e1a24bcaddf9963791aba9115fdfe5ecc1abcef1b4c9d8004396`

## Residual gates

- Production dispatch remains unchanged until the retained Fp12, G2, Miller,
  final-exponentiation, and pairing slices pass their own gates.
- The G1/RNS arithmetic and all source-level bound transitions remain formal
  proof obligations.
- The complete-add spill pattern may be revisited after the first proof
  milestone, when source transformations can be checked against the model.
- Sanitizers, fuzzing, broad suites, Apple CI, CUDA, BTx/BTE, and portable-AVX
  Release were outside this focused packet.
