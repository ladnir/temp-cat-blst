# MatrixNoK Fixed-Parameter Evidence

## Scope

This record covers only the isolated BLS12-381 MatrixNoK parameter import. The
tables are not connected to production arithmetic or dispatch in this packet.

- Upstream commit: `d6ecd534c54a41702ea80ae75b52250b6d51b372`
- Upstream path: `scripts/bls12_381_intrns4_params.hpp`
- Upstream Git blob: `d0e0437007f4dec6f193d528ce618f8b516e4fbe`
- CatBLST path: `blst/avx/generated_matrixnok_constants.hpp`
- Imported-header SHA-256:
  `d1bd0045e8d7f223cb332842bd4c9023daefc715079c0649fbce953cabdfea43`
- Namespace: `blst::avx::matrixnok::generated`

The data retains VROOM's MIT license and Simon Langowski attribution. The C++
translation uses immutable inline constants and requires no runtime GMP,
initialization, search, allocation, or mutable global state. Cache-line
alignment annotations added by CatBLST change no numeric data.

## Independent validation

`scripts/audit/validate_matrixnok_constants.py` uses only Python's standard
library and parses the checked-in C++ header directly. It does not import the
VROOM generator. It verifies:

- the 8-lane dimensions, 50-bit RNS bound, and 50-bit integer input radix;
- distinct odd pairwise-coprime moduli and exact `M1`, `M2`, and `M1M2`
  products;
- the required CRT, target-field inverse, QR premultiplier, and total-rotation
  congruences;
- all 16 Montgomery inverse factors modulo `2^52`;
- independent regeneration of every r1, r2, convert-to, and convert-from
  matrix, correction vector, and 64-bit fixed-point quotient estimate; and
- exact equivalence of the cyclic r1 permutation table to the base matrix; and
- 64-byte alignment declarations for every generated hot matrix, correction,
  and modulus array.

Primality is intentionally not claimed or required. These are RNS moduli; the
used algorithms require pairwise coprimality and oddness, not fields in each
lane.

The `--self-test` mode flips one bit independently in a product, rotation,
matrix entry, permutation entry, quotient estimate, and Montgomery factor. All
six corruptions must be rejected.

## Focused evidence

- MSVC 19.50 compiled the constants through the existing
  `blst_avx512_compat_tests` target with one build job.
- `python scripts/audit/validate_matrixnok_constants.py --self-test` passed.
- CTest registers the same validation as `audit.matrixnok_constants` when a
  Python 3 interpreter is available; this adds no test executable.
- `generate_matrixnok_lean.py --check` validates that the committed Lean
  snapshot contains every current C++ literal and the exact current source
  hash. CTest registers this as `audit.matrixnok_lean_constants`; it also adds
  no executable.
- The initial Lean P1 proof establishes dimensions, modulus ranges, oddness,
  pairwise coprimality, exact basis products, dynamic range, rotation and
  target inverses, Montgomery factors, multiplier relations, fixed-table
  shapes/ranges, and the cyclic R1 permutation relation.
- The second Lean P1 packet independently specifies CRT factors from checked
  inverse witnesses and regenerates every r1, r2, convert-to, and convert-from
  matrix, correction vector, and quotient-estimate vector. It proves the
  ceiling estimate's exact one-denominator error enclosure and proves that
  eight 52-bit inputs multiplied by 64-bit estimates cannot wrap the C++
  `uint128_t` accumulator. Encoded tower-constant derivations and the complete
  per-kernel intermediate-bound ledger remain open.

Native execution is not required for constant regeneration. The candidate
arithmetic that consumes these tables will require portable differential tests
and native Peach evidence before selection.

`matrixnok_fixed.h` derives the remaining reducer multiples, Montgomery
factors, and shifted corrections at compile time from these checked constants.
Focused tests compare its stateless operations directly with the former generic
runtime-composed matrices and reducers. The derived-header SHA-256 is
`5938d841014795e97715107ab8a1345e036be2323755890c7377e03b63309a2d`.
