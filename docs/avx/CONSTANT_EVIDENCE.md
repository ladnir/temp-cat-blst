# MatrixNoK Fixed-Parameter Evidence

## Scope

This record covers only the isolated BLS12-381 MatrixNoK parameter import. The
tables are not connected to production arithmetic or dispatch in this packet.

- Upstream commit: `d6ecd534c54a41702ea80ae75b52250b6d51b372`
- Upstream path: `scripts/bls12_381_intrns4_params.hpp`
- Upstream Git blob: `d0e0437007f4dec6f193d528ce618f8b516e4fbe`
- CatBLST path: `blst/avx/generated_matrixnok_constants.hpp`
- Imported-header SHA-256:
  `b9b59ac2bfaace9bae2966fc746fae4f73a6276dbff01183c04d768f77ffa78a`
- Namespace: `blst::avx::matrixnok::generated`

The data retains VROOM's MIT license and Simon Langowski attribution. The C++
translation uses immutable inline constants and requires no runtime GMP,
initialization, search, allocation, or mutable global state.

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
- exact equivalence of the cyclic r1 permutation table to the base matrix.

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

Native execution is not required for constant regeneration. The candidate
arithmetic that consumes these tables will require portable differential tests
and native Peach evidence before selection.
