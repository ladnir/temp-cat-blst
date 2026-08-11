# Deleted `blst/avx/rns/inversion.h` Review

## Scope and disposition

This record covers the complete deleted VROOM addition-chain header through
removal commit `81f21ac`. The file had no CatBLST includes or instantiations.
Active AVX base-field and Fp12 inversion is owned by
`blst/avx/fields/conversion_inversion.h`.

## Arithmetic review

For BLS12-381 parameter magnitude `x = 0xd201000000010000`, independent
symbolic tracing establishes:

- `exp_x(a) = a^x`;
- `exp_x_m3(a) = a^(x-3)`;
- the complete schedule uses 425 modular products and returns
  `a^(x^6 + 2*x^5 - 2*x^3 - x - 3)`;
- the exponent equals `3p-4`, which is congruent to `p-2` modulo `p-1`.

The chain therefore computed the inverse for nonzero Fp inputs. It returned
zero for zero, for which the mathematical inverse is undefined. Existing loop
counts were fixed and introduced no input-dependent control flow, assuming the
ring multiplication itself was constant time.

## Finding and remediation

`AUD-258` records that a standalone-test implementation was retained as a
production file despite having no callers. The unused generic interface also
stored one `Ring` reference while accepting and ignoring another, exposed all
chain helpers publicly, accepted a signed runtime shift count, and lacked an
include guard, provenance, and an explicit nonzero precondition.

The file was deleted rather than repaired. The production dispatch validator
rejects its restoration. No active symbol, call graph, code generation, public
API, or runtime behavior changed.

## Evidence and residual risk

- repository-wide symbol/include search found zero CatBLST consumers;
- VROOM provenance shows the header used by its standalone
  `test_inversion.cpp`, not the production conversion/inversion owner;
- independent integer exponent tracing verifies the complete chain;
- `validate_matrixnok_dispatch.py` passes after deletion.

The active scalar/AVX conversion and inversion bridge is a separate review
owner. No build was required for a file absent from every translation unit.
