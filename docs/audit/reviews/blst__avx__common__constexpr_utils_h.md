# `blst/avx/common/constexpr_utils.h` Review

## Scope

This record covers the wide-integer compatibility layer at commit `33f6a8b`.
The remaining constexpr parsing, division, mask, inverse, and literal helpers
remain pending for the broader foundational RNS utility review.

## Contracts

- On Clang and GCC, `int128_t` and `uint128_t` are the compiler's native
  signed and unsigned 128-bit integer types.
- On native MSVC, both carriers are exactly 16 bytes. The signed carrier uses
  a two-limb two's-complement representation and the unsigned carrier uses
  low/high limbs.
- Addition, subtraction, comparison, shifts, and multiplication retain the
  low 128 bits. Every current signed producer is separately bounded to the
  representable signed range.
- Constant evaluation uses portable limb arithmetic. Native MSVC runtime
  64-by-64 multiplication uses `_umul128` and returns the exact low/high limbs.
- MatrixNoK quotient accumulation multiplies two 64-bit inputs. Its optimized
  production path must contain direct machine multiplies, not calls to a
  generic wide-integer operator.

## Findings

- `AUD-239`: MSVC previously defined `int128_t` as `int64_t`, invalidating the
  stated overflow-safe bound arithmetic and truncating Debug wide-value
  reconstruction.
- `AUD-240`: the portable 32-bit decomposition used for MSVC unsigned products
  was emitted as eight out-of-line `operator*` calls in one production Fp
  multiplication.

## Validation

- Compile-time assertions cross both sides of `INT64_MIN`/`INT64_MAX`, multiply
  signed boundary values, and exercise the portable constexpr unsigned path.
- Runtime tests compare all 64 pairs from eight signed boundary values with
  MSVC `_mul128`, and cover unsigned high/low multiplication and 128-by-64
  modulo behavior.
- The serial portable-Debug target passes all 32 focused MatrixNoK/C API tests.
- Optimized MSVC emits one `mul` and no branch or call for the standalone
  64-by-64 helper. Production MatrixNoK Fp multiplication emits eight inline
  `mul` instructions and no wide-integer call or conditional jump; before the
  repair it emitted eight calls to the custom `uint128_t operator*`.
- The structural dispatch validator rejects the former 64-bit signed alias and
  requires the 16-byte signed carrier plus constant/runtime split.

## Residual work

The rest of `constexpr_utils.h` is not complete. In particular, its independent
`ceil_div`, mask construction, modular inverse, modulo-zero behavior, and
384-bit literal parser need their own contract and caller review. Native
compiler lowering remains an external trust boundary; the inspected MSVC
assembly is evidence, not a formal compiler proof.
