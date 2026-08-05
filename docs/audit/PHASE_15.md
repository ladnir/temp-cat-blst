# Phase 15: Scalar Curves and Protocols

## Scope

Phase 15 performs the file-by-file hand audit of scalar curve arithmetic,
serialization, hash-to-curve, pairing, signing, MSM, and related protocol
logic. C, selected native C++, and safe Rust contracts are updated as each
implementation family is reviewed. BTx/BTE is frozen, CUDA is experimental,
and neither is part of this phase's release gate.

Performance-sensitive arithmetic retains narrow preconditions where they
enable faster formulas. Those preconditions must be explicit and consistent
across supported language surfaces; Release validation is added only when the
API is explicitly fallible or checked.

## Scalar point review

The first complete file packet reviewed `blst/scalar/curve/point.h` and
`point.cpp` at commit `07d79a48cc7d82560f5194f76a137e4ee2527d06`.
It traced native and C callers, generated/safe Rust exposure, scalar GT Miller
fallbacks, fixed-window MSM consumers, tests, and the pre-refactor union/macro
lineage.

### Findings and decisions

- `AUD-104` fixed a 96-byte G2 selection from a 48-byte `Fp::one()` object by
  restoring the field-generic `Field::one()` operand. The older C macro had
  selected the larger `radix384.p2` union view; the typed translation lost it.
- `AUD-105` replaced projective-prefix-to-affine references in generators,
  serializers, and the scalar Miller fallback with genuine affine constants
  and typed coordinate helpers. Projective serialization retains one
  normalization buffer and adds no affine copy.
- `AUD-106` removed the unused, undefined internal `mult_w5` declaration.
- `AUD-107` made native high-level decoding transactional on the final G1
  x-zero rejection. The C boundary was already transactional.
- `AUD-108` corrected point timing, aliasing, decoding, XYZZ, generator, and
  `add_n_dbl` contracts and removed reserved/stale header constructs.
- `AUD-096` is accepted under the inherited API: unsuffixed `add` is the fast
  incomplete formula for finite distinct, nonopposite points, while
  `add_or_double`/`dadd` is the complete operation for equality or opposites.
- `AUD-109` remains open in the pending `bulk_addition.h` review. Its scratch
  accumulator casts live projective scratch objects to affine references and
  needs a zero-copy, code-generation-reviewed redesign.

### Validation

- MSVC 19.50 Release scalar focused target: 6/6 passed.
- MSVC 19.50 Debug scalar focused target: 6/6 passed.
- Linux GCC 13 ASan/LSan Debug focused target: 6/6 passed.
- The focused target covers typed generators, projective/affine serialization
  parity, G2 mixed identity selection including exact output aliasing, complete
  equal/opposite addition, failed-output preservation, and valid publication.
- Optimized MSVC inspection shows projective G1/G2 serialization and
  compression pass x/y addresses directly to coordinate helpers. Frames remain
  the existing one-Point normalization size; no affine copy, conversion, or
  second point buffer appears.
- Public documentation and group-contract integrity scripts pass.

No benchmark was run because the arithmetic fix changes only the selected
constant address and the serialization rewrite preserves the prior call/frame
shape. No full suite, CUDA, BTx/BTE, Apple CI, or portable-AVX Release build was
run.

## Status and next action

`point.h` and `point.cpp` are complete at the recorded commit. `GT.cpp` and the
public C/Rust containing files remain pending except for their linked focused
families.

Next, review the performance-sensitive scalar accumulation/MSM files beginning
with `bulk_addition.h`, resolving `AUD-109` before proceeding through
`multi_scalar.h`/`.cpp`. Keep builds single-job and use the focused point/MSM
targets; obtain approval before benchmarks or broad validation.
