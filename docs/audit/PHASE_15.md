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
- `AUD-109` was opened for the then-pending `bulk_addition.h` review. Its
  zero-copy, code-generation-reviewed resolution is recorded in the next
  packet below.

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

## Scalar bulk-addition review

The second complete file packet reviewed
`blst/scalar/curve/bulk_addition.h` and its selected wrapper at commit
`bd45afa58630e7d235d4bbe048e8fc88bc5d5a09`. Supported fixed-table MSM and
C/Rust callers were traced far enough to repair their scratch and routing
boundaries without declaring the larger MSM files reviewed.

### Findings and decisions

- `AUD-109` replaces projective-prefix casts with `BatchAffineNode<Field>`, a
  real affine subobject followed by field workspace. Static assertions pin its
  size, alignment, offsets, trivial lifetime, and limb-storage compatibility.
- Raw C fixed-table scratch is treated as storage: each node lifetime begins at
  its byte offset before gather, and the index-based accumulator never assumes
  that separately constructed objects form a C++ array.
- `AUD-110` restores public C and safe Rust batch-add routing. The scalar batch
  path starts at 16 public points; smaller inputs remain complete serial adds.
- `AUD-111` replaces unchecked allocation with one nonthrowing allocation
  capped at 144 KiB. Failure uses complete serial addition, so the legacy void
  C boundary cannot throw, publish a partial result, or silently return the
  identity.
- `AUD-112` confines internal tree helpers to `scalar::detail`, prefixes the
  two scratch-limit macros, and replaces the generic scratch-size macro with a
  typed constexpr helper.
- `AUD-113` records—but does not change—the separate fixed-table precompute
  `alloca` in `multi_scalar.cpp`. Its 144 KiB normal bound (45 KiB in the
  retained WebAssembly branch) requires an explicit stack-budget decision in
  that file's review.

### Validation

- MSVC 19.50 scalar Debug and Release focused executions pass 11/11.
- G1/G2 coverage includes 0/1/15/16/17 points, exact batch and 144 KiB chunk
  boundaries, equality, opposites, infinity, explicit and compact pointers,
  invalid heads, and injected allocation failure against complete serial
  oracles.
- A raw `blst_limb_t` fixed-table C regression proves that typed lifetimes are
  established in caller-owned scratch and matches individual multiplication.
- Optimized MSVC inspection shows no construction stores, affine conversion,
  or extra x/y copy. The accumulator uses direct 144-byte G1 and 288-byte G2
  addressing; the public wrapper performs one checked allocation and branches
  directly to the serial helper on failure.
- Structural checks find no projective-to-affine scratch cast, unprefixed
  scratch macro, or public re-export of the tree helpers.

No benchmark was run: the layout/copy/code-generation structure is preserved,
and benchmark execution remains permission-gated. No broad suite, sanitizer,
CUDA, BTx/BTE, Apple CI, or portable-AVX Release build ran.

## Status and next action

`point.h`, `point.cpp`, and scalar/selected bulk-addition headers are complete
at their recorded commits. `GT.cpp` and the public C/Rust containing files
remain pending except for their linked focused families.

Next, review `multi_scalar.h`/`.cpp`. Start with `AUD-113`, the native
fixed-table dispatch's allocation/error paths, sizing overflow, raw scratch
contracts, and every exposed helper before changing the established kernels.
Keep builds single-job and use focused point/MSM targets; obtain approval before
benchmarks or broad validation.
