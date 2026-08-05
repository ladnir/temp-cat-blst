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

## Scalar multi-scalar review

The third complete packet reviews `blst/scalar/curve/multi_scalar.h` and
`multi_scalar.cpp`, the selected native facade, and the safe Rust multipoint
wrapper at commit `af09f66ee9097fe35617af33fbcbc7dfd4229237`.

### Findings and decisions

- `AUD-113` replaces the 144 KiB input-dependent `alloca` with fast output-
  table reuse, one optional nonthrowing allocation capped at 144 KiB, and a
  fixed 8 KiB stack fallback. The fallback is a cold non-inlined helper, so
  primary precompute entry points no longer carry its stack probe.
- `AUD-114` stops compact point/scalar descriptor traversal at exactly
  `npoints`; exact n-entry arrays no longer require an undocumented sentinel.
- `AUD-115` defines the extractor-derived Pippenger window range as `[1,24]`
  and makes explicit scratch shifts overflow-safe.
- `AUD-116` defines zero-point/zero-bit identity behavior before input access
  and consistently rejects more than 256 scalar bits.
- `AUD-117` makes tiled Pippenger honor the exact public `bit0` rather than
  rounding it down to a window multiple.
- `AUD-118` removes the remaining Point/Affine/Field prefix-object violations.
  General conversion uses real affine x subobjects; precompute uses aligned
  byte storage with explicit trivial Point/Field lifetime start.
- `AUD-119` removes the unused native dispatcher and its unchecked allocation.
- `AUD-120` removes production scalar/pointer/bucket printing and limits the
  retained hash helper to experimental CUDA builds.
- C, selected native C++, and safe Rust contracts now agree on subgroup/public
  data, scalar layout/width, windows, scratch, compact arrays, empty results,
  overlap, allocation, and failure behavior.

### Validation

- MSVC 19.50 scalar Debug focused packets pass 15/15, including maximum-window
  G1/G2 tables, packed/pointer/aliased layouts, exact descriptor arrays,
  nonaligned tile offsets, a heap-sized case, and cap-forced fallback cases.
- MSVC 19.50 scalar Release focused packets pass 15/15.
- Full-Rust-binding layout and G1/G2 multipoint tests pass 3/3.
- Optimized MSVC has no `alloca` or production debug symbols. The only
  `__chkstk` calls are the two 8,344-byte cold fallback instantiations; primary
  Fp/Fp2 precompute functions have no stack probe.

No benchmark was run because arithmetic kernels and the fast in-table schedule
are unchanged, and object inspection confirms the cold frame is isolated. No
broad suite, sanitizer, CUDA, BTx/BTE, Apple CI, or portable-AVX Release build
ran.

## Status and next action

`point.h`, `point.cpp`, scalar/selected bulk addition, scalar/selected MSM, and
the safe Rust multipoint wrapper are complete at their recorded commits.
`GT.cpp` and the larger public C
containing file remain pending except for linked focused families.

Next, review `multi_scalar_batch_affine.cpp`. Preserve its fixed-size hot
kernels, trace every scratch/layout contract into C and Rust, and batch its
findings before proposing repairs. Keep builds single-job and focused; obtain
approval before benchmarks or broad validation.

## Scalar batch-affine MSM review

The fourth packet reviews `blst/scalar/curve/multi_scalar_batch_affine.cpp`
and its explicit C/native/raw-Rust boundaries. The implementation is a local
Hypercat/CatBLST extension derived from upstream Pippenger extraction and
complete point formulas; the pinned upstream BLST tree has no matching
batch-affine implementation file.

### Findings and decisions

- `AUD-121` narrows the explicit window contract from the incorrect `[1,32]`
  claim to the shared extractor's real `[1,24]` domain. Production selection
  remains 10 through 15.
- `AUD-122` replaces typed array casts over limb scratch with aligned,
  byte-indexed lifetime-aware views. The compile-time batch sizes, array
  layout, and direct address arithmetic remain unchanged.
- `AUD-123` rejects more than 256 scalar bits and more than `UINT32_MAX`
  scheduled points before input conversion or extraction.
- `AUD-124` adds deliberate equal-point, opposite-point, and projective-
  collision differential cases for both G1 and G2.
- `AUD-125` makes the expert C/native contract precise about public data,
  subgroup/nonidentity inputs, windows, scratch, overlap, empty output, and
  failure atomicity; generated raw Rust inherits the C contract.

### Validation

- MSVC scalar Debug and Release focused execution each pass 7/7.
- Strict MSVC `/W4 /WX /O2` compilation passes.
- Optimized MSVC contains no lifetime-start calls or construction stores; the
  placement operations are abstract-machine-only. It also contains no
  scratch-wrapper or allocation calls.
- Public documentation, group-contract, manifest, and structural checks pass.

No benchmark ran because optimized inspection found no changed hot-loop work.
No broad suite, sanitizer, CUDA,
portable-AVX Release, Apple CI, or BTx/BTE audit is included.

## Scalar fixed-window multiplication review

The fifth packet reviews `blst/scalar/curve/ec_mult.h` and its direct scalar
point callers. It compares the C++ templates with pinned upstream
`54e6e55674722fc2797ebb4bbb71b26d881eb4b8:src/ec_mult.h`, including the
low-order correction in upstream commit `7c535f1`.

### Findings and decisions

- `AUD-126` constrains the fixed-window templates to their actual `[2,8]`
  domain. Production remains width four/five, so the hot layout is unchanged.
- `AUD-127` makes the allocation-hardened status non-discardable and records
  the precomputed GLV/GLS success invariant without a Release branch.
- `AUD-128` restores direct order-11 G1 and order-13 G2 evidence for the two
  inherited exact-scalar `unchecked_mult` exceptions.
- `AUD-129` masks unused source bits for future non-byte-aligned wNAF
  instantiations. Current 64/128-bit production forms compile unchanged.
- `AUD-130` removes volatile erasure only from explicitly public variable-time
  decomposition. Secret-capable constant-time and signing paths retain it.
- The internal contracts now state extraction masking, fixed-window bounds,
  table ownership, compact descriptors, aliasing, and timing. The reserved
  double-underscore scalar and forwarding-header guards are replaced.

### Focused validation

- MSVC scalar Release window/allocation tests pass 6/6.
- The retained non-subgroup oracle and the new low-order G1/G2 packet pass 2/2.
- Selected CT/VT agreement, exact-alias, full-width boundary, and independent
  crossover-oracle tests pass 4/4, for a focused Release total of 12/12.
- MSVC `/W4 /WX /O2` compilation of the owning point translation unit passes.
- Optimized 64/128-bit wNAF recoders return without erasure loops. Both
  precomputed GLV/GLS functions inline `mult_window_n` and contain no
  out-of-line helper call or allocation reference.

No benchmark, broad suite, sanitizer, CUDA, BTx/BTE, Apple CI, or portable-AVX
Release build is included.
