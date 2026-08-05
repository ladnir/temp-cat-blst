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

## Scalar GT and Miller-loop review

The sixth packet completes `blst/scalar/curve/GT.h`,
`blst/scalar/curve/GT.cpp`, and the selected target-group facade. It also
applies the required semantic migrations to retained host consumers without
auditing frozen BTx/BTE or benchmark implementations.

### Findings and decisions

- `AUD-099` is accepted with a narrower residual risk: native C++ now has a
  distinct `MillerResult`, raw-field GT construction is explicit, and C/Rust
  representation-only `blst_fp12` storage is documented. Repeated Release
  membership checks remain intentionally absent; decoder policy stays in
  `AUD-085`.
- `AUD-131` makes the AVX generator the compile-time native representation of
  the canonical generator pairing and removes runtime initialization.
- `AUD-132` preserves constant-time conditional negation while selecting only
  the Fp12 half changed by conjugation.
- `AUD-133` replaces the misleading native `final_mult` name with
  `final_exp` and requires a `MillerResult` input. No compatibility overload
  weakens the domain boundary.
- `AUD-134` replaces faux affine-point storage for `(-2x, 2y)` with
  `PreparedG1LineScale` throughout scalar prepared-line evaluation and middle
  product workspaces.
- Scalar/AVX GT boundary conversion uses direct field conversion rather than a
  576-byte serialization round trip. Empty Miller products return identity;
  malformed nonempty C pointer heads fail closed to the zero sentinel.

### Focused validation

- MSVC scalar Debug production and linked test targets compile serially.
- The focused scalar GT, pairing, C boundary, prepared-line, and middle-product
  packet passes 24/24.
- Portable-AVX Debug production plus five directly affected test translation
  units compile serially.
- Native AVX-512/IFMA production compiles serially.
- Prior generator evidence at `f710040` proves scalar/AVX equality and no
  dynamic-initializer section.

No benchmark was run. Benchmark sources received only mechanical native API
migration and were not audited. Frozen BTx/BTE received only the mechanical
type/name migration required to keep compiling. No full suite, sanitizer,
CUDA, Apple CI, or portable-AVX Release build ran. Experimental CUDA still uses
the old native C++ spelling/domain and remains outside this production gate.

## Scalar map-to-group review

The seventh packet completes `blst/scalar/curve/map_to_G.h`, its selected
facade, and the temporary AVX compatibility wrapper at implementation commit
`a1ef4b38e406a78a95acc5550745c73305760ce5`. The hash-to-field writer, C
boundary, and public header received focused repairs but remain pending their
own complete file reviews.

### Findings and decisions

- `AUD-135` replaces the G2 `Fp2[]`-to-`Fp*` flattening cast with typed field
  subobject access and a compile-time destination-extent proof.
- `AUD-136` restores immutable host static coefficient storage. Experimental
  CUDA retains automatic device constants and remains outside this gate.
- `AUD-137` removes affine inversions from the AVX projective bridge and routes
  scalar-storage C encode/hash directly to the scalar implementation.
- `AUD-138` makes the retained four RFC 9380 G1/G2 NU/RO outputs exact tests,
  with exceptional-zero and alias regressions beside them.
- `AUD-139` records precise C, selected native C++, and generated Rust
  contracts and removes raw map helpers from the selected namespace.
- `AUD-140` constrains internal templates to G1/G2, documents alias behavior,
  and asserts the nonaliasing cofactor-clear precondition in Debug.

### Focused validation

- MSVC scalar Debug and Release map targets pass 6/6 each.
- Portable-AVX Debug passes 7/7 map tests and 5/5 focused AVX C API contracts.
- Native AVX-512/IFMA Debug compiles serially; the unsupported local MSVC
  native binary was not executed.
- Optimized G1/G2 isogeny stack frames fall from 4,240/2,272 to 1,696/1,040
  bytes; SSWU frames fall from 640/1,024 to 448/736 bytes. G1 no longer emits
  `__chkstk`, and per-call coefficient initialization is absent.
- Structural checks find no field-flattening cast, affine AVX bridge, selected
  raw-helper export, or redundant C encode/hash AVX round trip.

No benchmark, full suite, sanitizer sweep, CUDA build, BTx/BTE work, Apple CI,
or portable-AVX Release build ran.

## SHA-256 review

The eighth packet completes `blst/sign/sha256.h` at implementation commit
`5c69ad63c76d9e63e9cf75a708ac84ede45b5eb6`. The fixed SHA-256 compression and
padding schedule matches the pinned BLST lineage, and all
message-content-dependent work remains branchless. Five findings are validated:

- `AUD-141` restores the supported C/C++ declaration for the retained
  `blst_sha256` export.
- `AUD-142` adds independent end-to-end one-shot and streaming known-answer
  evidence, including padding boundaries.
- `AUD-143` separates fast functional buffer zeroing from secure lifecycle
  erasure.
- `AUD-144` securely clears the public one-shot wrapper's locally owned
  context after digest publication.
- `AUD-145` makes context state, invariants, lifecycle, and integer-width
  assumptions explicit while preserving the `SHA256_CTX` alias required by
  frozen BTx/BTE.

### Focused validation

- Scalar Debug passes 5/5 SHA audit tests, including context-invariant death
  tests; scalar Release passes the four applicable functional tests.
- Standard empty, `abc`, multi-block, and million-`a` vectors pass. Independent
  deterministic vectors cover both sides of the 56-byte padding boundary and
  the 64-byte block boundary over one-shot and several streaming schedules.
- The preexisting keygen, EIP-2333, selected hash-to-group, and exact RFC XMD
  packet passes 13/13 after the context migration.
- The public header compiles as strict MSVC C11 with `/W4 /WX`; the public API
  documentation and repository manifest checkers pass.
- MSVC 19.50 Release emits eight 64-bit stores for functional buffer zeroing,
  and the one-shot wrapper emits a fixed 112-byte volatile erase after digest
  publication. Its object has exactly one bare C SHA symbol, `blst_sha256`.

No benchmark, full suite, sanitizer sweep, CUDA build, Apple CI, or
portable-AVX Release build ran. Frozen BTx/BTE compiled through the retained
internal type alias but received no source edit or audit.

## Key-generation review

The ninth packet completes `blst/sign/keygen.cpp` at implementation commit
`f36e35b3831ec7b527e767e8221c00b53ff739a7`. The six retained upstream C
symbols, signatures, version selection, valid-input outputs, null-info/default-
salt behavior, and child exact-alias behavior remain unchanged.

### Findings and decisions

- `AUD-146` replaces inactive-union-member HMAC access with active limb storage
  and a permitted unsigned-byte object-representation view.
- `AUD-147` adds direct RFC HMAC/HKDF, output-limit, version-compatibility, and
  EIP-2333 evidence.
- `AUD-148` restores internal linkage to implementation-only helpers while
  retaining all six required cross-translation-unit and bare C symbols.
- `AUD-149` documents precise C/C++/safe-Rust contracts, generated Rust
  inheritance, Debug invariants, and nonthrowing C++ boundaries.
- `AUD-150` separates functional initialization from final secure erasure and
  initializes only the two scalar limbs not overwritten by the 48-byte import.
- `AUD-151` removes the reserved double-underscore HKDF test-mode spelling.

### Focused validation

- Scalar MSVC Debug passes 7/7 keygen tests; Release passes 6/6.
- RFC 4231 HMAC, RFC 5869 HKDF, the 255-block limit, all versioned keygen
  outputs, upstream null behavior, child exact aliasing, and four published
  EIP-2333 master/child vectors pass.
- Both safe-Rust signature variants pass the focused empty-v5-salt regression.
- MSVC 19.50 Release uses eight ordinary 64-bit HMAC zero stores, two high-limb
  scalar initialization stores, and fixed volatile final erasures over the
  complete 416-byte and 8,576-byte secret scratch objects.
- Internal helper symbols are Static, the C boundary emits exactly the six
  upstream names, and the keygen object has no heap-allocation reference.

No benchmark, full suite, sanitizer sweep, CUDA build, Apple CI, or portable-
AVX Release build ran. Frozen BTx/BTE was neither edited nor audited.

## Status and next action

Scalar point, bulk addition, MSM, fixed-window multiplication, GT/Miller, and
map-to-group, SHA-256, and key generation are complete at their recorded
commits. Continue Phase 15 with `blst/sign/aggregate.cpp` and `aggregate.h`.
Preserve the now-documented keygen and group-domain contracts while tracing
signature aggregation and verification. `AUD-085` remains the explicit GT
deserialization decision.
