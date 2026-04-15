# Production Readiness Audit Plan

## Purpose

This document defines a structured audit plan to move the current `blst`
implementation from active development state to production-ready state.

The immediate trigger is the class of bugs we just found:

1. implementation bugs caused by invalid memory contracts, and
2. public-surface bugs where the implementation works but the exposed API is
   incorrect or unusable.

The goal is not to do ad hoc cleanup. The goal is to systematically reduce the
probability of correctness, safety, and integration failures in the code paths
we expect downstream users to trust.

## Why We Are Doing This

The implementation has gone through substantial development. That means there
is now real value in the code, but it also means the code carries development
history:

- rewrites from upstream C into C++
- split headers and compatibility layers
- scalar, AVX, CUDA, and binding surfaces
- internal layout assumptions
- performance-oriented code that makes aliasing and initialization rules easy to
  get wrong

The keygen bug is a representative example. The algorithm was correct, but the
scratch-buffer contract was wrong. That is exactly the kind of bug that escapes
casual testing and only shows up when an implementation is treated as a system
that must obey precise invariants.

So this plan exists to replace opportunistic review with a deliberate audit.

## Goals

- identify and fix correctness bugs caused by invalid memory, lifetime,
  aliasing, or layout assumptions
- identify and fix public API regressions and header-surface regressions
- reduce divergence risks introduced by upstream-to-fork rewrites
- add regression coverage for each bug class we find
- define a repeatable audit process so future changes do not reintroduce the
  same classes of problems

## Non-Goals

- no broad style cleanup for its own sake
- no benchmark-driven optimization work unless it directly supports audit
  coverage
- no large architectural rewrite unless an audit finding proves it is required
- no claim of formal verification

## Definitions

- `memory contract`: the required relationship between a buffer's full logical
  contents and the code that reads or writes it. Example: a reducer that reads
  all 512 bits of a widened object requires all 512 bits to be initialized.
- `API contract`: the required behavior and declarations promised to external
  users. Example: including `blst.h` must expose the documented public surface.
- `differential audit`: a review method that compares forked code against the
  upstream implementation or an earlier trusted version to find semantic drift.
- `invariant`: a property that must always hold if the implementation is
  correct. Example: the same key-generation inputs must produce the same output.
- `gate`: a condition that must be satisfied before the next phase may begin.
- `trusted reference`: a more trusted implementation or behavior baseline used
  to check correctness. In this repository that usually means upstream `blst`,
  upstream VROOM where relevant, or a scalar path that is already well covered.

## Audit Principles

- every serious bug class must become an explicit review checklist item
- every bug fix must add a regression test when practical
- high-risk code is reviewed by contract first and by style second
- upstream-derived code must be audited against upstream semantics, not just
  local intention
- public headers and bindings are part of the implementation, not an afterthought
- no phase is complete without both code review output and validation output

## Current Baseline

The audit is already in progress. The plan should be read as the next-pass
guide, not as if the repository were still at an untouched baseline.

Work already completed:

- fixed zero-length scalar decode handling in the C and C++ scalar helpers
- fixed `blst_uniq_sizeof(0)` underflow behavior
- fixed HKDF expand allocation and failure handling in key generation paths
- fixed `miller_loop_n` empty-input and malformed leading-null handling in
  scalar, AVX, and C API wrappers
- fixed public `wbits` range validation and malformed pointer-array rejection
  in the multipoint C API
- fixed public header and binding drift, including `BLST_BAD_SCALAR` and the
  missing scalar helper declarations
- fixed an AVX unaligned typed-store site
- added regression coverage for the findings above across binding-header,
  C API, and backend-alias tests
- enabled Linux sanitizer builds through `BLST_ASAN` and `BLST_UBSAN` without
  requiring dedicated presets
- used the new sanitizer path to catch and fix a signed-overflow bug in GT
  scalar splitting, with the two's-complement assumption now documented at the
  conversion sites

Implications for the next audit pass:

- do not spend time rediscovering the already-fixed bug classes above unless a
  new surface appears to reintroduce them
- prefer looking for adjacent classes of bugs in code with similar contracts
- treat the sanitizer build as available infrastructure and use it to drive the
  next pass rather than treating sanitizer support as future work

## Risk Categories

The audit will treat the following categories as first-class:

### 1. Memory-Contract Risk

Examples:

- partial writes into widened buffers
- reducers consuming bytes or limbs that were never initialized
- stale stack data influencing arithmetic
- zeroing the wrong object after splitting a former in-place buffer into input
  and output objects

### 2. Layout and Aliasing Risk

Examples:

- invalid `reinterpret_cast`
- ABI assumptions without `sizeof` and `alignof` proof
- semantic type punning where two types share layout but not invariants
- overlap-sensitive code with no aliasing tests

### 3. Rewrite Drift Risk

Examples:

- upstream C code rewritten into C++ with different scratch behavior
- split helper functions changing ordering or lifetime assumptions
- manual “equivalent” rewrites that silently changed boundary handling

### 4. Public-Surface Risk

Examples:

- umbrella headers that no longer expose documented APIs
- bindings out of sync with implementation
- C API declarations or layout expectations that no longer match implementation

### 5. Validation Coverage Risk

Examples:

- deterministic functions without determinism tests
- decode and conversion functions without boundary tests
- code that only has benchmark coverage but no correctness coverage

## Final Intended State

At the end of this plan, the intended state is:

- every high-risk memory-contract site has been reviewed explicitly
- every ABI-sensitive cast site has proof or has been removed
- every public header documented as supported compiles and exposes the intended
  surface
- every major deterministic cryptographic entry point has deterministic
  regression tests
- scalar and AVX code paths have parity coverage where they are expected to
  agree
- sanitizer and fuzzing runs exist for the code that can realistically be
  instrumented

## Execution Rules

- do not audit files in random order; follow the phase and priority lists below
- each finding must be tagged to one of the risk categories above
- each fixed finding must either add a regression test or record why a test is
  impractical
- if upstream has a directly corresponding implementation, differential review
  is mandatory
- if a fix changes public behavior, header and binding coverage must be updated
  in the same change

## Validation Rules

For any completed audit phase, validation must answer:

1. What was reviewed?
2. What findings were made?
3. What fixes were applied?
4. What tests were added or extended?
5. What validation commands passed?
6. What residual risks remain?

## Phase 0: Freeze Scope And Baseline

### Deliverables

- record the code areas in scope for the first audit pass
- record the validation environments we will rely on
- record the risk categories and review checklist

### Required Environments

- local Windows debug build for development and API-surface checks
- Linux sanitizer build for instrumented scalar-path checks
- AWS AVX host for native AVX validation and parity checks

### Initial In-Scope Areas

- scalar arithmetic and serialization
- scalar keygen and EIP-2333
- public C API and binding headers
- AVX base-field and RNS conversion glue

### Exit Criteria

- this plan is present
- the initial file-priority list is recorded
- the review checklist is fixed enough to use in subsequent phases

## Phase 1: Memory-Contract Audit

### Scope

Review all code sites where full logical objects may be only partially written
before full consumption.

### Primary Search Targets

- `limbs_from_be_bytes`
- `limbs_from_le_bytes`
- `vec_zero`
- `bytes_zero`
- `memset`
- `UFr`
- `UFp`
- `UFp2`
- local `scratch` structs

### Review Checklist

For each site:

1. What is the full logical size of the destination object?
2. How many bytes or limbs are actually written by the producer?
3. Is the remainder explicitly initialized first?
4. Does any consumer later read the full object?
5. Is this a rewrite of an upstream in-place buffer?

### Highest-Priority Files

- `blst/sign/keygen.cpp`
- `blst/common/bytes.h`
- `blst/scalar/fields/Fr.h`
- `blst/scalar/fields/Fp.h`
- `blst/avx/fields/conversion_inversion.h`
- `blst/avx/cpu/vector/changebase.hpp`

### Exit Criteria

- all grep hits in this category are triaged
- each confirmed bug is fixed or recorded
- each fix has a regression test where practical

## Phase 2: Layout And Aliasing Audit

### Scope

Review all ABI-sensitive and alias-sensitive code.

### Primary Search Targets

- `reinterpret_cast`
- `std::bit_cast`
- `memcpy`
- `union`
- `static_assert(sizeof`
- `static_assert(alignof`
- `standard_layout`
- `trivially_copyable`

### Review Checklist

For each cast or layout bridge:

1. Is the destination layout proven compatible?
2. Is alignment proven compatible?
3. Is the source type standard-layout and trivially copyable?
4. Do the two types share not only bytes but also semantic invariants?
5. Is aliasing permitted by the function contract?
6. Is there an explicit test for the bridged behavior?

### Highest-Priority Files

- `blst/c_api.cpp`
- `blst/avx_c_api.cpp`
- `bindings/blst.hpp`
- `blst/avx/fields/fp_pack.h`
- `blst/avx/cpu/vector/changebase.hpp`

### Exit Criteria

- all public ABI bridges have proof via code or added `static_assert`
- any semantically unsafe type punning is replaced or tightly documented
- alias-sensitive sites have direct tests if the API permits overlap

## Phase 3: Differential Audit Against Upstream

### Scope

Review forked code paths against the upstream implementation or prior trusted
version to detect semantic drift.

### Method

For each chosen function:

1. identify the upstream or trusted baseline version
2. compare the algorithm, scratch use, ordering, and boundary conditions
3. flag all semantic changes, not just large textual changes
4. confirm which changes were intended and which were accidental

### First-Pass Targets

- key generation and EIP-2333 derivation
- scalar serialization and deserialization
- scalar field reduction and Montgomery conversion
- public C API wrappers
- AVX conversion and reduction glue

### Exit Criteria

- each target function has a reviewed upstream comparison
- semantic deltas are documented as intended or fixed as bugs

## Phase 4: Public-Surface Audit

### Scope

Audit the declared public surface, not just the implementation bodies.

### Areas

- `bindings/blst.h`
- `bindings/blst_scalar.h`
- `bindings/blst_avx.h`
- `bindings/blst_middle_product.h`
- `bindings/blst.hpp`
- language bindings that rely on these declarations

### Review Checklist

1. Does the documented umbrella header compile?
2. Does it expose the documented symbols?
3. Do declarations match implementation and layout assumptions?
4. Are binding generators or handwritten bindings still in sync?
5. Are public examples or README guidance still correct?

### Exit Criteria

- documented public headers compile in standalone smoke tests
- representative C and C++ consumer tests pass
- obvious declaration and documentation drift is fixed

## Phase 5: Determinism And Invariant Test Expansion

### Scope

Add tests for pure and deterministic cryptographic operations and critical
conversion invariants.

### Required Test Families

- deterministic key generation across repeated calls
- deterministic master and child key derivation
- serialize/deserialize round trips
- scalar and field conversion round trips
- boundary-value decoding and rejection tests
- scalar versus AVX parity for shared operations
- aliasing tests where overlap is allowed

### Priority APIs

- `blst_keygen`
- `blst_keygen_v3`
- `blst_keygen_v4_5`
- `blst_keygen_v5`
- `blst_derive_master_eip2333`
- `blst_derive_child_eip2333`
- scalar and field encoding helpers
- public C API conversion helpers

### Exit Criteria

- deterministic functions have explicit determinism tests
- conversion-heavy APIs have round-trip and boundary tests
- each bug fixed in previous phases has regression coverage

## Phase 6: Sanitizer And Instrumented Validation

### Scope

Run instrumented builds on the code that can be meaningfully sanitized.

### Tools

- AddressSanitizer
- UndefinedBehaviorSanitizer
- MemorySanitizer where practical

### Strategy

- use Linux sanitizer builds as a standard audit tool, not a one-off exercise
- `BLST_ASAN=ON` and `BLST_UBSAN=ON` are now expected to work for Linux test
  builds without requiring dedicated presets
- prefer scalar or otherwise instrumentable configurations when the goal is to
  find C++ and API-contract bugs, even if faster asm or AVX paths exist
- disable asm or other opaque shortcuts when necessary so instrumented C++ code
  executes
- focus sanitizer runs on scalar, C API, and conversion-heavy code first
- if a sanitizer finding appears in a low-level arithmetic routine, audit the
  mirrored scalar and AVX copies together before declaring the issue fixed

### Validation Targets

- unit tests
- standalone API smoke tests
- deterministic regression tests added in previous phases

### Exit Criteria

- sanitizer baseline commands are documented
- sanitizer runs are green for the chosen instrumentable configurations
- remaining unsanitized surfaces are explicitly called out

### Current Status

- Linux ASan+UBSan test builds have been exercised successfully on the AWS host
- sanitizer support itself already found one real bug in GT scalar splitting
- sanitizer coverage for AVX expert C API tests still depends on whether the
  chosen build exposes that surface
- MemorySanitizer remains aspirational rather than part of the current baseline

## Phase 7: Fuzzing And Adversarial Input Coverage

### Scope

Add or wire fuzzing for parser, decoder, serializer, and public API surfaces
that consume untrusted bytes.

### Priority Targets

- scalar byte decoders
- field decoders
- point deserialize and uncompress functions
- public C API decode wrappers
- pairing and aggregate setup functions that consume serialized input

### Exit Criteria

- at least the highest-risk decode surfaces have fuzz targets
- crashes and sanitizer findings are triaged and fixed

## Phase 8: Native AVX Validation And Parity

### Scope

Validate that native AVX behavior matches expectations on real AVX hardware.

### Required Checks

- native AVX build succeeds on the AWS AVX host
- direct RNS and reduction smoke tests pass
- public API tests pass under native AVX defaults
- scalar and AVX parity checks exist for shared functionality

### Exit Criteria

- AVX validation commands are documented
- parity or equivalence checks pass for the selected surfaces
- known pre-existing unrelated failures are isolated from audit findings

## Plausible Next Directions

The next audit passes should prioritize work that is adjacent to the findings
already made and can make direct use of the validation infrastructure now in
place.

### 1. Decoder And Parser Fuzzing

Start with byte-driven public surfaces:

- scalar byte decoders
- field decoders
- point deserialize and uncompress entry points
- public C API wrappers that accept serialized or externally supplied buffers

### 2. Static Analysis For C++ Contract Bugs

Wire a practical static analysis pass such as `clang-tidy` over the C++ targets
most likely to contain lifetime, narrowing, signedness, and contract mistakes:

- `blst/c_api.cpp`
- `blst/sign/keygen.cpp`
- `blst/scalar/fields/Fp12.cpp`
- binding wrappers and split public headers

### 3. Signedness And Conversion Audit

The sanitizer-discovered GT split bug suggests a broader search for:

- signed/unsigned boundary conversions
- balanced-digit decompositions
- integer narrowing in scalar and AVX helper code
- assumptions about two's-complement or wraparound that are currently implicit

### 4. Binding Validation Expansion

Rust already has meaningful CI coverage. The next likely gap is:

- Go smoke or round-trip tests for the updated keygen and empty-slice behavior
- C# smoke coverage for keygen preconditions and scalar decode failure paths

### 5. Portable Dynamic Analysis

As a secondary dynamic check after sanitizers:

- run a portable scalar Linux build under Valgrind
- keep asm, CUDA, and default-AVX disabled for that pass so the results stay
  focused on instrumentable code

## Work Queue Priority

The audit should begin in this order:

1. `blst/sign/keygen.cpp`
2. `blst/common/bytes.h`
3. `blst/scalar/fields/Fr.h`
4. `blst/c_api.cpp`
5. `blst/avx_c_api.cpp`
6. `bindings/blst.h`
7. `bindings/blst.hpp`
8. `blst/avx/fields/conversion_inversion.h`
9. `blst/avx/fields/fp_pack.h`
10. `blst/avx/cpu/vector/changebase.hpp`

This order is chosen because it attacks the highest-impact correctness and
public-surface risks first.

## Required Deliverables Per Finding

Every confirmed finding should produce:

- file and line reference
- risk category
- concise root-cause description
- fix description
- regression test or stated reason no regression test was added

## Status Tracking

If this plan turns into active execution, maintain a companion status file:

- `docs/production-readiness-audit-plan_status.md`

That status file should record:

- reviewed files
- findings made
- fixes landed
- tests added
- remaining open risks
