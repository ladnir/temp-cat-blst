# File Review: `bindings/blst_scalar.h`

## Record

| Field | Value |
| --- | --- |
| Path | `bindings/blst_scalar.h` |
| Roadmap phase / queue | Phase 18 public C/C++/Rust surface; Wave 1 Set B |
| Reviewed commit | `16e5b49` |
| Reviewer | Codex audit team |
| Date | 2026-08-12 |
| Disposition | in-scope production |
| Subsystem | plain scalar C ABI and C++-callable compatibility surface |
| Profiles | C99-or-later, C++20, scalar, selected AVX-512, Rust FFI |
| Generation/provenance | handwritten split of upstream BLST's public header, extended by CatBLST APIs |
| Upstream references | `54e6e55674722fc2797ebb4bbb71b26d881eb4b8:bindings/blst.h` |
| Related files | `bindings/blst_types.h`, `blst/common/error_codes.h`, `blst/exports.cpp`, `blst/c_api.cpp`, selected field/curve/sign implementations, `bindings/rust/src/bindings.rs`, `bindings/rust/src/primitives.rs`, `bindings/rust/src/pairing.rs`, `bindings/rust/src/pippenger.rs`, `bindings/rust/src/lib.rs`, `test/BindingHeader_tests.cpp` |

This record covers the exact implementation commit above. All 253 function
declarations and the two opaque context types were inspected. It does not
substitute for the implementation-file reviews referenced by each family.

## Purpose and trust boundary

This is the authoritative plain C declaration surface for scalar CatBLST. It
is also directly callable from supported C++, and its layout and calling
conventions define the raw Rust FFI. It crosses all important trust boundaries:
untrusted byte decoding, scalar and secret-key import, subgroup validation,
pairing verification, hash-to-curve, caller-owned scratch, and caller-owned
opaque context storage.

The header deliberately exposes low-overhead operations with narrow contracts.
Ordinary field and point arithmetic does not repeatedly validate trusted
in-memory operands. Optimized G1, G2, and GT multiplication and MSM require
previously established subgroup operands. Point decoding establishes canonical
encoding and on-curve status, not subgroup membership. Only the two inherited
`blst_p1_unchecked_mult` and `blst_p2_unchecked_mult` operations promise exact
scalar action for on-curve points outside the prime-order subgroup.

Secret inputs include secret-key scalars and canonical scalars supplied to
constant-time multiplication. Public inputs include verification messages,
encodings, pairing weights, FFT/MSM scalars, counts, windows, and scratch
shapes. The declarations distinguish constant-time secret-capable paths from
explicitly variable-time public paths.

## Contracts and invariants

### Common boundary rules

- Output pointers name complete writable objects; input pointers name complete
  readable objects unless a declaration explicitly permits null.
- Exact aliasing is allowed only where stated. Partial overlap is not
  supported. Multi-output butterflies require distinct `x0` and `x1`.
- Fallible operations state whether output or context is unchanged on failure.
- No C++ exception may cross the C ABI. Legacy `void` boundaries rely on
  allocation-free/nonthrowing implementation or documented preconditions.
- `blst_pairing` and `blst_uniq` have distinct incomplete tags. Their queried
  sizes and alignments govern caller-owned storage; the pointer types are not
  interchangeable.
- C++-only optional null/zero defaults are provided through the locally scoped
  `BLST_CXX_DEFAULT` macro, which expands empty in C and is undefined before
  the header ends.

### Complete family inventory

1. **SHA-256:** one-shot hashing accepts a null message only for zero length,
   permits output/input overlap, publishes 32 bytes after consuming input, and
   clears the local context.
2. **Scalar and secret-key helpers:** `blst_scalar` is arbitrary 256-bit
   little-endian storage, not proof of reduction. Fixed-width imports preserve
   all bits. `blst_scalar_fr_check` establishes `[0,r)`;
   `blst_sk_check` establishes `[1,r)`. Variable-length imports reduce modulo
   `r` and report whether the result is nonzero. Secret-key arithmetic assumes
   canonical nonzero inputs and publishes the modular result even when it is
   zero.
3. **Fr:** reduced Montgomery operands, exact result/input aliasing, public
   exponent for `pow`, constant-time conditional negation and selected inverse,
   explicit Montgomery conversions, and distinct butterfly outputs.
4. **Fp and Fp2:** reduced Montgomery operands, exact aliasing, variable-time
   Euclidean inverse versus fixed-chain inverse, explicit square-root failure
   semantics, and constant-time conditional negation.
5. **General Fp12 and GT:** general arithmetic remains distinct from
   cyclotomic and GT-only operations. Frobenius accepts only powers 1--3 and is
   nodiscard/failure-atomic. Final exponentiation consumes a nonzero general or
   Miller value. GT arithmetic requires GT membership without Release checks;
   `mult` is constant-time for canonical scalars at fixed public width, while
   `mult_vartime` is public-only.
6. **G1 and G2 primitives:** basic arithmetic, representation conversion,
   observation, and serialization accept on-curve operands. Incomplete `add`
   requires distinct, nonopposite finite inputs; `add_or_double` is the complete
   alternative. Membership predicates check curve and subgroup. Conditional
   negation is constant-time.
7. **G1 and G2 multiplication:** ordinary `mult` and `mult_vartime` require
   subgroup points and accept scalars modulo `r`; only canonical inputs receive
   the constant-time guarantee. Legacy `unchecked_mult` accepts on-curve points
   and applies the exact scalar without subgroup endomorphisms.
8. **Batch conversion, fixed-table MSM, Pippenger, tiled Pippenger, and native-Fr
   MSM:** advanced operations require subgroup points and public scalars.
   Compact descriptors use a non-null head followed by optional contiguous-null
   entries. Counts, widths, windows, scratch bytes/alignment, empty behavior,
   disjointness, and failure publication are stated. Zero-point ordinary
   Pippenger scratch queries return zero.
9. **Point serialization and decoding:** serialization requires an on-curve
   point. Decode validates flags, canonical coordinates, infinity, and the
   curve equation but not complete subgroup membership; output is unchanged on
   failure. Pointer-only aggregate decoding states its first-byte-dependent
   readable extent.
10. **Miller and pairing:** Miller results are pre-final Fp12 values, not GT.
    Pairing contexts borrow only DST; messages and augmentations are consumed
    during the call. Orientation, nullable contributions, validation flags,
    commit/merge lifecycle, failure atomicity, alignment, ownership, and lack
    of internal synchronization are explicit.
11. **Aggregate/core verification:** ordinary aggregation validates signatures
    but assumes prior PK subgroup validation; `chk_n` variants transfer checks
    according to their flags. Weighted operations use public variable-time
    scalars of at most 256 bits. Core verification checks both subgroups and
    rejects infinite public keys.
12. **Keygen, EIP-2333, hash/encode/map, and signing:** versioned keygen
    documents IKM, salt, info, retry, and invalid-input behavior. Hash/map
    produces subgroup points. Signing and public-key derivation require a
    canonical nonzero key and subgroup hash without Release validation.
13. **Uniqueness and XMD:** uniqueness storage is caller-owned, externally
    synchronized, and borrows inserted messages. XMD states its RFC 9380 output
    bound, SHA-256 transcript bound, pointer/overlap rules, oversized-output and
    oversized-input behavior, and allocation-free exception contract.
14. **ABI queries:** point, affine-point, Fp12, pairing, and uniqueness queries
    expose exact storage requirements without constructing objects.

### Supported declaration surfaces

The C declarations in this file are authoritative. Supported C++ may call the
same `extern "C"` functions and also uses the selected native C++ APIs in
`blst/fields`, `blst/curve`, and `blst/sign`; the domain and timing vocabulary
was compared against those declarations. Rust raw declarations cover every C
symbol. Safe Rust documents and enforces slice lengths and owned pairing/DST
lifetime where wrappers exist; expert field, scratch, and native-Fr operations
remain raw FFI rather than being presented as automatically safe abstractions.

## Dependency and call-path map

- `blst_types.h` supplies all POD layouts, fixed-width aliases, `BLST_ERROR`,
  `BLST_NODISCARD`, and C/C++ portability macros.
- `blst/exports.cpp` implements scalar, Fr, Fp, Fp2, parsing, conversion, and
  SHA entry points.
- `blst/c_api.cpp` implements Fp12/GT, point, MSM, Miller, pairing adapters,
  keygen, map/hash/sign, uniqueness, XMD, and size queries.
- Selected aliases route production point/GT work to scalar or AVX-512 at
  compile time while retaining scalar boundary layouts. This header exposes no
  runtime dispatch and creates no conversion itself.
- `bindings/rust/src/bindings.rs` mirrors all raw symbols. `primitives.rs`,
  `pairing.rs`, `pippenger.rs`, and `lib.rs` provide the reviewed safe or typed
  wrappers for group, pairing, MSM, key, and protocol operations.
- `blst.h` reaches this file through `blst_aux.h`; direct inclusion remains
  supported.

## Upstream/specification comparison

The exact pinned upstream comparison is
`54e6e55674722fc2797ebb4bbb71b26d881eb4b8:bindings/blst.h`. All 176 pinned
upstream function names remain present. CatBLST adds 77 declarations for native
field utilities, explicit GT operations, variable-time multiplication,
caller-scratch MSM, EIP-2333, uniqueness, XMD, and ABI queries.

The upstream C++ optional defaults were restored at the reviewed commit.
Renaming upstream's generic `byte` and `limb_t` aliases to project-prefixed C
types is ABI-neutral. The intentional `blst_fp12_frobenius_map` return change
from `void` to `bool` is the previously approved `AUD-080` safety repair; it
rejects unsupported powers without output mutation. The split header and added
CatBLST declarations are intentional; no pinned upstream symbol was removed.

Relevant specifications include C99 declaration and enum rules, the platform C
ABI, RFC 9380 for XMD/hash-to-curve, EIP-2333 key derivation, and the BLS12-381
prime-order group definitions recorded in the group-contract ledger.

## Checklist evidence

### 1. Identity and purpose

Full header, include guard, two opaque types, result enum, 253 declarations,
and public-family comments were inspected at the exact reviewed commit.

### 2. Contracts and invariants

Every declaration is covered by an immediately preceding family or per-function
contract. Scalar validity, field representation, point domain, subgroup
requirements, timing, aliasing, scratch, ownership, failure, and output
publication were traced to implementations and supported wrappers.

### 3. Upstream and specification lineage

The pinned upstream symbol-set comparison found 176 retained names, zero
removed names, and 77 CatBLST additions. The only material retained-symbol
signature divergence is the documented Frobenius status result.

### 4. Correctness

Declarations were matched to implementation references. Group-domain contracts
agree with the existing complete group ledger. Scalar fixed-width imports were
corrected to exact, unreduced storage semantics under `AUD-188`.

### 5. Memory and lifetime safety

Distinct context tags prevent accidental context cross-use (`AUD-187`). Wire
aggregate extents are explicit (`AUD-190`). Native-Fr MSM scratch and overlap
requirements are explicit (`AUD-191`). Pairing retains only DST; uniqueness
retains inserted message bytes.

### 6. Undefined and implementation-defined behavior

The header uses C99-compatible declarations and fixed-width types. Opaque
storage requires queried alignment. Partial overlaps and invalid scratch are
outside contract. Enum width/future-value concerns are owned by `AUD-189` and
the `blst_types.h`/`error_codes.h` review; the exact commit enforces four-byte C
storage and forward-safe Rust wrappers.

### 7. Side channels and cryptographic hygiene

Secret-capable multiplication promises constant time only for canonical scalar
inputs at fixed public `nbits`; reduction-required timing is disclosed.
Variable-time paths require public inputs. Field and point conditional negation
are explicitly constant-time. Narrow optimized subgroup contracts are retained
without Release validation overhead.

### 8. Randomness and secret lifecycle

No randomness is generated by this declaration file. Keygen is deterministic.
Secret-key and scalar temporary erasure is an implementation obligation already
tracked in the corresponding implementation reviews; the header does not claim
erasure of caller-owned objects.

### 9. API and ABI

All pinned upstream symbols remain. C++ optional defaults are restored
(`AUD-192`), opaque pointer tags are distinct without layout change, and
Frobenius is nodiscard (`AUD-193`). Raw Rust covers every declaration.

### 10. Concurrency

Stateless arithmetic is reentrant for disjoint caller storage. Pairing and
uniqueness contexts are not internally synchronized. Parallel pairing uses one
context per worker followed by committed merge.

### 11. Resource exhaustion and denial of service

Scratch-size queries check overflow; windows and scalar widths are bounded.
Fixed-table precomputation uses bounded nonthrowing allocation with stack/serial
fallback. XMD bounds output at 8160 bytes. Caller-supplied lengths still govern
public work and readable storage.

### 12. Portability and build behavior

The file is valid for supported C99-or-later and C++ consumers, contains no
compiler intrinsic, and does not expose AVX types. Selected AVX operation is an
implementation choice behind scalar ABI layouts. Apple and non-MSVC packaging
remain milestone evidence rather than declaration defects.

### 13. Performance-sensitive structure

This file adds no runtime abstraction. Narrow contracts preserve direct field,
endomorphism, Pippenger, and pairing paths. Documentation and type-tag changes
generate no instructions. The zero-point scratch fix removes unnecessary
allocation without changing positive-count sizing or the hot kernel.

### 14. Tests and evidence

The 17-case binding-header suite covers opaque type distinction, enum width,
C++ defaults, parsing, and transitive declarations. Existing focused suites
cover scalar conversion, Fp12 Frobenius, native-Fr MSM/scratch capacity,
point-domain behavior, pairing lifecycle, keygen/signing, SHA/XMD, and
aggregate failure atomicity. This record adds no build or test run.

### Specialized additions

The group API contract ledger and public API documentation ledger were used to
cross-check C, selected native C++, generated Rust, and safe Rust terminology.
CUDA, frozen BTx/BTE, Go, benchmarks, and the removed legacy C-over-ABI C++
wrapper are outside this review.

## Findings

- `AUD-187` -- validated: pairing and uniqueness use distinct opaque tags.
- `AUD-188` -- validated: scalar storage and canonicality are distinguished.
- `AUD-190` -- fixed: aggregate wire readable extents are documented.
- `AUD-191` -- fixed: native-Fr MSM failure, timing, overlap, and scratch
  contracts are documented.
- `AUD-192` -- validated: pinned upstream C++ optional defaults are restored.
- `AUD-193` -- validated: fallible Frobenius is nodiscard.
- `AUD-194` -- validated: zero-point Pippenger scratch is zero.
- `AUD-337` -- validated: XMD rejects SHA-256 transcript lengths outside the
  fewer-than-2^64-bit domain before reading input.

Related cross-file `AUD-189` is owned by the public-types/error-code reviews;
its four-byte/forward-safe result representation is present at this commit.
Earlier declaration findings listed in the manifest remain linked to their
own implementation and evidence records.

## Validation

Read-only review commands used for this record:

```powershell
git rev-parse aa4c7ed
git show aa4c7ed:bindings/blst_scalar.h
git show 54e6e55674722fc2797ebb4bbb71b26d881eb4b8:bindings/blst.h
rg -n "blst_[A-Za-z0-9_]+\(" bindings/blst_scalar.h blst bindings/rust/src
git diff 7cc0c0dbde812d2390c34ae50445c3e78676bf7f..aa4c7ed -- bindings/blst_scalar.h
```

Results:

- exact reviewed hash resolved to
  `aa4c7ed187f5a63496f0e94f9284108cbd46c7f1`;
- 253 unique public function declarations were inventoried;
- every declaration has an implementation reference and a raw Rust declaration
  or deliberate matching Rust-local extern;
- all 176 pinned upstream function names remain, with 77 additive CatBLST
  declarations;
- the production diff at the reviewed commit implements the seven findings
  listed above without arithmetic or hot-path changes.

The finding ledger records strict C11/C++20 probes, the 17-case binding-header
suite, targeted Rust layout/unknown-value checks, focused Fp12 tests, and focused
native-Fr/scratch tests as passing. No build, benchmark, sanitizer, CUDA,
portable-AVX Release, or Apple job was launched while creating this record.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Existing repository C/C++ security baseline | Prior implementation findings linked through the manifest | Declaration wording, mathematical domains, ABI source compatibility, and constant-time promises require hand review |
| Clang Static Analyzer / clang-tidy | Existing focused scalar/C-boundary baseline | No new exact-header diagnostic; `AUD-187`--`AUD-194` came from hand review | Cannot prove cross-language ownership, subgroup provenance, or complete C ABI compatibility |
| Clippy | Retained Rust binding compilation and focused wrapper checks | Raw declarations accounted for; related `AUD-189` validated | Raw FFI cannot enforce C pointer, alias, or subgroup preconditions |
| ASan/UBSan | Existing focused boundary/MSM evidence | Implementation evidence remains linked per family | Header-only exact-commit changes do not execute; broad sanitizer milestone remains pending |
| Fuzzing | N/A for declarations | N/A | Decoder, aggregate, XMD, and scratch implementations require their own bounded fuzz targets |

## Residual risk and decision

The declaration file is complete at the reviewed commit: every supported C
declaration has a precise contract and was reconciled with implementation,
selected native C++, raw Rust, and safe Rust where present. `AUD-190` and
`AUD-191` are documentation fixes whose enforcement remains caller-side by
design; the inherited pointer-only ABI cannot validate readable extents or
overlap dynamically without changing the surface or hot paths.

Residual release evidence belongs to underlying implementations: sanitizer and
fuzz expansion, Apple and other supported toolchains, native AVX-512 codegen,
and full milestone matrices. The open pairing-policy decision `AUD-097` remains
an explicit semantic policy item, not a missing declaration contract. These
gaps do not prevent this header's exact file review from being complete.

## Phase 14 field-contract follow-up

The exact header was reconciled again at `314c186` for the reviewed field
families. Fr/Fp/Fp2 shift counts are explicitly public, and the sparse Fp12
operand is documented as packed `{x,y,z}` for `((x,y,0),(0,z,0))` rather than
an ordinary Fp6 multiplier. No signature, ABI layout, symbol, or runtime path
changed.

## Phase 18 XMD-domain follow-up

At `dfd0740`, the public XMD declaration records the complete SHA-256 domain,
including the 17-byte oversize-DST prefix, 64-byte Z_pad, and fixed b0 framing.
Pinned bindgen/Rust 1.65 regeneration is byte-identical and retains exact
253/253 declaration parity. Scalar Debug XMD tests pass 7/7, safe Rust
primitives pass 9/9 in Release, and optimized MSVC inlines the public guard
before hashing without an out-of-line helper call. No signature, ABI layout, or
valid-input digest changed.

## Phase 18 Rust support alias follow-up

At `16e5b49`, the complete G1/G2 projective add declarations state their
existing exact output/input alias support and partial-overlap prohibition. This
lets Rust tests exercise the C contract with raw pointers instead of overlapping
Rust references. Pinned generation remains exact at 253/253 declarations; no
signature, implementation, or code generation changed.
