# File Review: `blst/exports.cpp`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/exports.cpp` |
| Roadmap phase / queue | Phase 18 public C boundary, Wave 1 |
| Reviewed commit | `02f5e5ff078b251ca3b5f86f4f83ea64af664bc7` (implementation); public declaration/header companion `aa4c7ed187f5a63496f0e94f9284108cbd46c7f1` |
| Reviewer | Codex audit team |
| Date | 2026-08-06 |
| Disposition | in-scope production, complete |
| Subsystem | scalar, Fr/Fp/Fp2, SHA-256, and diagnostic conversion C exports |
| Profiles | scalar, C, C++, Rust; selected AVX callers use the same scalar storage contracts |
| Generation/provenance | Hand-maintained C++ translation of pinned upstream `src/exports.c` |
| Upstream references | BLST `54e6e55674722fc2797ebb4bbb71b26d881eb4b8:src/exports.c` |
| Related files | `bindings/blst_scalar.h`, `bindings/blst_types.h`, `blst/fields/Fr.h`, `Fp.h`, `Fp2.h`, `blst/sign/sha256.h`, `blst/exports.h`, `test/audit_field_boundaries.cpp`, `test/BindingHeader_tests.cpp` |

## Purpose and trust boundary

This file exposes direct C-callable scalar and base-field arithmetic,
serialization, validity checks, Fp2 operations, short SHA-256, and retained hex
test helpers. Callers provide all object and byte storage. Imported bytes and
scalars are untrusted; secret-key operations can consume secrets. The boundary
must preserve exact public representations, alias rules, constant-time behavior
where promised, failure publication, and C exception safety without inserting
copies or dynamic dispatch into field hot paths.

## Contracts and invariants

- `blst_scalar` is arbitrary 32-byte little-endian integer storage. It is not
  implicitly canonical. `blst_scalar_fr_check` accepts `[0,r)`; `blst_sk_check`
  accepts `[1,r)`.
- `blst_fr` stores an Fr value in the public representation documented by the
  header. `blst_fr_from_scalar` requires a scalar in `[0,r)` and does not
  validate or reduce in Release. The reverse conversion publishes canonical
  scalar storage.
- Exact scalar/Fr output-input aliasing is supported. At exact alias, the live
  source is copied to a typed local, the destination object's lifetime is
  placement-started, and the local is wiped after publication. Distinct paths
  remain direct. Partial overlap is unsupported.
- Field inputs must be valid representations unless a decoder/check operation
  states otherwise. Import/export byte lengths and endianness are fixed by the
  declaration. Zero-bit Fr exponentiation returns one without reading an
  exponent pointer.
- `blst_fr_ct_bfly` and `blst_fr_gs_bfly` require distinct `x0` and `x1`;
  `twiddle` may exactly alias either. Constant-time conditional negation does
  not branch on the flag.
- Fp square root clears output on failure. The accepted `AUD-084` Fp2 behavior
  publishes an unspecified candidate on `false`; callers must discard it, and
  exact input/output alias may overwrite the input.
- Fallible Frobenius mapping leaves output unchanged for unsupported powers and
  is `BLST_NODISCARD` at the public declaration.
- SHA-256 supports its documented pointer/length pairs and length domain; it
  owns no storage after return. Hex helpers are retained test facilitators, not
  parsers for adversarial production input.
- Functions allocate no heap memory, throw no exception, retain no caller
  pointer, and use no runtime backend selection. No exception path crosses C.
- Authoritative C and C++-callable declarations are in
  `bindings/blst_scalar.h`. Core C++ field types and operations are documented
  in their owning `blst/fields` headers. Raw Rust declarations are generated
  in `bindings/rust/src/bindings.rs`; safe Rust field/scalar wrappers document
  validity and result behavior in `primitives.rs`. The declaration ledger marks
  these supported surfaces documented.

## Dependency and call-path map

C/C++ consumers and Rust wrappers call these exports. Arithmetic forwards
directly to scalar `Fr`, `Fp`, and `Fp2` static operations and backend assembly
through those owning types. Conversion uses common byte helpers. SHA-256 calls
the reviewed internal short-message interface. There is no dependency on the
selected AVX vector type, no runtime dispatch, no allocation, and no mutable
global state.

`blst/exports.h` is not the public declaration surface. Experimental CUDA still
includes that path, so `AUD-195` reduced it to a documented empty compatibility
shim rather than deleting it. Supported host declarations remain solely in the
binding header.

## Upstream/specification comparison

The direct wrapper structure and symbols follow pinned upstream
`src/exports.c`. CatBLST expresses public/native type identity in C++, replaces
C union-punning and lifetime assumptions with standard C++20 operations,
defines boundary cases, narrows secret erasure, and documents aliases and
failure output. Field arithmetic remains BLS12-381-specific; SHA-256 behavior
is compared to FIPS 180-4 vectors. The Fp2 square-root failure candidate is an
intentional documented compatibility behavior under `AUD-084`.

## Checklist evidence

### 1. Identity and purpose

Every definition, declaration, direct callee, representation conversion,
alias path, error result, include, and history entry was inspected. The stale
neighboring `exports.h` declarations are no longer mistaken for a public API.

### 2. Contracts and invariants

Scalar/Fr domains, field validity, endianness, exponent widths, exact aliasing,
butterfly overlap, conditional flag behavior, sqrt/Frobenius failure
publication, and SHA pointer/length contracts are documented across C,
C++-callable, and Rust surfaces.

### 3. Upstream and specification lineage

The file is compared against pinned `src/exports.c`; intentional C++ lifetime,
fixed-width ABI, documentation, and boundary-hardening divergences are linked
to findings. Underlying field proofs and assembly bounds remain recorded in
their owning file reviews.

### 4. Correctness

Focused field tests cover scalar import/export and validity, Fr/Fp/Fp2
arithmetic, zero exponent, exact aliases, sqrt success/failure, Frobenius
domains, and SHA known answers. Exact scalar/Fr alias results match distinct
storage and preserve canonical conversion semantics.

### 5. Memory and lifetime safety

`AUD-070` removed invalid union/whole-object punning from field boundaries. At
commit `02f5e5f`, exact scalar-to-Fr and Fr-to-scalar conversion now captures
the source, placement-starts the destination type, publishes the result, and
wipes the temporary. Distinct paths add no copy. Every routine remains
caller-storage-only with no retained pointer or allocation.

### 6. Undefined and implementation-defined behavior

Public/native exact types, fixed-width booleans and integers, size/alignment,
and standard-layout/trivial-copy properties are checked by the header/type
packet. Zero exponent avoids underflow and null reads. Shift counts and import
sizes are bounded by owning routines. Exact alias is handled as an object-
lifetime transition, not type punning; partial overlap is excluded.

### 7. Side channels and cryptographic hygiene

Field arithmetic and `cneg` preserve constant-time selection where promised.
Validity checks and Euclidean inverse are documented according to their actual
timing. No data-dependent diagnostic output exists on production paths. The
boundary adds no heap, dispatch, or representation conversion to hot same-type
arithmetic wrappers.

### 8. Randomness and secret lifecycle

N/A for randomness: this file has no RNG. Narrow secret-derived temporaries in
scalar/Fr exact-alias conversion are erased through the retained compiler
barrier. Generic field outputs are not indiscriminately zeroed, avoiding
unnecessary hot-path overhead.

### 9. API and ABI

Definitions use the exact public `blst_*` types. C and C++ declarations share
one authoritative header, while Rust uses generated fixed-layout declarations.
Enum and opaque-context changes in the companion packet do not alter these
symbols. `exports.h` contains no competing declarations. `BLST_NODISCARD`
protects the fallible Frobenius result.

### 10. Concurrency

All state is local or caller-owned. Constants are immutable. There is no lazy
initialization, shared mutable cache, thread pool, or callback.

### 11. Resource exhaustion and denial of service

All operations use fixed object sizes. SHA work is linear in caller-supplied
length and does not allocate. There are no input-sized stack arrays or heap
paths in this file.

### 12. Portability and build behavior

The implementation uses C++20 object-lifetime tools and fixed-width public
types on the supported 64-bit profiles. Standalone C11 and C++20 headers and
MSVC/GCC/Clang strict compilation cover the boundary. CUDA is experimental and
its compatibility include was not audited as a production implementation.

### 13. Performance-sensitive structure

Same-type arithmetic wrappers are one direct static call and retain upstream
assembly routing. The lifetime-safe scalar/Fr local is taken only for exact
aliasing; the common distinct case stays direct. No runtime polymorphism,
allocation, conversion layer, or generalized temporary was introduced.

### 14. Tests and evidence

Existing focused field-boundary tests pass the exact-alias, zero, boundary,
sqrt, and conversion cases under supported scalar configurations. Strict
header/type tests pass at `aa4c7ed`. The host public/scalar build passed after
the `exports.h` cleanup. This documentation-only review ran no build or test.

### Specialized additions

Optimized wrapper inspection previously confirmed Fr/Fp/Fp2 calls contain no
boundary copy or conversion on normal paths. Constant-time and assembly bounds
are reviewed in the scalar field/backend records, not inferred from wrapper
source alone.

## Findings

`AUD-041`, `AUD-063`, `AUD-069` through `AUD-073`, `AUD-084`, `AUD-141`,
`AUD-142`, `AUD-144`, and `AUD-195`. The final exact-alias lifetime repair is
part of validated `AUD-070`; the neighboring stale-header cleanup is
`AUD-195`.

## Validation

- Commit `02f5e5ff078b251ca3b5f86f4f83ea64af664bc7`: focused scalar field
  boundary coverage passed, including exact scalar/Fr aliasing in both
  directions; structural inspection confirms destination lifetime starts
  precede typed writes.
- Prior MSVC/GCC/Clang focused translation units and the 18-case scalar boundary
  probe passed for `AUD-070`; optimized reciprocal and ordinary field wrappers
  retained their recorded frames/call shapes.
- Focused Fp square-root residue and exact-aliased nonresidue tests passed;
  `AUD-084` remains an explicitly documented accepted result contract for Fp2.
- Retained SHA-256 known-answer and boundary tests passed under the recorded
  scalar packet.
- Commit `aa4c7ed187f5a63496f0e94f9284108cbd46c7f1`: strict C11/C++20 declaration
  probes and the 17-case binding-header suite passed; the host public-header and
  scalar build passed with `exports.h` reduced to its compatibility shim.
- No command was executed for this documentation-only record.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Existing scalar C-boundary baseline | `AUD-041`; stale-comment hygiene closed as `AUD-195` | Exact C++ lifetime and constant-time code generation required manual inspection |
| Clang Static Analyzer / clang-tidy | Existing MSVC/GCC/Clang strict baselines | No remaining mapped blocker | Algebraic correctness and assembly bounds are outside generic static analysis |
| Clippy | Existing supported Rust binding runs | No wrapper blocker mapped here | Does not analyze native implementation or raw FFI object lifetime |
| ASan/UBSan | Existing scalar boundary baselines; not rerun for final packet | No mapped runtime report | Exact-alias lifetime transition was validated structurally and by focused behavior |
| Fuzzing | Planned import/serialization campaign | Coverage gap | Hex test helpers and malformed byte domains are not yet comprehensively fuzzed |

## Residual risk and decision

The file is complete at the recorded implementation commit with its public
contracts reviewed at the companion commit. Remaining release work is the
broader sanitizer/fuzz/coverage matrix and cross-toolchain milestone evidence.
The accepted Fp2 square-root failed-output candidate remains intentionally
documented under `AUD-084`. Experimental CUDA's use of the empty compatibility
header is outside production readiness; no supported host declaration depends
on it.
