# File Review: `blst/scalar/fields/Fp2.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/scalar/fields/Fp2.h` |
| Roadmap phase / queue | Phase 14 scalar extension fields, after Fp representation/lifetime remediation |
| Reviewed commit | `0ee24de0da80b256979832ea872272f635599853` |
| Reviewer | Codex audit team |
| Date | 2026-08-03 |
| Disposition | in-scope production |
| Subsystem | scalar quadratic extension field |
| Profiles | scalar, C, C++, Rust; experimental CUDA definitions are not release evidence |
| Generation/provenance | Hand-maintained C++ refactor with upstream BLST `src/fp12_tower.c` and `src/exp.c` lineage |
| Upstream references | Pinned upstream head `54e6e55674722fc2797ebb4bbb71b26d881eb4b8` |
| Related files | `bindings/blst_types.h`, `blst/scalar/fields/Fp.h`, `blst/scalar/fields/Fp12.cpp`, `blst/exports.cpp`, `test/Fp2_tests.cpp`, `test/C_api_tests.cpp` |

## Purpose and trust boundary

This header implements Montgomery-form arithmetic in
`Fp2 = Fp[u]/(u^2+1)`, including sign, square-root, inverse, Frobenius, and
the low-level scalar/assembly selection used by G2, Fp12, pairing, and public
C/Rust field operations. Exponents passed to `exp` are public; values passed
through ordinary arithmetic may be secret. No operation allocates or throws.

## Contracts and invariants

- Each value has two reduced Montgomery Fp coefficients.
- Full-object exact aliasing is supported by the selected add, subtract,
  multiply, square, negate, conjugate, inverse, and square-root paths.
- `reciprocal` is variable-time; `flt_reciprocal` is the fixed-chain choice.
- Nonzero public exponents declare their most-significant bit; zero bits must
  produce one.
- Internal square-root functions deliberately produce candidates on failure;
  C fallible wrappers must not publish those candidates.

## Dependency and call-path map

Fp2 calls fixed-width Fp operations and selected scalar assembly primitives.
Fp12 and scalar G2 consume it directly. This header owns the native C++ type;
the C++ binding name `blst_fp2` is an exact alias, while C sees the matching
POD layout. `blst/exports.cpp` exposes retained C and Rust operations without
a representation copy. There is no allocation, callback, global mutation,
runtime dispatch, or exception path.

## Upstream/specification comparison

The multiply, square, inverse, square-root, sign, and Frobenius formulas retain
the pinned upstream tower arithmetic. CatBLST adds native C++ types, supported
public methods, and host/device selection. The inherited exponent helper
assumed a nonzero declared top bit; unlike repaired Fp/Fr, its zero boundary
was not carried forward. RFC 9380 governs the extension-field sign convention.

## Checklist evidence

### 1. Identity and purpose

Every line and direct supported caller was inspected at the recorded commit.
The file is handwritten/derived, not generated.

### 2. Contracts and invariants

Zero/one, representation, aliasing, exponent, inverse, square-root, and sign
contracts were traced. `AUD-079` and `AUD-084` record violated boundaries.

### 3. Upstream and specification lineage

Arithmetic was compared with upstream `src/fp12_tower.c` and `src/exp.c`.
The Apache-2.0/SPDX header is present.

### 4. Correctness

Fixed-width component arithmetic and temporary ordering are alias-safe. The
zero-bit exponent underflow is a direct correctness and bounds defect.

### 5. Memory and lifetime safety

Storage is fixed and typed. There is no heap or retained pointer. The exponent
underflow can read beyond caller storage until remediated.

### 6. Undefined and implementation-defined behavior

No active union punning remains. Fixed limb widths and little-endian support
are enforced globally. Macro leakage and dead definitions are in `AUD-087`.

### 7. Side channels and cryptographic hygiene

Fp2 arithmetic and fixed-chain square-root paths are branchless with respect
to field values. `exp` is a public-exponent routine. Euclidean reciprocal is
correctly documented variable-time.

### 8. Randomness and secret lifecycle

N/A: no randomness or secret ownership. General temporaries contain field
values and are not boundary-created long-lived secret copies.

### 9. API and ABI

The shared C/C++ representation is exact. Transitive C operations are
allocation-free and nonthrowing. The Fp2 square-root wrapper violates
failure-atomic publication until `AUD-084` is fixed.

### 10. Concurrency

Host constants are immutable internal-linkage values. No lazy or mutable
global state exists.

### 11. Resource exhaustion and denial of service

Work and stack are fixed except for caller-declared public exponent length.
The underflowed zero length is the only identified unbounded path.

### 12. Portability and build behavior

The host path parses under MSVC, GCC, and Clang through the Fp12 translation
unit. Experimental CUDA is outside this release review.

### 13. Performance-sensitive structure

The selected assembly calls, widened Fp12 consumers, fixed temporaries, and
in-place operation ordering should remain intact. Approved boundary fixes need
no heap, virtual dispatch, or runtime provider selection.

### 14. Tests and evidence

Existing tests cover basic arithmetic, inverse, successful square root,
Frobenius, and zero shifts. They omit zero exponent and failed public C square
root; both gaps map to open findings.

### Specialized additions

N/A: this is scalar extension-field code. Its assembly callees are covered by
their generator reviews; experimental CUDA is not a release gate.

## Findings

`AUD-079`, `AUD-084`, `AUD-087`.

## Validation

- Full source/caller/test structural review at `f0a9d9e`.
- MSVC 19.50 strict syntax through `Fp12.cpp`; local downstream warnings are
  recorded in `AUD-087`.
- Clang 18 C++20 strict syntax through `Fp12.cpp` completed with warnings
  dispositioned to existing common-file reviews or `AUD-087`.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Existing C/C++ baseline | No separate new alert | Does not prove exponent preconditions, constant time, or arithmetic |
| Clang Static Analyzer / clang-tidy | Not rerun for this packet | Pending remediation | Source review and strict compiler diagnostics only |
| Clippy | N/A | N/A | Rust only calls the C ABI |
| ASan/UBSan | Prior focused scalar baseline | Does not cover zero exponent | Regression required after fix |
| Fuzzing | Not yet present | Coverage gap | Fp2 arithmetic/oracle target belongs to the later fuzz phase |

## Residual risk and decision

### Post-review remediation

`AUD-079` is fixed: zero exponent bits return one without reading the exponent,
nonzero inputs retain the declared-top-bit Debug assertion, and reverse
iteration is width-safe. `AUD-084` is accepted with an explicit C contract:
failed square-root output is unspecified, must be ignored, and exact aliasing
may overwrite the input. The unused generic macro from `AUD-087` was removed.
The focused zero-bit regression and strict MSVC test-unit compilation pass.

The file remains in review for independent arithmetic-oracle coverage,
sanitizers, and retained Rust/Apple evidence.

## Generated-addition-chain follow-up

`AUD-182` records that fixed Fp2 inversion deliberately maps zero to zero via
the fixed Fp p-2 chain, which branchless G2 infinity normalization relies on.
The declaration now states this postcondition, the redundant direct chain
include is removed, and focused scalar Debug/Release tests cover the Fp2 zero
case. The implementation is recorded at
`529a79377b6ef848ec90fce35e7580a9f451fc4f`; the file's broader residual
evidence remains pending. No generic arithmetic scratch erasure was added.
