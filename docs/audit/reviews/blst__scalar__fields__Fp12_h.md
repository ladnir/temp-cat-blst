# File Review: `blst/scalar/fields/Fp12.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/scalar/fields/Fp12.h` |
| Roadmap phase / queue | Phase 14 scalar extension fields |
| Reviewed commit | `314c186149b9cc95d85a0c54e2dd651f137e4236` |
| Reviewer | Codex audit team |
| Date | 2026-08-03 |
| Disposition | in-scope production |
| Subsystem | scalar Fp6/Fp12 types and internal Fp4 support |
| Profiles | scalar C, C++, and Rust shared representation |
| Generation/provenance | Hand-maintained native C++ type declarations with upstream tower lineage |
| Upstream references | `54e6e556:src/fp12_tower.c` for the internal Fp4 scratch shape |
| Related files | `bindings/blst_types.h`, `blst/scalar/fields/Fp2.h`, `blst/scalar/fields/Fp12.cpp`, `blst/fields/Fp12.h` |

## Purpose and trust boundary

The header owns the native C++ Fp6 and Fp12 representations, GLV schedule
types, method contracts, typed Fp12 identity, and the two-Fp2 `Fp4` scratch
aggregate used by scalar cyclotomic squaring. The C++ branch of
`bindings/blst_types.h` aliases these exact types; it does not redefine them.
The header performs no dynamic allocation or mutable global initialization.

## Contracts and invariants

Fp4, Fp6, and Fp12 are fixed typed aggregates. Indexing is unchecked and valid
indices are a caller precondition. GLV schedules have fixed storage, while
their currently mutable length contracts remain covered by `AUD-083`.
`Fp12::from_bendian` is `[[nodiscard]]`, requires 576 readable bytes, rejects
noncanonical coefficients, permits input/output overlap through temporary
decoding, and leaves output unchanged on failure.

## Dependency and call-path map

`Fp12.cpp`, pairing/GT, FFT, and the C boundary consume these declarations.
`blst/fields/Fp12.h` re-exports the native names, and `bindings/blst_types.h`
maps the C++ `blst_fp6`/`blst_fp12` names directly to them.

## Upstream/specification comparison

The Fp4 storage mirrors upstream `typedef vec384x vec384fp4[2]`; Fp6/Fp12
retain the tower layout used by the C ABI. Remaining provenance and schedule
hygiene work is tracked by `AUD-087`.

## Checklist evidence

### 1. Identity and purpose

The original header and the declarations formerly housed in
`bindings/blst_types.h` were inspected at the recorded review commit. Their
ownership move is mechanical; the open Fp12 API findings remain attached here.

### 2. Contracts and invariants

Only fixed internal indices are used. No external length or pointer enters.

### 3. Upstream and specification lineage

The Fp4 scratch layout matches pinned upstream. Apache/SPDX provenance is now
present.

### 4. Correctness

The aggregate layouts and accessors are structurally correct and compile
against every out-of-line definition. The GT scalar domain is now explicit and
the obsolete missing GLV declaration has been removed.

### 5. Memory and lifetime safety

Typed aggregate storage only; no union, allocation, or retained reference.
Executors that trust caller-mutable schedule lengths remain covered by
`AUD-083`.

### 6. Undefined and implementation-defined behavior

Unchecked indexing has the same internal precondition as a raw two-element
array and is only called with compile-time-valid indices.

### 7. Side channels and cryptographic hygiene

The declarations distinguish general, GT-specific, constant-time, and
explicitly variable-time GLV operations. `AUD-082` is resolved by the safe
unsuffixed default and named public-scalar fast path.

### 8. Randomness and secret lifecycle

N/A: no randomness or ownership.

### 9. API and ABI

The C++ `blst_fp6`/`blst_fp12` names are exact aliases of these native types.
C and Rust expose the matching POD/`repr(C)` storage, so layout is ABI-relevant
even though member functions are visible only in C++.

### 10. Concurrency

N/A: no state.

### 11. Resource exhaustion and denial of service

All represented objects and schedules have fixed storage. Schedule executors
have an accepted low-level precondition that their input is unmodified compiler
output (`AUD-083`).

### 12. Portability and build behavior

Uses only standard C++ aggregate/reference operations plus project host/device
attributes.

### 13. Performance-sensitive structure

The fixed tower layouts are appropriate for hot pairing arithmetic. The move
adds no wrapper, virtual dispatch, allocation, or representation conversion.

### 14. Tests and evidence

Exercised transitively by scalar Fp12 cyclotomic-square and final-exp tests.

### Specialized additions

N/A.

## Findings

`AUD-081`, `AUD-082`, `AUD-083`, `AUD-085`, `AUD-086`, `AUD-087`, and
`AUD-088`, `AUD-198`, `AUD-204`, `AUD-205`, and `AUD-207`.

## Validation

Repository-wide reference search and strict host parsing through `Fp12.cpp`.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Existing C/C++ baseline | No alert | Header hygiene is review-driven |
| Clang Static Analyzer / clang-tidy | Not rerun | N/A | No executable logic |
| Clippy | Retained raw binding baseline | No header-specific diagnostic | Rust shares the C storage layout but not this C++ implementation |
| ASan/UBSan | Transitively covered | No header-specific failure | Does not assess provenance |
| Fuzzing | N/A | N/A | Storage-only header |

## Residual risk and decision

### Post-review remediation

The header now identifies `cyclotomic_exp` and all GLV operations as GT-only,
describes schedules as variable-time trusted bytecode, restores Apache/SPDX
provenance, and removes the obsolete undefined `glv_split_4` duplicate. These
changes resolve the header portions of `AUD-083`, `AUD-086`, and `AUD-087`.

Full-width GT scalar behavior is now explicit and covered by focused comparison
against basic exponentiation. The GT timing decision is implemented, selected,
and benchmarked under `AUD-082`/`AUD-088`. Canonical, failure-atomic Fp12
decoding resolves this header's part of `AUD-085`. The file remains in review
pending a decision whether the Fp4 scratch type should remain re-exported by the
selected header.

## Wave 2 exact-contract follow-up

The complete declaration surface was reread at `314c186`. General Fp12
exponentiation and raw GLV schedule helpers now say when exponent data is
public and variable-time; CT GT helpers state the canonical-input guarantee and
the permitted reduction-required timing distinction. `in_group` is documented
as nonzero prime-order GT membership, while `is_cyclotomic` is the weaker
relation that also accepts zero. Sparse `{x,y,z}` storage is defined as Fp12
`((x,y,0),(0,z,0))`, not an ordinary Fp6 operand.
