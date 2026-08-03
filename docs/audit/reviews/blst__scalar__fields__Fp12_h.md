# File Review: `blst/scalar/fields/Fp12.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/scalar/fields/Fp12.h` |
| Roadmap phase / queue | Phase 14 scalar extension fields |
| Reviewed commit | `f0a9d9e0c5e04c5e17e549554469baae94bf67ad` |
| Reviewer | Codex audit team |
| Date | 2026-08-03 |
| Disposition | in-scope production |
| Subsystem | scalar Fp12 internal support types |
| Profiles | scalar C++ |
| Generation/provenance | Hand-maintained local wrapper around the public scalar types |
| Upstream references | `54e6e556:src/fp12_tower.c` for the internal Fp4 scratch shape |
| Related files | `bindings/blst_types.h`, `blst/scalar/fields/Fp2.h`, `blst/scalar/fields/Fp12.cpp`, `blst/fields/Fp12.h` |

## Purpose and trust boundary

The header supplies the two-Fp2 `Fp4` scratch aggregate used only by scalar
cyclotomic squaring. Public Fp6/Fp12 representations and methods live in
`bindings/blst_types.h`. It performs no arithmetic, allocation, I/O, or state
management.

## Contracts and invariants

`Fp4` is fixed typed storage for exactly two Fp2 values. Indexing is unchecked
and only internal constant indices zero and one are valid.

## Dependency and call-path map

`Fp12.cpp` is the sole supported production consumer. `blst/fields/Fp12.h`
re-exports the name, although no supported non-Fp12 production caller uses it.

## Upstream/specification comparison

The storage mirrors upstream `typedef vec384x vec384fp4[2]`. The local header
lacks the normal license/provenance header and contains orphaned schedule
comments; both are tracked by `AUD-087`.

## Checklist evidence

### 1. Identity and purpose

All 40 lines and every active reference were inspected.

### 2. Contracts and invariants

Only fixed internal indices are used. No external length or pointer enters.

### 3. Upstream and specification lineage

The Fp4 scratch layout matches pinned upstream; provenance text is missing.

### 4. Correctness

The aggregate and accessors are structurally correct. Orphaned comments do not
describe any declaration.

### 5. Memory and lifetime safety

Typed automatic storage only; no union, allocation, or retained reference.

### 6. Undefined and implementation-defined behavior

Unchecked indexing has the same internal precondition as a raw two-element
array and is only called with compile-time-valid indices.

### 7. Side channels and cryptographic hygiene

N/A: storage/accessors only, with constant indices.

### 8. Randomness and secret lifecycle

N/A: no randomness or ownership.

### 9. API and ABI

This should remain an internal implementation type; it is not part of the C or
Rust ABI.

### 10. Concurrency

N/A: no state.

### 11. Resource exhaustion and denial of service

N/A: fixed 192-byte automatic aggregate.

### 12. Portability and build behavior

Uses only standard C++ aggregate/reference operations plus project host/device
attributes.

### 13. Performance-sensitive structure

The fixed two-element layout is appropriate for hot cyclotomic squaring.

### 14. Tests and evidence

Exercised transitively by scalar Fp12 cyclotomic-square and final-exp tests.

### Specialized additions

N/A.

## Findings

`AUD-087`.

## Validation

Repository-wide reference search and strict host parsing through `Fp12.cpp`.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Existing C/C++ baseline | No alert | Header hygiene is review-driven |
| Clang Static Analyzer / clang-tidy | Not rerun | N/A | No executable logic |
| Clippy | N/A | N/A | No Rust implementation |
| ASan/UBSan | Transitively covered | No header-specific failure | Does not assess provenance |
| Fuzzing | N/A | N/A | Storage-only header |

## Residual risk and decision

The file remains in review pending provenance/hygiene cleanup and a decision
whether the scratch type should remain re-exported by the selected field header.
