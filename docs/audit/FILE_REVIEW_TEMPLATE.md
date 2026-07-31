# File Review: `<path>`

## Record

| Field | Value |
| --- | --- |
| Path | `<path>` |
| Roadmap phase / queue | `<phase and dependency position>` |
| Reviewed commit | `<full commit>` |
| Reviewer | `<name>` |
| Date | `<YYYY-MM-DD>` |
| Disposition | `<manifest disposition>` |
| Subsystem | `<subsystem>` |
| Profiles | `<scalar / avx512 / C / C++ / Rust / support>` |
| Generation/provenance | `<handwritten, generated from ..., vendored from ...>` |
| Upstream references | `<commit rows or not_applicable: reason>` |
| Related files | `<headers, implementation, tests, generators, callers>` |

## Purpose and trust boundary

Describe what the file does, what data crosses its boundary, which inputs are
secret, and which guarantees its callers depend on.

## Contracts and invariants

List preconditions, postconditions, valid representations, ranges, alignment,
aliasing, lifetime, ownership, and error behavior.

## Dependency and call-path map

List direct callers and callees. Explain the scalar/AVX-512 relationship and any
public C/C++/Rust exposure.

## Upstream/specification comparison

Record the exact upstream file/commit comparison, applicable standards or
papers, and a rationale for every intentional divergence.

## Checklist evidence

Use the section numbers from `CHECKLIST.md`. Each section must contain concrete
evidence or `N/A: <reason>`.

### 1. Identity and purpose

### 2. Contracts and invariants

### 3. Upstream and specification lineage

### 4. Correctness

### 5. Memory and lifetime safety

### 6. Undefined and implementation-defined behavior

### 7. Side channels and cryptographic hygiene

### 8. Randomness and secret lifecycle

### 9. API and ABI

### 10. Concurrency

### 11. Resource exhaustion and denial of service

### 12. Portability and build behavior

### 13. Performance-sensitive structure

### 14. Tests and evidence

### Specialized additions

## Findings

List finding IDs from `FINDINGS.md`, or `none`.

## Validation

Record exact commands and results. Link logs or artifacts when useful.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | `<query suite and build>` | `<clean / AUD-...>` | `<limits>` |
| Clang Static Analyzer / clang-tidy | `<checker set and compile database>` | `<clean / AUD-...>` | `<limits>` |
| Clippy | `<features and lint groups>` | `<clean / AUD-...>` | `<limits>` |
| ASan/UBSan | `<build and focused tests>` | `<result / AUD-...>` | `<limits>` |
| Fuzzing | `<target, corpus, seed/budget, max input/work>` | `<result / AUD-...>` | `<limits>` |

Map shared subsystem results to this file's exact functions and contracts. A
generic green workflow is not sufficient. Record documented false-positive
references and explain how assembly, FFI, constant-time, mathematical, or
platform-only behavior was covered separately.

## Residual risk and decision

State remaining risk, untested configurations, and why the final status is
appropriate. Do not use "none" unless the evidence supports it.
