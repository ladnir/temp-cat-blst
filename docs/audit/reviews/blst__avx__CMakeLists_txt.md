# File Review: `blst/avx/CMakeLists.txt`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/avx/CMakeLists.txt` |
| Roadmap phase / queue | Phase 19, Packet 19A native AVX compile contract |
| Reviewed commit | `20c075a5000b8f7c29bb06a714112b1c8f0f9742` |
| Reviewer | Codex audit team |
| Date | 2026-08-12 |
| Disposition | `in_scope_production` |
| Subsystem | `avx512` |
| Profiles | native AVX-512/IFMA and developer emulation |
| Generation/provenance | Handwritten local interface target |
| Upstream references | VROOM source pins in `docs/avx/README.md` |
| Related files | `blst/CMakeLists.txt`, AVX vector headers, selected native facades, AVX tests |

## Purpose and trust boundary

The interface target carries compile requirements for header-driven AVX code.
Incorrect flags either fail compilation or silently widen the deployment CPU
contract.

## Contracts and invariants

Ordinary/emulated builds require C++20 and no native ISA flag. Native builds
request AVX2, AVX-512F, AVX-512DQ, and AVX-512IFMA for C++ compilation. The
selected native C++ facade propagates those requirements; scalar-default C/C++
consumers do not.

## Dependency and call-path map

`blst_static` consumes the target privately unless AVX is the selected native
C++ facade. AVX tests explicitly consume it. Vector headers supply the intrinsic
operation surface.

## Upstream/specification comparison

The feature set was derived from the retained intrinsic inventory rather than
the former broad `-march=native` setting.

## Checklist evidence

### 1. Identity and purpose

All four target features/options reviewed.

### 2. Contracts and invariants

Native versus emulated requirements are compile-time-only.

### 3. Upstream and specification lineage

VROOM lineage is recorded by Phase 16.

### 4. Correctness

Structural intrinsic inventory requires AVX2/F/DQ/IFMA; native compilation is
pending supported Clang/GCC execution.

### 5. Memory and lifetime safety

N/A.

### 6. Undefined and implementation-defined behavior

No host autodetection or compiler-default ISA assumption remains.

### 7. Side channels and cryptographic hygiene

Compile flags add no dispatch or runtime path; arithmetic evidence remains in
Phase 16.

### 8. Randomness and secret lifecycle

N/A.

### 9. API and ABI

Only C++ compile requirements propagate; C translation units are excluded by
generator expressions.

### 10. Concurrency

N/A.

### 11. Resource exhaustion and denial of service

N/A.

### 12. Portability and build behavior

Compiler/version and architecture gates live in root CMake.

### 13. Performance-sensitive structure

Exact ISA flags replace `-march=native` without changing implementation or
introducing wrappers.

### 14. Tests and evidence

Preset parsing and scalar no-leakage pass; Peach native execution remains.

### Specialized additions

The retained intrinsic names were inventoried against the requested feature
families.

## Findings

`AUD-340`, with `AUD-342` and `AUD-345` at the root/preset boundary.

## Validation

Local configuration proves no native AVX flags in scalar mode and correct
target relationships. Supported native Clang/GCC compilation is pending.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CMake parser | local MSVC host | passed scalar configuration | MSVC is intentionally unsupported for native AVX |
| Peach/native AVX | deferred | pending | Required before completion |

## Residual risk and decision

Status remains `in_review` until the exact flag set compiles and runs the
focused native AVX target on supported hardware and compiler versions.
