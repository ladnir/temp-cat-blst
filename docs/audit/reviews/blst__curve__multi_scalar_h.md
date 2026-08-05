# File Review: `blst/curve/multi_scalar.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/curve/multi_scalar.h` |
| Roadmap phase / queue | Phase 15 selected scalar MSM facade |
| Reviewed commit | pending remediation commit |
| Reviewer | Codex audit team |
| Date | 2026-08-04 |
| Disposition | in-scope production |
| Subsystem | selected native C++ MSM surface |
| Profiles | native C++ scalar; experimental CUDA compatibility conditional |
| Generation/provenance | Local selected-backend facade |
| Upstream references | not_applicable: CatBLST selected native facade |
| Related files | scalar MSM header/implementation, C declarations, Rust multipoint wrapper |

## Purpose and trust boundary

The file selects the native scalar MSM surface. It contains no algorithm or
state; its security role is preventing internal/debug helpers from becoming an
accidental supported API and publishing the common contracts of retained
operations.

## Contracts and invariants

The header documents subgroup/public-scalar, 256-bit, window, scratch,
compact-descriptor, exact-array, and empty-input rules. Only fixed-table and
Pippenger operations plus their sizing helpers are exported. CUDA-only helper
aliases remain conditional and are outside the release gate.

## Dependency and call-path map

Native C++ consumers include this wrapper, which includes the scalar header.
C uses explicit wrappers instead. Safe Rust uses C FFI. No runtime dispatch or
backend conversion occurs here.

## Upstream/specification comparison

N/A as a local facade. The surface was compared with actual supported callers;
dead allocator dispatch and ordinary-build debug hash exposure were removed.

## Checklist evidence

### 1. Identity and purpose

Every alias and consumer was inspected.

### 2. Contracts and invariants

The complete retained family contract is adjacent to the aliases.

### 3. Upstream and specification lineage

N/A: local facade; underlying arithmetic retains upstream lineage.

### 4. Correctness

Scalar native/C/Rust focused results validate the selected functions.

### 5. Memory and lifetime safety

No storage or object access is implemented here.

### 6. Undefined and implementation-defined behavior

No casts, arithmetic, or layout assumptions are present.

### 7. Side channels and cryptographic hygiene

The facade identifies MSM as variable-time/public-data. Debug hashing is not
exposed in ordinary builds.

### 8. Randomness and secret lifecycle

N/A.

### 9. API and ABI

Only reviewed operations are aliased. Removed helpers were unsupported native
internals; public C ABI is unchanged.

### 10. Concurrency

No state.

### 11. Resource exhaustion and denial of service

Documented sizing helpers and caps belong to the underlying implementation.

### 12. Portability and build behavior

Repeated/standalone selected header compilation passes in scalar Debug/Release.

### 13. Performance-sensitive structure

Aliases are compile-time only and introduce no code or conversion.

### 14. Tests and evidence

Covered by the focused scalar and Rust packet.

### Specialized additions

Experimental CUDA aliases are preserved behind `BLST_CUDA` without release
claims. No CUDA build ran.

## Findings

`AUD-119`, `AUD-120`.

## Validation

Structural scans and MSVC Debug/Release compilation pass; normal builds expose
no dead dispatcher or hash alias.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Existing baseline | Header guard and stale-debug findings resolved | Alias contract required hand review |
| Clang Static Analyzer / clang-tidy | Existing baseline | No blocker | No implementation |
| Clippy | N/A | Rust wrapper reviewed separately | N/A |
| ASan/UBSan | N/A | No runtime code | N/A |
| Fuzzing | N/A | No runtime code | N/A |

## Residual risk and decision

Complete pending the final remediation commit hash. CUDA compatibility aliases
remain explicitly experimental.
