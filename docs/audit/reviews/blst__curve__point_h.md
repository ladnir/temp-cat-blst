# File Review: `blst/curve/point.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/curve/point.h` |
| Roadmap phase / queue | Phase 15 selected native C++ facades |
| Reviewed commit | `5d18e42ad96dbba5d2ce6fbda8af0ce160520e85` |
| Reviewer | Codex audit team |
| Date | 2026-08-06 |
| Disposition | in-scope production |
| Subsystem | selected native C++ G1/G2 point facade |
| Profiles | scalar, AVX-512, portable AVX Debug, native C++ |
| Generation/provenance | Hand-maintained CatBLST selected-backend facade |
| Upstream references | not_applicable: CatBLST selected native C++ facade |
| Related files | `blst/scalar/curve/point.h`, `blst/avx/curve/point.h`, `bindings/blst_scalar.h`, `test/avx/Avx_c_api_contract_tests.cpp` |

## Purpose and trust boundary

The file chooses the named native C++ G1/G2 types for the configured default
backend. It owns no arithmetic. Callers rely on the selected type having the
same public operation names and contracts in scalar- and AVX-default builds.
The generic templates, XYZZ types, and PointZZ aliases intentionally remain
visible scalar implementation types and are not described as selected types.

## Contracts and invariants

- `G1Point`, `G1AffinePoint`, `G2Point`, and `G2AffinePoint` select the default
  scalar or AVX implementation without encoding validation state.
- Decoding establishes canonical on-curve values, not subgroup membership.
- Basic point primitives accept on-curve points. Optimized multiplication and
  advanced algorithms require previously validated subgroup points.
- The inherited `unchecked_mult` exceptions apply an exact scalar to an
  on-curve point without requiring subgroup membership.
- The concrete selected class declaration is authoritative for timing,
  aliasing, byte sizes, and failure behavior.
- Generic point templates and implementation helpers are deliberately visible;
  this file does not promise backend selection for them.

The C contracts remain in `bindings/blst_scalar.h`; safe Rust uses its own
documented wrappers, and generated raw Rust inherits the C contract. This file
owns only the selected native C++ spelling.

## Dependency and call-path map

`BLST_DEFAULT_AVX` selects either `blst/avx/curve/point.h` or
`blst/scalar/curve/point.h`. Native C++ consumers include this facade directly
or transitively through other selected curve headers. Frozen BTx/BTE code uses
the retained generic scalar aliases; no BTx/BTE source was changed or audited.

## Upstream/specification comparison

This facade is CatBLST-specific. The public C names and the two exact-scalar
exceptions retain upstream BLST compatibility. Backend selection and native
C++ documentation are local policy.

## Checklist evidence

### 1. Identity and purpose

Every declaration and direct repository consumer was inspected. The facade
contains aliases and selected helper exports only.

### 2. Contracts and invariants

The file now distinguishes four selected named types from intentionally scalar
generic aliases. Concrete backend declarations carry the operation contracts.

### 3. Upstream and specification lineage

N/A: the selected C++ facade is local. Compatible operation spellings are
traced to the scalar and AVX declarations.

### 4. Correctness

Scalar selected types already supplied `unchecked_mult`. `AUD-164` adds a true
AVX fixed-window implementation so the selected surface no longer depends on
the configured backend.

### 5. Memory and lifetime safety

Aliases preserve the owning type's layout and lifetime. No reinterpretation,
allocation, ownership, or boundary copying occurs in this facade.

### 6. Undefined and implementation-defined behavior

`AUD-168` replaces the reserved double-underscore include guard. Supported
builds remain 64-bit little-endian as recorded by the platform audit.

### 7. Side channels and cryptographic hygiene

The facade adds no code. The AVX exact-scalar implementation uses fixed-window
Booth recoding, complete point operations, and full constant-time table scans
for a fixed public bit count; it does not use scalar-point conversions or
subgroup endomorphisms.

### 8. Randomness and secret lifecycle

N/A: no randomness or storage is owned here. Concrete multiplication timing
and secret-key preconditions are documented by the selected classes.

### 9. API and ABI

The named selected API has operation parity for the two inherited
`unchecked_mult` methods. Generic/internal visibility is an intentional design
decision under `AUD-167`, not a selected-backend promise.

### 10. Concurrency

N/A: aliases own no mutable state or execution.

### 11. Resource exhaustion and denial of service

N/A: aliases allocate no resources. The added exact-scalar implementation uses
fixed stack tables of 8 or 16 points and accepts at most 256 bits.

### 12. Portability and build behavior

The facade compiles in the focused portable-AVX Debug configuration. Native
AVX runtime and other milestone platforms remain Phase 16/release evidence.

### 13. Performance-sensitive structure

Backend choice remains compile-time. Native AVX operations stay in AVX point
storage; the exact-scalar exception introduces no runtime dispatch, allocation,
GLV/GLS decomposition, or scalar point conversion.

### 14. Tests and evidence

Portable-AVX Debug passes all 8 focused contract tests. Exact multiplication is
compared with the audited scalar oracle at 1, 128, 175, 176, and 256 bits and
with order-11 G1/order-13 G2 fixtures; zero bits and exact output aliasing are
also covered. Raw G1/G2 decode failures leave sentinels unchanged.

### Specialized additions

Low-order fixtures force the complete addition path during table construction
and accumulation, preventing recurrence of the upstream low-order fixed-window
defect recorded by `AUD-128`.

## Findings

`AUD-164`, `AUD-165`, `AUD-166`, `AUD-167`, `AUD-168`.

## Validation

```text
cmake --build out/build/audit-phase1-avxportable --target blst_avx_contract_tests --config Debug -j 2
out/build/audit-phase1-avxportable/test/Debug/blst_avx_contract_tests.exe
```

Result: build passed; 8/8 tests passed. Structural inspection confirms the AVX
exact-scalar implementation contains no scalar-point conversion or GLV/GLS
call and scans every fixed-window table entry.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Prior supported baseline | `AUD-165` was found by hand review | Header-only AVX timing and mathematics require direct review |
| Clang Static Analyzer / clang-tidy | Prior baseline | No packet-specific result | Does not prove group arithmetic or constant time |
| Clippy | N/A | N/A | Native C++ facade |
| ASan/UBSan | Not run in this focused packet | milestone evidence pending | User-gated broad run |
| Fuzzing | Not run in this focused packet | milestone evidence pending | Decoder fuzzing belongs to the planned campaign |

## Residual risk and decision

The selected facade review is complete. The owning AVX implementation file
remains pending its Phase 16 full-file review; native AVX-512 hardware codegen,
sanitizer/fuzz, Apple, and full-matrix evidence remain milestone gates. The
friendly/expert batch subgroup-check separation is fixed under `AUD-161`, with
supported native runtime and performance evidence still pending.
