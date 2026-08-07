# File Review: `blst/fields/Fp2.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/fields/Fp2.h` |
| Roadmap phase / queue | Phase 14 selected field facades |
| Reviewed commit | `aa4c7ed187f5a63496f0e94f9284108cbd46c7f1` |
| Reviewer | Codex audit team |
| Date | 2026-08-06 |
| Disposition | in-scope production |
| Subsystem | fields |
| Profiles | selected native C++; scalar and assembly-selected helpers |
| Generation/provenance | handwritten CatBLST facade over upstream-derived scalar implementation |
| Upstream references | not_applicable: CatBLST selected-type facade |
| Related files | `blst/fields/Fp.h`, `blst/scalar/fields/Fp2.h`, `blst/backend/backend.h` |

## Purpose, contracts, and dependency map

This header selects `scalar::Fp2` and the exact scalar/assembly helper names
needed by native C++ callers. It adds no wrapper or alternate layout. The two
assembly helpers are visible only when host assembly is selected and device
compilation is not active. Representation, Montgomery-domain, aliasing,
failure-output, and timing rules are those of the owning scalar declarations.

## Complete checklist evidence

All 47 lines were inspected. Includes are direct; every exported name resolves
to its owner; the host-assembly condition agrees with the backend selector;
there is no storage, conversion, allocation, exception, state, concurrency,
I/O, or arithmetic body. `AUD-196` replaced the implementation-reserved guard
with `BLST_FIELDS_FP2_H`. `AUD-070` covers typed identity use through the
included Fp facade. Checklist sections 1-3, 9, 12, and 13 pass directly;
behavioral sections 4-8, 10, and 11 are N/A here and remain owned by the scalar
and backend reviews. Generated-code additions are N/A.

## Findings

`AUD-070`, `AUD-196` (validated).

## Validation and residual risk

Strict MSVC C++20 parsing and the standalone Release field-boundary executable
pass. The rename changes no token visible after preprocessing and emits no
code. Fp2 arithmetic, assembly-provider, sanitizer, Apple, and native release
evidence remains with the owning implementation and subsystem gates.
