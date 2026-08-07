# File Review: `blst/fields/Fr.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/fields/Fr.h` |
| Roadmap phase / queue | Phase 14 selected field facades |
| Reviewed commit | `aa4c7ed187f5a63496f0e94f9284108cbd46c7f1` |
| Reviewer | Codex audit team |
| Date | 2026-08-06 |
| Disposition | in-scope production |
| Subsystem | fields |
| Profiles | selected native C++; scalar implementation |
| Generation/provenance | locally composed CatBLST facade with upstream lineage |
| Upstream references | not_applicable: CatBLST selected-type facade |
| Related files | `blst/scalar/fields/Fr.h`, `bindings/blst_types.h`, `blst/exports.cpp` |

## Purpose and contracts

This header selects the Montgomery scalar-field type and retained constants,
inverse, and secret-key helpers into the main native C++ namespace. `Fr` is
not interchangeable with byte-oriented `blst_scalar`: exact conversion is an
explicit boundary and canonical/secret-key predicates retain their documented
ranges. The facade adds no wrapper, copy, allocation, or runtime selection.

## Complete checklist evidence

All 43 lines were inspected. The header has a direct owning include, a
non-reserved guard, an Apache-2.0/SPDX header with Category Labs attribution,
and exact `using` declarations only. No executable arithmetic, pointer cast,
object storage, exception, state, concurrency, or I/O is introduced.
Checklist sections 1-3, 9, 12, and 13 pass directly; sections 4-8, 10, and 11
remain on the scalar implementation and C conversion boundary. `AUD-070` now
placement-starts the destination type for exact scalar/Fr aliasing.

## Findings

`AUD-070` (validated). Missing provenance was remediated during review without
a separate finding.

## Validation and residual risk

Strict MSVC C++20 parsing, the rebuilt scalar Release library, and the
standalone field-boundary executable pass. The exact-alias conversion erases
its temporary and adds no general field-operation cost. Scalar arithmetic,
constant-time, sanitizer, Apple, and cross-compiler release evidence remains
with the owning implementation and phase gates; no benchmark was needed.
