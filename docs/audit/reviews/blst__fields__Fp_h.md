# File Review: `blst/fields/Fp.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/fields/Fp.h` |
| Roadmap phase / queue | Phase 14 selected field facades |
| Reviewed commit | `aa4c7ed187f5a63496f0e94f9284108cbd46c7f1` |
| Reviewer | Codex audit team |
| Date | 2026-08-06 |
| Disposition | in-scope production |
| Subsystem | fields |
| Profiles | selected native C++; scalar implementation |
| Generation/provenance | handwritten CatBLST facade over upstream-derived scalar implementation |
| Upstream references | not_applicable: CatBLST selected-type facade |
| Related files | `blst/scalar/fields/Fp.h`, `blst/backend/backend.h`, `test/audit_field_boundaries.cpp` |

## Purpose and trust boundary

This 33-line header selects the scalar base-field type into the main native
C++ namespace and re-exports the constants and inversion helpers required by
retained callers. It owns no representation, arithmetic, storage, or runtime
dispatch. All secret, canonicality, aliasing, and timing contracts remain on
`scalar::Fp` and the backend operations it invokes.

## Complete checklist evidence

The include guard is project-prefixed; includes are direct and acyclic; every
declaration is an exact `using` declaration; no object, pointer conversion,
allocation, exception, mutable state, concurrency, I/O, or data-dependent
operation is introduced. `Fp` is deliberately scalar-selected in the current
public type policy. `AUD-070` verified that its identity consumers now address
live typed subobjects rather than punning C storage. The visible inverse helper
retains its owning scalar contract and is not wrapped by a runtime abstraction.

Sections 1-3, 9, 12, and 13 of `CHECKLIST.md` are satisfied by direct source,
dependency, lineage, API, portability, and zero-overhead inspection. Sections
4-8, 10, and 11 are delegated to the owning scalar/backend records because
this facade executes no behavior. Specialized generated/assembly checks are
N/A.

## Findings

`AUD-070` (validated).

## Validation and residual risk

Strict MSVC C++20 parsing, the scalar Release library, and the standalone
field-boundary executable pass at the reviewed commit. The facade itself needs
no benchmark or sanitizer target because it emits no code. Arithmetic and
constant-time assurance remains with `blst/scalar/fields/Fp.h` and the selected
backend; Apple and broader release evidence remain phase-wide gates.
