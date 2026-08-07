# File Review: `blst/fields/Fp12.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/fields/Fp12.h` |
| Roadmap phase / queue | Phase 14 selected field facades |
| Reviewed commit | `aa4c7ed187f5a63496f0e94f9284108cbd46c7f1` |
| Reviewer | Codex audit team |
| Date | 2026-08-06 |
| Disposition | in-scope production |
| Subsystem | fields |
| Profiles | selected native C++; scalar implementation |
| Generation/provenance | handwritten CatBLST facade over upstream-derived scalar implementation |
| Upstream references | not_applicable: CatBLST selected-type facade |
| Related files | `blst/fields/Fp2.h`, `blst/scalar/fields/Fp12.h`, scalar GT and pairing code |

## Purpose and contracts

The header selects `Fp4`, `Fp6`, and `Fp12` plus the retained scalar GLV
executors/constants into the main native C++ namespace. It explicitly does not
collapse general Fp12, cyclotomic, and prime-order GT domains; their narrower
preconditions remain on the owning declarations. The aliases create no new
representation, validation state, or runtime dispatch.

## Complete checklist evidence

All 29 lines were reviewed. The guard is non-reserved, includes are direct,
and every name is an exact scalar declaration. There is no pointer operation,
conversion, object lifetime, allocation, exception, mutable state,
concurrency, I/O, or executable arithmetic. Checklist sections 1-3, 9, 12,
and 13 pass directly. Correctness, side-channel, secret-lifecycle, and resource
checks remain on `scalar::Fp12`, GT, and pairing records; the corresponding
facade checks are N/A because no code is emitted. `AUD-070`, `AUD-085`,
`AUD-099`, and `AUD-100` establish the relevant type/domain boundaries.

## Findings

No new finding. Related validated findings: `AUD-070`, `AUD-085`, `AUD-099`,
and `AUD-100`.

## Validation and residual risk

Strict MSVC C++20 parsing and the standalone Release field-boundary executable
pass at the reviewed commit. Arithmetic correctness, timing, sanitizer, Apple,
and native backend evidence remains with the owning scalar/GT/pairing files;
this zero-code facade needs no benchmark.
