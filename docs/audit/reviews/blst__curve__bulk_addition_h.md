# File Review: `blst/curve/bulk_addition.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/curve/bulk_addition.h` |
| Roadmap phase / queue | Phase 15 scalar curves and protocols, second packet |
| Reviewed commit | `bd45afa58630e7d235d4bbe048e8fc88bc5d5a09` |
| Reviewer | Codex audit team |
| Date | 2026-08-04 |
| Disposition | in-scope production |
| Subsystem | selected native scalar bulk-addition facade |
| Profiles | selected native C++ scalar/default build |
| Generation/provenance | hand-maintained CatBLST selector wrapper |
| Upstream references | not applicable: local native namespace facade |
| Related files | `blst/scalar/curve/bulk_addition.h`, selected curve headers |

## Purpose and trust boundary

This wrapper exposes the selected scalar `points_add` operation in namespace
`blst`. It owns no representation, state, allocation, or arithmetic.

## Review result

The prior wrapper also re-exported `point_head`, `point_tail`, and
`points_accumulate`, making scratch-tree implementation details look like
supported native API. `AUD-112` removes those aliases. The wrapper now has one
include, one namespace, and one selected public alias. Its contract is owned by
the scalar declaration and reviewed C++ documentation.

## Checklist summary

- API: only `points_add` is selected; no C ABI is defined here.
- Memory, lifetime, arithmetic, timing, allocation, and failure behavior:
  delegated unchanged to the reviewed scalar implementation.
- Concurrency: no state.
- Portability: ordinary C++ namespace alias; scalar Debug/Release compile.
- Performance: compile-time alias only, with no wrapper call or conversion.

## Findings and validation

`AUD-112` is validated. Structural search finds no selected alias for the three
detail helpers. The focused scalar Debug/Release target compiles through this
header's containing graph and passes 11/11.

## Residual risk and decision

The wrapper is complete at the recorded commit. Future changes to bulk
semantics belong to the owning scalar implementation review.
