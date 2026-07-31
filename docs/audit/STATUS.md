# Audit Status

Last updated: 2026-07-30

## Current phase

Phase 3B Win64 portable-mode and fallback-symbol policy is implemented,
validated, and committed. The existing 100-file assembly tree has a
reproducible, allowlisted source-to-output path, and MSVC now exposes `$1`
fallback symbols only when `BLST_PORTABLE` requires runtime dispatch.
Upstream manifest translation and file-by-file review remain the controlling
workstreams.

## Ordered gates

| Gate | State | Exit criterion |
| --- | --- | --- |
| 1. Audit infrastructure | complete | Inventory and integrity checks pass |
| 2. Upstream manifest | in progress | Every upstream-only commit has a documented decision |
| 3. Core scalar audit | not started | Scalar dependency chain complete with evidence |
| 4. AVX-512 audit | not started | AVX-512 chain and scalar differential coverage complete |
| 5. FFT/middle-product audit | not started | Algorithms, limits, and differential tests complete |
| 6. Supported API/bindings | not started | C, C++, and Rust surfaces complete |
| 7. Planned removals | not started | CUDA, Go, BTx/BTE, and other bindings removed cleanly |
| 8. Rename and ABI mode | not started | CatBLST rename and generated two-prefix ABI validated |
| 9. Release readiness | not started | CI, packaging, security policy, and release evidence complete |

## Scope decisions

- Primary production profiles: scalar and AVX-512.
- Supported language surfaces: C, C++, and Rust.
- C++ support is native; the `bindings/blst.hpp` wrapper over the C ABI is a
  planned removal.
- Pre-production breaking API/ABI changes are acceptable when they remove
  unsafe legacy contracts or materially improve correctness and performance.
- Retained algorithm areas: FFT and middle product.
- Benchmarks remain inventoried but are excluded from file-by-file audit.
- Planned removal: CUDA, Go, BTx/BTE, C#, Java, Node.js,
  Python/SWIG, and Emscripten binding support.
- Future project name: CatBLST.
- Future C ABI: build-time selection of exactly one generated prefix,
  `blst_` or `catblst_`; no runtime dispatch.

## Immediate upstream queue

The first three hardening phases resolved the following escalated
security/correctness and provenance candidates:

1. `1172c18` — Rust `Pairing` DST ownership translated and validated.
2. `96a46e6` — Rust scalar-length validation translated and validated.
3. `6271918` — Rust Pippenger large-system routing translated and validated.
4. `dd5e80b` — stronger local equivalent confirmed and closed.
5. `a3c15d6` — zero-length-salt assembly-path fix validated.
6. `84926f2` — local heap workaround still requires bounded translation.
7. `AUD-009` — deterministic assembly generation is validated.
8. `41b511a` and `14cd09a` — Win64 fallback-symbol guards translated and
   validated through the verified generator path.
9. The remaining accepted assembly translations can now proceed through that
   path in Phase 3C.

## File-review dependency order

1. Build/configuration and public ABI selection.
2. Fundamental types, constants, and byte conversions.
3. Portable backend field arithmetic.
4. Scalar/curve/group operations and serialization.
5. Hash-to-curve, key generation, signature, pairing, and MSM.
6. AVX-512 primitives, RNS/reduction, fields, curve, pairing, and dispatch.
7. FFT and middle product.
8. Public C and C++ headers/wrappers.
9. Rust FFI and safe wrappers.
10. Tests, benchmarks, CI, packaging, and release controls.
11. Reviewed removal of non-production subsystems.

## Metrics

Current generated inventory:

- 617 repository files accounted for;
- 155 in-scope production;
- 114 generated production;
- 130 in-scope support;
- 59 benchmark files explicitly excluded;
- 123 planned-removal files;
- 36 documentation/archive files.

Upstream progress:

- 63 of 63 commits have decisions;
- 32 are closed as already equivalent, deliberate rejection, or verified
  removal-only changes;
- 8 accepted ports/translations are validated;
- 23 remain reviewing while accepted translations, policy decisions, or
  validation are completed;
- 0 remain pending.

Findings: 16 total. Twelve are validated (`AUD-001` through `AUD-004`,
`AUD-007` through `AUD-010`, and `AUD-012` through `AUD-015`). Four remain
open: 2 medium, 1 low, and 1 informational planned-removal item.


## Next action

Begin Phase 3C by applying the remaining accepted upstream x86 generator and
translator changes through the verified output path. Preserve each accepted
instruction schedule exactly, regenerate once, explain every output change,
and rerun scalar plus AVX-512 differential validation.

No production file should be marked complete before its related implementation
and validation evidence closes the corresponding upstream rows.

Phase evidence is recorded in `PHASE_01.md`, `PHASE_02.md`, `PHASE_03A.md`,
and `PHASE_03B.md`; the remaining dependency order and exit gates are in
`HARDENING_PLAN.md`.
