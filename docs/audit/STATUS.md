# Audit Status

Last updated: 2026-07-31

## Current phase

Phase 6C removed SWIG, C#, Emscripten JavaScript, Java, Node.js, and Python
bindings while retaining shared C/Rust inputs and vectors. BTx/BTE is frozen,
not a removal candidate: do not modify, audit, extract, or remove it until its
required future merge is ready. Scalar C/native C++, AVX-512 build, Rust test,
package-content, and reference checks pass. Upstream manifest translation and
file-by-file review remain the controlling workstreams.

## Ordered gates

| Gate | State | Exit criterion |
| --- | --- | --- |
| 1. Audit infrastructure | complete | Inventory and integrity checks pass |
| 2. Upstream manifest | in progress | Every upstream-only commit has a documented decision |
| 3. Core scalar audit | in progress | Scalar dependency chain complete with evidence |
| 4. AVX-512 audit | not started | AVX-512 chain and scalar differential coverage complete |
| 5. FFT/middle-product audit | not started | Algorithms, limits, and differential tests complete |
| 6. Supported API/bindings | in progress | C, C++, and Rust surfaces complete |
| 7. Scope cleanup | in progress | CUDA removed; unsupported bindings removed; BTx/BTE remains frozen pending merge |
| 8. Rename and ABI mode | not started | CatBLST rename and generated two-prefix ABI validated |
| 9. Release readiness | not started | CI, packaging, security policy, and release evidence complete |

## Scope decisions

- Primary production profiles: scalar and AVX-512.
- Supported language surfaces: C, C++, and Rust.
- C++ support is native; the legacy `bindings/blst.hpp` wrapper over the C ABI
  was removed in Phase 6A.
- Pre-production breaking API/ABI changes are acceptable when they remove
  unsafe legacy contracts or materially improve correctness and performance.
- Retained algorithm areas: FFT and middle product.
- Benchmarks remain inventoried but are excluded from file-by-file audit.
- Go support was removed in Phase 6B.
- SWIG, C#, Emscripten JavaScript, Java, Node.js, and Python bindings were
  removed in Phase 6C.
- Planned removal: CUDA.
- Frozen: BTx/BTE pending its required future merge; see `BTX_FREEZE.md`.
- Future project name: CatBLST.
- Future C ABI: build-time selection of exactly one generated prefix,
  `blst_` or `catblst_`; no runtime dispatch.

## Immediate upstream queue

Completed hardening phases resolved the following escalated
security/correctness and provenance candidates:

1. `1172c18` — Rust `Pairing` DST ownership translated and validated.
2. `96a46e6` — Rust scalar-length validation translated and validated.
3. `6271918` — Rust Pippenger large-system routing translated and validated.
4. `dd5e80b` — stronger local equivalent confirmed and closed.
5. `a3c15d6` — zero-length-salt assembly-path fix validated.
6. `84926f2` — allocation-free bounded scalar MSM precompute translated and
   validated.
7. `AUD-009` — deterministic assembly generation is validated.
8. `41b511a` and `14cd09a` — Win64 fallback-symbol guards translated and
   validated through the verified generator path.
9. `525a87a`, `4a74ffd`, `e789f8a`, `e3d708b`, `33d7f80`, and `bf008ec` —
   accepted assembly removal, schedule, translator, contract, and operand
   changes translated and validated in Phase 3C.
10. `AUD-017` — fork-specific native C++ assembly mapping restored and
    validated across direct and portable dispatch.
11. `AUD-018` — portable ELF/COFF/Mach-O wrappers include both q dispatchers
    and x providers; Linux GCC/ASan validation passed.
12. `06a5711` / `AUD-011` — the retained C/Rust uniqueness tree now bounds
    traversal and storage, reports exact outcomes, and passes Phase 5
    cross-language validation.
13. `e61acc2`, `d4ca0f2`, and `AUD-006` — the unsupported C-over-ABI C++
    wrapper was removed with clean C, native C++, and Rust validation.
14. Go-only upstream rows and `AUD-019` — Go support and its orphaned stale
    amalgamation were removed with clean retained-language validation.
15. SWIG/language-binding upstream rows and `AUD-020` — all remaining
    unsupported language bindings were removed with clean retained-language
    validation.

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

- 585 repository files accounted for;
- 155 in-scope production;
- 113 generated production;
- 138 in-scope support;
- 59 benchmark files explicitly excluded;
- 32 frozen BTx/BTE files pending the required merge;
- 52 planned-removal files, all CUDA;
- 36 documentation/archive files.

Upstream progress:

- 63 of 63 commits have decisions;
- 32 are closed as already equivalent, deliberate rejection, or verified
  removal-only changes;
- 16 accepted ports/translations are validated;
- 15 remain reviewing while accepted translations, policy decisions, or
  validation are completed;
- 0 remain pending.

Findings: 20 total. Nineteen are validated. One low finding remains open.


## Next action

Do not touch BTx/BTE. Continue Phase 6 with CUDA removal as an isolated reviewed
change, or return to the remaining upstream/core audit if CUDA is deferred.
Regenerate the file manifest and run the C, native C++, and Rust matrices after
each implementation batch.

No production file should be marked complete before its related implementation
and validation evidence closes the corresponding upstream rows.

Phase evidence is recorded in `PHASE_01.md`, `PHASE_02.md`, `PHASE_03A.md`,
`PHASE_03B.md`, `PHASE_03C.md`, `PHASE_04.md`, `PHASE_05.md`,
`PHASE_06A.md`, `PHASE_06B.md`, and `PHASE_06C.md`; the remaining dependency
order and exit gates are in `HARDENING_PLAN.md`.
