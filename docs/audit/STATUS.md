# Audit Status

Last updated: 2026-07-30

## Current phase

Phase 5 is implemented, validated, and committed. The shared uniqueness tree
now has checked size/alignment, stored capacity, bounded traversal, safe length
comparison, explicit results, and a migrated Rust helper. Scalar, portable
dispatch, targeted AVX-512, Linux ASan, Rust, public-header, and 32-bit compile
validation pass. `blst/common/rb_tree.cpp` and its internal header are the first
production files to complete the full per-file checklist. Upstream manifest
translation and file-by-file review remain the controlling workstreams.

## Ordered gates

| Gate | State | Exit criterion |
| --- | --- | --- |
| 1. Audit infrastructure | complete | Inventory and integrity checks pass |
| 2. Upstream manifest | in progress | Every upstream-only commit has a documented decision |
| 3. Core scalar audit | in progress | Scalar dependency chain complete with evidence |
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

- 623 repository files accounted for;
- 156 in-scope production;
- 114 generated production;
- 135 in-scope support;
- 59 benchmark files explicitly excluded;
- 123 planned-removal files;
- 36 documentation/archive files.

Upstream progress:

- 63 of 63 commits have decisions;
- 32 are closed as already equivalent, deliberate rejection, or verified
  removal-only changes;
- 16 accepted ports/translations are validated;
- 15 remain reviewing while accepted translations, policy decisions, or
  validation are completed;
- 0 remain pending.

Findings: 18 total. Sixteen are validated. Two remain open: one low and one
informational planned-removal item.


## Next action

Begin Phase 6 with the smallest supported-surface removal batch: delete the
unsupported C-over-ABI `bindings/blst.hpp` wrapper and all package,
documentation, and build references (`AUD-006`), then remove Go in a separate
reviewed change. Keep CUDA and BTx/BTE removal isolated so retained FFT and
middle-product dependencies are proved before deletion. Regenerate the file
manifest and run the C, native C++, and Rust matrices after each batch.

No production file should be marked complete before its related implementation
and validation evidence closes the corresponding upstream rows.

Phase evidence is recorded in `PHASE_01.md`, `PHASE_02.md`, `PHASE_03A.md`,
`PHASE_03B.md`, `PHASE_03C.md`, `PHASE_04.md`, and `PHASE_05.md`; the remaining
dependency order and exit gates are in `HARDENING_PLAN.md`.
