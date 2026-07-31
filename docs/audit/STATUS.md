# Audit Status

Last updated: 2026-07-31

## Current phase

Phase 12 bootstraps automated analysis before the main file-by-file audit:
CodeQL, Clang Static Analyzer/CodeChecker, clang-tidy, Clippy, ASan/UBSan, and
initial fuzzing. Two upstream rows remain for retained CodeQL and immutable
action maintenance; WASI is not a supported platform. BTx/BTE remains frozen;
CUDA, non-Apple ARM64, and CHERI remain opt-in experimental.

## Ordered gates

| Gate | State | Exit criterion |
| --- | --- | --- |
| 1. Audit infrastructure | complete | Inventory and integrity checks pass |
| 2. Upstream manifest decisions | complete | Every upstream-only commit has a documented approved, rejected, or not-applicable decision |
| 3. Automated analysis baseline | not started | Static analysis, sanitizers, and initial fuzzers execute reproducibly with all results triaged |
| 4. Core scalar audit | in progress | Scalar dependency chain complete with evidence |
| 5. AVX-512 audit | not started | AVX-512 chain and scalar differential coverage complete |
| 6. FFT/middle-product audit | not started | Algorithms, limits, and differential tests complete |
| 7. Supported API/bindings | in progress | C, C++, and Rust surfaces complete |
| 8. Scope cleanup | complete | Unsupported bindings removed; CUDA and non-Apple ARM64/CHERI experimental; Apple Silicon supported; BTx/BTE frozen pending merge |
| 9. Rename and ABI mode | not started | CatBLST rename and generated two-prefix ABI validated |
| 10. Release readiness | not started | CI, packaging, security policy, and release evidence complete |

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
- Experimental: opt-in CUDA; see `CUDA_EXPERIMENTAL.md`.
- Supported scalar platform: Apple Silicon macOS; its native release matrix is
  still required before release.
- Experimental and blocked by default: other ARM64 and Morello/CHERI; see
  `ARM64_CHERI_EXPERIMENTAL.md`.
- Frozen: BTx/BTE pending its required future merge; see `BTX_FREEZE.md`.
- Future project name: CatBLST.
- Future C ABI: build-time selection of exactly one generated prefix,
  `blst_` or `catblst_`; no runtime dispatch.
- Runtime scheduling: C, native C++, and Rust execute entirely on the caller's
  thread; CatBLST creates no threads and owns no scheduler.
- WASI/WebAssembly is not a supported production platform or release gate.

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

16. `cf739f0` and `9c35e86` - allocation-free variable-length Fp decoders and
    their Rust FFI companion were translated and validated in Phase 7A.
17. `f262a6e` - the C2x `bool` ambiguity was translated into the split public
    type header and validated across available C modes in Phase 7B.
18. `b83544b` and `dece82e` - Rust aggregate subtraction and direct dependency
    cleanup were translated and validated in Phase 7C.
19. `07b109d` / `AUD-021` - the upstream bounded-channel optimization was
    superseded by removal of the global Rust pool, unsafe scoped-job lifetime
    erasure, and every internal scheduler path in Phase 8.
20. `2e1948b`, `ed865dc`, and `fe85dc3` - CHERI inverse scratch confinement,
    exact regeneration of 10 ARM artifacts, translator operand hardening, and
    a corrected inverse oracle were validated in Phase 9.
21. `ff66582` and `67fc54b` - pairing-friendly-curves draft-13 generator and
    pairing evidence plus draft-07 basic-signature evidence were translated
    into C++/Rust scalar and AVX tests with accurate release language in Phase
    11.
22. `85808fd` and `33d3f3a` - WASI/Wasmtime CI was rejected as outside the
    supported production matrix; stale Rust WebAssembly build guidance was
    removed in the Phase 12 planning amendment.

## Remaining upstream implementation queue

Two rows remain `reviewing` in one approved implementation packet. The full
constraints are recorded in `UPSTREAM_BULK_APPROVAL.md`.

| Packet | Rows | Required outcome |
| --- | --- | --- |
| Retained CI security and provenance | `68b58bd`, `f62244e` | Add CodeQL for retained C/C++/Rust and modernize retained actions with immutable full-SHA pins; no WASI |

The next phase is the automated-analysis bootstrap in
`AUTOMATED_ANALYSIS_PLAN.md`. Establish and triage the baseline first, then
begin the dependency-ordered hand audit while expanding sanitizer and fuzz
evidence with each subsystem. CodeQL and action provenance close the remaining
upstream rows as part of that bootstrap.

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

- 600 repository files accounted for;
- 154 in-scope production;
- 113 generated production;
- 154 in-scope support;
- 59 benchmark files explicitly excluded;
- 32 frozen BTx/BTE files pending the required merge;
- 52 retained experimental CUDA files;
- 0 planned-removal files;
- 36 documentation/archive files.

Upstream progress:

- 63 of 63 commits have final decisions and the decision gate is complete;
- 35 are closed as already equivalent, deliberate rejection, or verified
  removal-only changes;
- 26 accepted ports/translations are validated;
- 2 are bulk-approved and remain reviewing while their translations and
  validation are completed;
- 0 remain pending.

Findings: 23 total. Twenty-two are validated. One low finding remains open.


## Next action

Do not touch BTx/BTE, and do not include experimental CUDA, non-Apple ARM64, or
CHERI in the primary production-readiness gate. Apple Silicon macOS remains in
the supported scalar release matrix. Implement the automated-analysis baseline
in `AUTOMATED_ANALYSIS_PLAN.md`, including the two retained upstream CI rows.
After baseline triage, start the hand audit with fundamental types, constants,
and byte conversions. Apply `VALIDATION_POLICY.md`: run focused checks in the
change loop, batch related work for subsystem suites, and reserve the full
supported matrix for milestone gates.

No production file should be marked complete before its related implementation
and validation evidence closes the corresponding upstream rows.

Phase evidence is recorded in `PHASE_01.md`, `PHASE_02.md`, `PHASE_03A.md`,
`PHASE_03B.md`, `PHASE_03C.md`, `PHASE_04.md`, `PHASE_05.md`,
`PHASE_06A.md`, `PHASE_06B.md`, `PHASE_06C.md`, `PHASE_07A.md`,
`PHASE_07B.md`, `PHASE_07C.md`, `PHASE_08.md`, `PHASE_09.md`,
`PHASE_10.md`, and `PHASE_11.md`; the remaining dependency order and exit
gates are in `HARDENING_PLAN.md`.
