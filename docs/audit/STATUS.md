# Audit Status

Last updated: 2026-08-02

## Current phase

Phases 12-20 are the main production-audit program defined in
`PRODUCTION_AUDIT_ROADMAP.md`. Phase 12 establishes the shared analysis
platform and initial baseline; automation then remains active while Phases
13-19 perform the subsystem file audits, findings remediation, fuzz/sanitizer
expansion, and release qualification. Phase 20 performs the CatBLST rename and
selectable ABI prefix after the audited surfaces are under control.

The latest Phase 12 hygiene packet, `AUD-038`, removed 171 actionable
commented-code findings and one wholly commented-out unused test header. The
CodeQL rerun leaves exactly 750 generated addition-chain trace comments and no
commented-code result elsewhere; those traces remain as arithmetic provenance.

The clang-tidy baseline is also complete and triaged: 73 unique diagnostics,
72 documented non-defects, and one validated low-severity template-contract
finding (`AUD-039`). The remediation preserves the production width-four hot
loop and native AVX-512 code generation. No clang-tidy fix-it was applied.

The report-only Clippy baseline is complete across default, full-binding,
portable-scalar, Windows portable-AVX Debug, and full-plus-secret-serialization
feature families. Broad linting found one medium-severity unsafe scratch-storage
finding (`AUD-040`), now validated without buffer initialization or portable
Release compilation. The remaining diagnostics are triaged MSRV,
generated-ABI, documentation/style, or deliberate unsafe-boundary review items.
No Clippy fix-it was applied.

The focused scalar sanitizer baseline and `AUD-041` remediation are complete.
The runner now uses pure C/C++ scalar arithmetic with assembly disabled;
ASan/LSan and UBSan each passed 75 tests with 10 expected AVX-only skips and
no sanitizer report. A separate non-ADX generic-assembly profile passed 55/55
focused tests and references no x/ADX symbols. The no-assembly archive
references no assembly primitives. All 18 optimized default-ADX production
objects and the static archive are byte-identical before and after the fix.

Phase 13 hand review has started at common configuration and representation
boundaries. The opening packet records `AUD-042` through `AUD-047`: a high-
severity C ABI exception escape in `blst_miller_loop_n`, empty FFT sizing UB,
unsafe ceiling-division boundaries, an MSVC constant-time evidence gap, global
warning suppression, and partially exposed 32-bit configuration. Remediation
is ordered in `PHASE_13.md`; `config.h.in`, `defines.h`, the affected C API, and
the FFT adapter remain pending until implementation and exact-commit review
evidence are complete.

## Main program phases

| Phase | State | Exit criterion |
| --- | --- | --- |
| 12. Audit platform and initial baseline | in progress | CodeQL, Clang/CodeChecker, clang-tidy, Clippy, ASan/UBSan, initial fuzzers, corpora/artifacts, and triage run reproducibly |
| 13. Fundamental scalar trust chain | in progress | Representations, conversions, portable backend, assembly generators/outputs, and CPU selection complete |
| 14. Scalar and extension fields | not started | Field arithmetic/representation files complete with property, boundary, sanitizer, and oracle evidence |
| 15. Curves and protocols | not started | G1/G2, serialization, hash-to-curve, pairing, signatures, MSM, and uniqueness files complete |
| 16. AVX-512 | not started | AVX files complete with native hardware and scalar differential evidence |
| 17. FFT and middle product | in progress | Scalar/AVX algorithms and APIs complete with reference, sanitizer, and fuzz evidence |
| 18. C/C++/Rust surfaces | not started | Public ABI, native C++, Rust FFI/wrappers, ownership, errors, features, and package inputs complete |
| 19. Support and release evidence | not started | Remaining tests/build/CI/docs/package files complete and full supported release matrix passes |
| 20. CatBLST identity and ABI | not started | Rename plus separate `blst_` and `catblst_` build-time prefix modes validated |

## Scope decisions

- Primary production profiles: scalar and AVX-512.
- Portable-AVX Debug remains available for routine CI and automated analysis.
  Portable-AVX Release is excluded because its optimized compile cost is too
  high for ordinary audit cycles.
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

Phase 12 starts the standing automated-analysis workstream in
`AUTOMATED_ANALYSIS_PLAN.md`; it does not contain the whole audit. CodeQL and
action provenance close the remaining upstream rows during that phase. After
its initial baseline is triaged, Phase 13 begins while automation continues to
expand through Phase 19.

## File-review dependency order

1. Audit platform, build capture, and initial automated baseline.
2. Common types, constants, limbs, bytes/endian, scalar backend, assembly
   generators/outputs, CPU selection, and boundary conversions.
3. Scalar and extension fields.
4. Scalar curve/group operations, serialization, hash-to-curve, pairing,
   signatures, MSM, and uniqueness tree.
5. AVX-512 primitives, RNS/reduction, fields, curves, pairing, MSM, and
   dispatch.
6. FFT and middle product.
7. Public C/native C++ surfaces, Rust FFI and safe wrappers, and package inputs.
8. Tests, vectors, builds, CI, audit controls, documentation, examples,
   packaging, and release evidence.
9. CatBLST rename and selectable ABI prefix.

## Metrics

Current generated inventory:

- 622 repository files accounted for;
- 156 in-scope production;
- 113 generated production;
- 174 in-scope support;
- 59 benchmark files explicitly excluded;
- 32 frozen BTx/BTE files pending the required merge;
- 52 retained experimental CUDA files;
- 0 planned-removal files;
- 36 documentation/archive files.

Main audit work volume:

| Work group | Phases | Total | Complete | Pending |
| --- | --- | ---: | ---: | ---: |
| Scalar foundations, fields, curves, and protocols | 13-15 | 197 | 2 | 195 |
| AVX-512 | 16 | 42 | 0 | 42 |
| FFT and middle product | 17 | 10 | 0 | 10 |
| C/C++/Rust APIs and bindings | 18 | 28 | 0 | 28 |
| Support, tests, builds, CI, audit process, and release controls | 12 and 19 | 166 | 0 | 166 |
| **Total in-scope production/generated/support** | **12-19** | **443** | **2** | **441** |

The rollup is explained in `PRODUCTION_AUDIT_ROADMAP.md`; the per-file
`FILE_MANIFEST.csv` remains authoritative.

Upstream progress:

- 63 of 63 commits have final decisions and the decision gate is complete;
- 35 are closed as already equivalent, deliberate rejection, or verified
  removal-only changes;
- 26 accepted ports/translations are validated;
- 2 are bulk-approved and remain reviewing while their translations and
  validation are completed;
- 0 remain pending.

Findings: 38 total. Thirty-six are validated. Two remain unresolved: one low
finding and one informational finding are open. `AUD-035` is validated with preserved
allocated/caller-table/precomputed modes and 10/10 focused MSVC Release, Linux
ASan, and Linux UBSan tests. `AUD-036` is validated on native
Apple Silicon with AppleClang 17 at commit `0b7ec47`; all 336 C/C++ tests and
the Rust default/full-binding matrix passed. `AUD-037` is validated: pure-C/C++
`cc`, selected-production `vv`, and native intrinsic domains are separated,
and the 9-test direct compatibility target passed on Peach's AVX-512 IFMA CPU.


## Next action

Do not touch BTx/BTE, and do not include experimental CUDA, non-Apple ARM64, or
CHERI in the primary production-readiness gate. Apple Silicon macOS remains in
the supported scalar release matrix. Add the initial bounded public
byte-boundary fuzz platform while continuing the deeper Phase
17 generic FFT/interpolation/middle-product file audit from the boundary
established in `PHASE_17.md`. Coverage remains a late audit/release
metric by user decision; sanitizer, fuzz, and static-analysis work still
accompanies subsystem review. Expand
analyzers, sanitizers, fuzzers, and oracles alongside every subsystem. Apply
`VALIDATION_POLICY.md`: run focused checks in the change loop, batch related
work for subsystem suites, and reserve the full supported matrix for milestone
gates.

No production file should be marked complete before its related implementation
and validation evidence closes the corresponding upstream rows.

Phase evidence is recorded in `PHASE_01.md`, `PHASE_02.md`, `PHASE_03A.md`,
`PHASE_03B.md`, `PHASE_03C.md`, `PHASE_04.md`, `PHASE_05.md`,
`PHASE_06A.md`, `PHASE_06B.md`, `PHASE_06C.md`, `PHASE_07A.md`,
`PHASE_07B.md`, `PHASE_07C.md`, `PHASE_08.md`, `PHASE_09.md`,
`PHASE_10.md`, `PHASE_11.md`, `PHASE_12.md`, and `PHASE_17.md`; the remaining dependency order and exit
gates are in `HARDENING_PLAN.md`.
