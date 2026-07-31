# Audit Status

Last updated: 2026-07-31

## Current phase

Phase 7A translated the upstream arbitrary-length Fp decoder API and its Rust
FFI companion with allocation-free chunking and explicit null-safe empty-input
semantics. BTx/BTE is frozen, not a removal candidate: do not modify, audit,
extract, or remove it until its required future merge is ready. CUDA is
retained as opt-in experimental code outside the scalar/AVX-512
production-readiness claim. Focused C/native C++, Rust FFI, symbol, header, and
AVX-512 build checks pass. Upstream manifest translation and file-by-file
review remain the controlling workstreams.

## Ordered gates

| Gate | State | Exit criterion |
| --- | --- | --- |
| 1. Audit infrastructure | complete | Inventory and integrity checks pass |
| 2. Upstream manifest | in progress | Every upstream-only commit has a documented decision |
| 3. Core scalar audit | in progress | Scalar dependency chain complete with evidence |
| 4. AVX-512 audit | not started | AVX-512 chain and scalar differential coverage complete |
| 5. FFT/middle-product audit | not started | Algorithms, limits, and differential tests complete |
| 6. Supported API/bindings | in progress | C, C++, and Rust surfaces complete |
| 7. Scope cleanup | complete | Unsupported bindings removed; CUDA experimental; BTx/BTE frozen pending merge |
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
- Experimental: opt-in CUDA; see `CUDA_EXPERIMENTAL.md`.
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

16. `cf739f0` and `9c35e86` - allocation-free variable-length Fp decoders and
    their Rust FFI companion were translated and validated in Phase 7A.

## Remaining upstream queue

Thirteen rows remain `reviewing`, grouped into six implementation packets:

| Packet | Rows | Required outcome |
| --- | --- | --- |
| C2x public-header compatibility | `f262a6e` | Resolve the retained `bool`/`_Bool` ambiguity across supported C modes |
| Rust API and dependency hygiene | `b83544b`, `dece82e` | Add `sub_aggregate`; remove the unnecessary direct `rand` development dependency without raising MSRV |
| Rust parallel performance | `07b109d` | Decide and validate bounded grid/add queue translation without regressing the Phase 2 worker-result design |
| CHERI/ARM generator and artifacts | `2e1948b`, `ed865dc`, `fe85dc3` | Port `scbnds` source changes through the verified generator; obtain ARM/CHERI validation or record the hardware block |
| Retained CI/WASI/CodeQL | `85808fd`, `33d3f3a`, `68b58bd`, `f62244e` | Decide Rust WASI support, pin/verify tooling, and modernize retained C/C++/Rust CodeQL/actions |
| Conformance vectors and release claim | `ff66582`, `67fc54b` | Translate authoritative vectors into C++/Rust scalar+AVX tests before accepting the BLS V6 documentation claim |

The next recommended packet is `f262a6e`. It is a narrow public-header
compatibility change and should use standalone C-mode compilation plus the
closest header regression, not a full build matrix.

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

- 588 repository files accounted for;
- 155 in-scope production;
- 113 generated production;
- 141 in-scope support;
- 59 benchmark files explicitly excluded;
- 32 frozen BTx/BTE files pending the required merge;
- 52 retained experimental CUDA files;
- 0 planned-removal files;
- 36 documentation/archive files.

Upstream progress:

- 63 of 63 commits have decisions;
- 32 are closed as already equivalent, deliberate rejection, or verified
  removal-only changes;
- 18 accepted ports/translations are validated;
- 13 remain reviewing while accepted translations, policy decisions, or
  validation are completed;
- 0 remain pending.

Findings: 20 total. Nineteen are validated. One low finding remains open.


## Next action

Do not touch BTx/BTE, and do not include experimental CUDA in the primary
production-readiness gate. Continue with the 13 upstream rows still marked
`reviewing`, starting with `f262a6e`. Apply `VALIDATION_POLICY.md`: run focused
checks in the change loop, batch related work for subsystem suites, and reserve
the full supported matrix for milestone gates.

No production file should be marked complete before its related implementation
and validation evidence closes the corresponding upstream rows.

Phase evidence is recorded in `PHASE_01.md`, `PHASE_02.md`, `PHASE_03A.md`,
`PHASE_03B.md`, `PHASE_03C.md`, `PHASE_04.md`, `PHASE_05.md`,
`PHASE_06A.md`, `PHASE_06B.md`, `PHASE_06C.md`, and `PHASE_07A.md`; the
remaining dependency order and exit gates are in `HARDENING_PLAN.md`.
