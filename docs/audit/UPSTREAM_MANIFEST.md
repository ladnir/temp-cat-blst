# Upstream Commit Manifest

## Comparison boundary

- Upstream: `supranational/blst`
- Common base: `f48500c1fdbefa7c0bf9800bccd65d28236799c1`
- Common-base date: 2025-08-31
- Audited upstream head: `54e6e55674722fc2797ebb4bbb71b26d881eb4b8`
- Upstream-only commits at snapshot: 63
- Local branch at manifest creation: `master`
- Local head at manifest creation: `67a2dfc`

This is a decision ledger, not a cherry-pick list. Each commit must be inspected
in full, compared with the current local implementation, and assigned one final
decision:

- `port`: apply the upstream change substantially as written.
- `translate`: implement the same protection in the local architecture.
- `already_equivalent`: prove the current tree already provides the behavior.
- `not_applicable`: prove the affected surface is absent or being removed.
- `reject`: retain the local behavior with a documented technical reason.
- `defer`: accepted sequencing decision with an owner and exit condition.

Initial relevance is triage only:

- `direct`: touches a retained production or production-support surface.
- `conceptual`: affected surface may be removed, but the defect class or test
  still applies to retained code.
- `removal_candidate`: appears confined to a planned-removal surface; it still
  requires verification before closure.

No row is closed merely because its language binding is scheduled for removal.

## Ledger

| # | Commit | Date | Initial relevance | Subject | Decision | Status | Evidence / notes |
| ---: | --- | --- | --- | --- | --- | --- | --- |
| 1 | `e61acc2` | 2025-09-17 | direct | `bindings/blst.hpp`: smooth some edges | not_applicable | closed | The unsupported C-over-ABI wrapper was removed in Phase 6A; see `AUD-006` |
| 2 | `07b109d` | 2025-09-19 | direct | Rust: use `sync_channel` to optimize allocations | translate | reviewing | Phase 2 replaced the small-work queue with cache-separated result slots and made remaining queues disconnect safely; decide bounded grid/add queue translation during the full Rust performance review |
| 3 | `b83544b` | 2025-09-19 | direct | Rust: add `sub_aggregate` | port | reviewing | Retained additive Rust API; required p1/p2 cneg FFI already exists |
| 4 | `1519e46` | 2025-09-19 | direct | Fix `blst_256_t` right-shift operator typos | port | validated | Both operators use modular right shift with portable/assembly dispatch; independent oracle and zero-shift tests pass; see `AUD-007` and `PHASE_01.md` |
| 5 | `3aa4b71` | 2025-09-19 | removal_candidate | Go README documentation | not_applicable | closed | Go-only documentation removed in Phase 6B |
| 6 | `e7f90de` | 2025-09-19 | direct | Rust crate version bump | reject | closed | Do not claim upstream release identity; define CatBLST crate/version policy during rename |
| 7 | `374b7cd` | 2025-10-02 | removal_candidate | Refresh `.golangci.yml` | not_applicable | closed | Go lint configuration removed in Phase 6B |
| 8 | `17226b2` | 2025-10-10 | conceptual | Go infinity-uncompress test fix | translate | validated | Test intent translated to retained scalar C, native C++, Rust, and AVX public paths with exact 48/96-byte inputs; see `AUD-012` and `PHASE_01.md` |
| 9 | `097818a` | 2025-10-13 | removal_candidate | Go WASM diagnostic binding | not_applicable | closed | Go/WASM surface removed in Phase 6B |
| 10 | `cf739f0` | 2025-10-27 | direct | Add `blst_fp_from_{l,b}e_bytes` | translate | reviewing | API absent; translate into split C header/C++ export architecture with zero-length hardening |
| 11 | `dd5e80b` | 2025-10-27 | direct | Foolproof zero-length `blst_fr_from_{l,b}e_bytes` | already_equivalent | closed | Stronger local fix `d6a55c08`; zeroes output and returns false; C/header regressions passed in `x64-Debug` |
| 12 | `9c35e86` | 2025-10-27 | direct | Refresh generated build outputs | translate | reviewing | Generated companion to `cf739f0`; never apply standalone; regenerate Rust FFI/exports after tooling repair |
| 13 | `b79e2bc` | 2025-10-27 | removal_candidate | SWIG Python 4.4 compatibility | not_applicable | closed | Python/SWIG fix plus Java whitespace only |
| 14 | `6aa1fe7` | 2025-10-28 | conceptual | C++/SWIG Java 4.4 compatibility | not_applicable | closed | Diff is Java/SWIG compatibility only; retained C++ API unchanged |
| 15 | `d497208` | 2025-10-28 | removal_candidate | Harmonize Java Android build | not_applicable | closed | Java Android build only and later reverted |
| 16 | `d7ab1f0` | 2025-10-30 | conceptual | Rework C++/SWIG `bytes_t` | not_applicable | closed | Affects the unsupported C-over-ABI wrapper and SWIG surfaces; native C++ API design is separate |
| 17 | `0438d74` | 2025-10-30 | removal_candidate | Revert Java build modification | not_applicable | closed | Java-only revert; close with `d497208` |
| 18 | `d4ca0f2` | 2025-11-17 | direct | Correct C++ standard support claim | not_applicable | closed | Applied only to the C-over-ABI wrapper removed in Phase 6A |
| 19 | `85808fd` | 2025-11-17 | conceptual | Cache Wasmtime in CI | translate | reviewing | Retain Rust WASI testing only if supported; pin and verify Wasmtime rather than curl-to-shell |
| 20 | `33d3f3a` | 2025-11-17 | conceptual | Fix CI YAML syntax | translate | reviewing | Superseded detail of the retained Rust WASI workflow translation |
| 21 | `7f04257` | 2025-11-25 | removal_candidate | Add Zig bindings | not_applicable | closed | Zig binding/build files only |
| 22 | `68c1d90` | 2025-11-25 | removal_candidate | Exercise Zig tests in CI | not_applicable | closed | Zig CI only |
| 23 | `67fc54b` | 2025-11-25 | conceptual | Update README | defer | reviewing | Reject Zig listing; accept BLS V6 claim only after retained conformance evidence |
| 24 | `8065152` | 2025-11-26 | removal_candidate | Rearrange Zig generator | not_applicable | closed | Zig generator only |
| 25 | `ce49aad` | 2026-02-16 | removal_candidate | Resolve Go lint nit | not_applicable | closed | Go lint-only code removed in Phase 6B |
| 26 | `2594ed1` | 2026-02-16 | removal_candidate | Resolve Python 3.14 warning | not_applicable | closed | Python/SWIG compatibility only |
| 27 | `2e1948b` | 2026-02-16 | direct | Recognize CHERI `scbnds` in ARM translator | port | reviewing | CHERI source/artifacts remain tracked; Phase 3A proves the current translator/output baseline, so this source port can proceed through the verified generator path |
| 28 | `ed865dc` | 2026-02-16 | direct | Confine ARM inverse temporary storage with `scbnds` | port | reviewing | Bound CHERI inverse scratch capabilities and regenerate through the Phase 3A allowlist; ARM/CHERI execution remains required |
| 29 | `fe85dc3` | 2026-02-16 | direct | Refresh generated build outputs | translate | reviewing | Generated companion to accepted ARM/API sources; consolidate through the validated Phase 3A generator rather than applying artifacts directly |
| 30 | `84926f2` | 2026-02-16 | direct | Cap MSM precompute stack utilization | translate | validated | Translated bounded horizontal-row chunking without a heap fallback, retained checked arithmetic and accumulator padding, and validated G1/G2 threshold/layout/infinity behavior with allocation hooks and ASan; see `AUD-005` and `PHASE_04.md` |
| 31 | `68b58bd` | 2026-02-16 | direct | Resolve CodeQL workflow deprecation | translate | reviewing | Move to current CodeQL and retained C/C+++Rust scope; remove C# and private-repo skip |
| 32 | `dafa98f` | 2026-02-16 | removal_candidate | Fix Emscripten CI failure | not_applicable | closed | Emscripten workflow and Zig log text only |
| 33 | `f262a6e` | 2026-03-16 | direct | Avoid C2x ambiguity in `blst.h` | translate | reviewing | Split `blst_types.h` retains direct `bool -> _Bool` macro for pre-final C2x modes |
| 34 | `a3c15d6` | 2026-04-01 | direct | Accept zero-length salt in `blst_keygen_v5` | port | validated | Zero-length copy is skipped while NULL+0 rejection is preserved; exact-vector C/Rust tests pass in portable and assembly/ADX configurations; see `AUD-004` and `PHASE_01.md` |
| 35 | `aee68f7` | 2026-04-01 | conceptual | Accept empty salt in Go `KeyGenV5` | not_applicable | closed | Go wrapper removed in Phase 6B; retained C/Rust regression is linked to `a3c15d6`/`AUD-004` |
| 36 | `06a5711` | 2026-04-01 | removal_candidate | Resolve Go `gosec` lint warning | translate | validated | Applied the bounded-walk intent to the retained shared C/Rust tree with a stronger capacity-aware API, checked layout arithmetic, explicit results, and 32/64-bit depth analysis; see `AUD-011` and `PHASE_05.md` |
| 37 | `7c17bf5` | 2026-04-01 | removal_candidate | Pin Go lint version | not_applicable | closed | Go lint workflow removed in Phase 6B |
| 38 | `c1bb460` | 2026-04-01 | removal_candidate | Resolve TypeScript compiler failure | not_applicable | closed | TypeScript/Node CI only |
| 39 | `f62244e` | 2026-04-01 | conceptual | Resolve workflow warnings | translate | reviewing | Go portion removed in Phase 6B; apply retained action/CodeQL maintenance during CI hardening and drop Node with its binding |
| 40 | `4eb6250` | 2026-04-16 | removal_candidate | Update Zig build for 0.16 | not_applicable | closed | Zig build only |
| 41 | `f41f715` | 2026-04-16 | removal_candidate | Widen Zig CI test matrix | not_applicable | closed | Zig matrix only; carry verified-toolchain lesson into general CI hardening |
| 42 | `dece82e` | 2026-04-17 | direct | Resolve Rust dependabot warning | translate | reviewing | Remove unnecessary direct `rand` dev dependency; use `rand_core` traits and retain MSRV |
| 43 | `6271918` | 2026-05-17 | direct | Fix Rust Pippenger divide-by-zero on large systems | translate | validated | Oversubscribed routing and stronger planner width guards pass deterministic 1,024-CPU G1/G2 and large-bit property tests; see `AUD-003` and `PHASE_02.md` |
| 44 | `db3defd` | 2026-05-17 | conceptual | Fix Go Pippenger divide-by-zero on large systems | not_applicable | closed | Go code removed in Phase 6B; independently corroborates retained Rust `6271918`/`AUD-003` regression |
| 45 | `6056ba3` | 2026-06-04 | removal_candidate | Add experimental Node-API support | not_applicable | closed | Node-API/SWIG only |
| 46 | `9a5af53` | 2026-06-04 | removal_candidate | Exercise Node-API in CI | not_applicable | closed | Node-API CI only |
| 47 | `d29fef4` | 2026-06-04 | removal_candidate | Refine Node addon detection | not_applicable | closed | Node addon detection only |
| 48 | `525a87a` | 2026-07-21 | direct | Omit `sqr_mont_382x` | translate | validated | Removed the unused specialized generator routines, fork declarations/wrapper/test, and exact generated symbols while retaining production `sqr_382x`; see `PHASE_03C.md` |
| 49 | `4a74ffd` | 2026-07-21 | direct | Refine x86-64 Montgomery carry chains | port | validated | Upstream ADX schedules were preserved exactly in generator source and generated artifacts; scalar direct-ADX, portable-dispatch, and AVX-512 suites passed; see `PHASE_03C.md` |
| 50 | `e789f8a` | 2026-07-21 | direct | Refresh generated build outputs | translate | validated | Regenerated through the 100-file deterministic allowlist; all 20 changed outputs match upstream `14cd09a` after newline normalization; see `PHASE_03C.md` |
| 51 | `1172c18` | 2026-07-21 | direct | Bind Rust `Pairing` DST lifetime | translate | validated | Rust now owns DST in a relocatable trailing allocation, including safe reuse/growth and thread moves; see `AUD-001` and `PHASE_02.md` |
| 52 | `96a46e6` | 2026-07-21 | direct | Harden Rust `Pairing.mul_n_aggregate` | translate | validated | Overflow-safe scalar sizing rejects undersized buffers before FFI for both group orientations; see `AUD-002` and `PHASE_02.md` |
| 53 | `23e330f` | 2026-07-21 | conceptual | Remove Zig assumptions about DST lifetime | not_applicable | closed | Zig code removed; owned trailing-DST design informs retained Rust `1172c18`/`AUD-001` |
| 54 | `23cdd8a` | 2026-07-21 | conceptual | Harden Go weighted aggregation interface | not_applicable | closed | Go code removed in Phase 6B; retained Rust fix is `96a46e6`/`AUD-002`; document raw C/C++ contract |
| 55 | `e3d708b` | 2026-07-23 | direct | Update x86-64 translator from Cryptogams | port | validated | Translator now preserves braced directive tokens and contains the corrected CFI spelling; the local file matches upstream `41b511a` after newline normalization; see `PHASE_03C.md` |
| 56 | `33d7f80` | 2026-07-23 | direct | Clarify `mulx_mont_256` limitations | port | validated | Imported the exact operand-bound comment and audited native C++ call sites; new wide-input conversion coverage passes in direct and portable modes; see `PHASE_03C.md` |
| 57 | `bf008ec` | 2026-07-23 | direct | Align exports/reciprocal code with x86 assembly | translate | validated | Phase 1 fixed scalar `Fr::to_mont`; Phase 3C fixed the remaining native C++ conversion and reciprocal operand order and added boundary regressions; see `AUD-008`, `AUD-017`, and `PHASE_03C.md` |
| 58 | `41b511a` | 2026-07-23 | direct | Minimize x86 assembly symbol pollution on MSVC | translate | validated | Translated with one `__BLST_PORTABLE__` MASM guard policy; exact portable/non-portable object allowlists and focused scalar tests passed; see `AUD-010` and `PHASE_03B.md` |
| 59 | `14cd09a` | 2026-07-23 | direct | Refresh generated Windows build outputs | translate | validated | Six Win64 outputs regenerated through the deterministic 100-file allowlist; only the expected guard and unified-spelling changes occurred; see `PHASE_03B.md` |
| 60 | `3ede151` | 2026-07-24 | removal_candidate | Fix Go README typos | not_applicable | closed | Go documentation removed in Phase 6B |
| 61 | `ff66582` | 2026-07-24 | conceptual | Add pairing-friendly-curves test vectors via Python runner | translate | reviewing | Retain authoritative vectors in C++/Rust scalar+AVX tests; drop Python harness |
| 62 | `8aaaef0` | 2026-07-24 | removal_candidate | Bump Zig version | not_applicable | closed | Zig package version only |
| 63 | `54e6e55` | 2026-07-24 | direct | Bump Rust crate version | reject | closed | Do not claim upstream release identity; replace with CatBLST crate/version policy |

Detailed evidence for the first parallel review wave is recorded in
[`upstream/WAVE_01.md`](upstream/WAVE_01.md).
Evidence for the remaining commits is recorded in
[`upstream/WAVE_02.md`](upstream/WAVE_02.md).

## Per-commit review record

For each row, append concise evidence in the table and use this structure below
when the decision needs more explanation:

### `<commit>` — `<subject>`

- Full commit:
- Files changed:
- Defect/feature intent:
- Local equivalent and exact comparison:
- Applicability to scalar:
- Applicability to AVX-512:
- Applicability to C/C++/Rust:
- Decision and rationale:
- Local implementation commit, if any:
- Tests/validation:
- Performance impact:
- Residual risk:

## Completion gate

The upstream phase is complete only when all 63 rows have:

- a final decision other than `pending`;
- status `closed`;
- concrete comparison evidence;
- a linked implementation commit when the decision is `port` or `translate`;
- recorded validation; and
- an explicit residual-risk statement.
