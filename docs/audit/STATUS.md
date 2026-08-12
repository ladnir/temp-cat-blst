# Audit Status

Last updated: 2026-08-11

Waves 1 and 2 of the parallel file audit are integrated through implementation
commit `535c4c6`. Wave 1 closes the scalar C-boundary lifetime and public
type/contract packets. Wave 2 records `AUD-198` through `AUD-207`, repairs the
approved scalar/extension-field findings, and extends `AUD-060` so portable
conditional negation matches assembly for every nonzero flag. Focused scalar
Release evidence passes 7/7 and optimized small-shift/cneg inspection passes.
The carryover pass validates `AUD-199` with retained Debug REDC bound tests and
validates `AUD-201` with MSVC/Clang/GCC inverse code-generation evidence. No
benchmark, CUDA, BTx, sanitizer, Apple job, or broad runtime sweep ran.

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

The first Phase 13 hand-review packet is complete. `config.h.in` and
`defines.h` are reviewed at exact commit `7719303ad809dd7a738421c721be8ebe91bf6ce5`.
`AUD-042` through `AUD-049` are validated: the C multi-Miller boundary is
allocation-free, scalar infinity handling is corrected, sizing and FFT-root
edge cases are hardened, the MSVC/Clang/GCC value barriers have optimized
code-generation evidence, broad warning suppression is gone, and configuration
is explicitly 64-bit-only. The larger C API and FFT files remain pending their
own full reviews. The retained-C-ABI inventory opened `AUD-050` because legacy
`void` keygen still allocates an HKDF information buffer; it does not throw but
must become allocation-free in a dedicated keygen packet.

The portable-backend review subsequently upgraded `AUD-050` to high severity
and opened `AUD-051` through `AUD-054`. Remediation is in progress for
allocation-free key derivation, defined zero-count public shifts, exact ADX
output typing, portable division borrow propagation, and MSVC zero-input bit
scanning.

`AUD-050` through `AUD-054` are now validated. Key derivation is allocation-
free, zero-count public shifts are identity, portable division and MSVC bit
scanning are corrected, and ADX outputs have exact intrinsic-compatible types.
The complete Phase 13 file reviews for `backend_portable.h` and
`backend_portable.cpp` are recorded at commit
`5606e84262401225c4e02b3ec786b849c85adfcc`. Broader field, keygen, and public-
surface reviews remain pending in their scheduled phases.

The `backend.h` hand review reopened `AUD-051` after finding that the initial
zero-shift remediation missed supported native C++ template operators. It also
opened `AUD-055`: production host SHA-256 compression is incorrectly pinned to
the portable backend. Both findings are now validated: zero identity is
centralized at every selected shift boundary, production SHA compression again
selects assembly, and focused behavior plus code-generation evidence passes.
The exact `backend.h` review is recorded at
`dc07be802c39a3c04a17347a5cec64f74d95ecd7`; no broad build or benchmark was
required.

The common-foundation remediation is locally complete. `AUD-059`, `AUD-060`,
`AUD-061`, `AUD-064`, and `AUD-065` are validated; `AUD-058` and `AUD-062` are
fixed pending Apple CI, and `AUD-063` is fixed pending the next retained Rust
compile. Exact per-file records now cover `cpuid.cpp`, `vect.h`, `bytes.h`,
`errors.h`, and `timer.h`; the completed `defines.h` record was refreshed for
the canonical `uintptr_t` representation. Three iterator/mdspan helpers have
no supported production consumer and remain in experimental CUDA/parallel-MSM
scope. Final local evidence used one MSVC worker and no benchmark, CUDA job,
full suite, sanitizer sweep, or parallel build.

The existing Apple workflow run `30826190888` was inspected without rerunning
or pushing. Its arm64 architecture check and configuration passed, but both
the native C/C++ and Rust Release jobs failed because
`blst/scalar/curve/point.h` used `printf` through a Debug-only transitive
include. `AUD-066` records the regression. Commit `79010a1` adds the owning
header's direct dependency, and focused strict Release-style `GT.cpp` object
compiles pass with MSVC 19.50, Clang 18.1.3, and GCC 13.3. The finding remains
fixed rather than validated until a later, intentionally scheduled Apple run;
the old run does not validate the newer CPU-capability or endian fixes.

The first Phase 14 Fr/Fp wave is implemented. `AUD-067`, `AUD-068`,
`AUD-072`, `AUD-073`, and `AUD-074` are validated; `AUD-071` is fixed pending
retained Rust execution. Full-width x86 repeated-square generation is committed
at `7c0e399`, and the remaining straightforward field boundary packet is at
`b6a51bb`. Commit `61c8cd4` validates the `AUD-070` object-lifetime packet with
an explicit little-endian platform requirement, typed Fp12 identity storage,
raw-limb reciprocal scratch, and documented scalar/Fr representation and alias
contracts. Commit `a720184` retires the isolated generic `blst_t.hpp` shim and
its three dedicated tests as `AUD-078`; supported native C++ remains on the
canonical field, curve, and protocol types. Complete host reviews now cover
`Fr.h`, `Fp.h`, `Fr.cpp`, and `Fp.cpp`; the headers remain in review for their
listed cross-platform, sanitizer, and independent-oracle evidence. Focused
evidence used one build worker, three
selected public tests, standalone C/C++ probes, strict GCC/Clang syntax checks,
100-output assembly verification, and MSVC optimized code inspection. No full
test sweep, benchmark, sanitizer, CUDA, CI push, or Apple run was performed.

The next Phase 14 wave completed full hand review of `Fp2.h`, `Fp12.h`, and
`Fp12.cpp` at `f0a9d9e`. Commit `0ee24de` fixes `AUD-079` through `AUD-082`,
`AUD-086`, `AUD-087`, and the cross-cutting `AUD-088`/`AUD-089`; it accepts and
documents the trusted-schedule and Fp2 failed-output contracts in
`AUD-083`/`AUD-084`. The later Phase 15 follow-up closes `AUD-085` with a
canonical, membership-checking default GT decoder and a canonical
`deserialize_unchecked` path that defers only membership.

The primary 256-bit group packet restores ordinary subgroup multiplication and
the two inherited generic `unchecked_mult` exceptions, prepares scalars without
allocation, and adds genuinely variable-time fast paths across C, native C++,
and Rust. Canonical G1/G2/GT scalars retain the fixed-public-bit-length
constant-time contract; noncanonical inputs are accepted modulo r with a
documented reduction-required timing distinction. Scalar G1
and G2 select wNAF-5; native AVX selects wNAF-5 for G1 and wNAF-4 for G2. GT
uses the fixed unsigned radix-z splitter for CT, retains the balanced scalar VT
splitter, and selects the faster unsigned native-AVX VT splitter. FFT and middle
product and weighted pairing explicitly use public-scalar operations, while
ordinary Pippenger and its one-point fallback use the subgroup contract.

The subsequent API-wide group-domain inventory defines the shared contract
vocabulary and classifies 241 public C symbols into 52 disjoint families.
The recalibration validates `AUD-091`, `AUD-092`, `AUD-093`, and `AUD-095`:
ordinary scalar/AVX MSM has one subgroup/modulo-r contract independent of
scalar width and point count; the experimental unchecked MSM surface is gone;
the GT-only GLV4 MSM has only GT-named C/Rust entries; and public membership
predicates compose the curve equation with the subgroup relation. The
allocation-free full-width AVX path uses caller-owned split scratch. A focused
no-BTx portable-AVX Debug target passes 5/5 locally, including the Fp2 affine
membership regression discovered during this packet. The remaining contract
queue covers decode-versus-validate policy, incomplete addition, pairing
validation, signing/SK preconditions, GT construction, and cyclotomic wording.

The supported GCC 15.2 native build also opened low-severity `AUD-102`: GCC
repeats `-Wignored-attributes` for `__m512i` aliases used as `std::array`
template arguments in the selected AVX vector/change-base path. The build and
focused execution pass, but the exact alignment and code generation must be
recorded before choosing a narrow suppression or any type/layout remedy.
`AUD-102` is validated at `d222e05`. `AVXScalarArray` remains an exact
`std::array` alias and confines the GCC-only diagnostic suppression to that
single alias declaration. Compile-time checks bind the 64-byte element size,
aggregate size, and element alignment. GCC 15.2 accepts the construct with
`-Werror=ignored-attributes` and emits a direct aligned ZMM load plus `vpaddq`;
a full local GCC 13.3 MatrixNoK production probe also passes the warning as an
error and retains direct IFMA/ZMM code. No wrapper, copy, runtime code, broad
build, or benchmark was introduced. Exact correspondence rows `AVX-P0-001`
and `AVX-P2-003` are reaccepted at the new source hashes.

The focused group target now passes 10/10 in scalar MSVC Debug/Release and
Peach native AVX-default Release; the dedicated AVX contract target passes 5/5
under portable-AVX Debug and Peach native AVX-512/IFMA Release. The focused
Rust method test and full-Rust-binding AVX type-check also pass. Sequential
CPU-pinned Peach benchmarks used one process at a time and compare the exact
`2744efc` baseline with the final candidate. The largest primary regression is
`+1.193%` (one-point GT MSM), every other primary regression is at most
`+0.458%`, one-point G1/G2 MSM improves by 39%/53%, and GT CT multiplication
improves by about 2.9%. GCC 15.2 inspection finds the approved one-branch
canonical fast path, retained secret erasure, and no out-of-line
selected-vector wrappers. `AUD-088` and `AUD-089` are validated. No full suite,
sanitizer, CUDA, BTx/BTE, portable-AVX Release, CI push, or Apple run was
performed for this packet.

The focused `AUD-090` correctness follow-up now erases the AVX G1/G2 CT
division scratch, split components, Booth schedules, and executor scalars while
leaving public-scalar VT/MSM paths unchanged. Portable-AVX Debug and Peach
native Release group tests pass. GCC 15.2 retains every expected fixed-size
volatile erase region in the optimized call sites; the cross-compiler erase
primitive itself was already validated by `AUD-045`. No benchmark or broad
build was started for this narrow hardening change.

The line-by-line `blst/avx/curve/scalar_mult.h` review is complete with a
permanent independent oracle. The focused eight-test target passes under scalar
Release, portable-AVX Debug, and Peach native AVX-512/IFMA Release; native CT
symbol sizes are unchanged. The file remains `in_review` because caller tracing
opened high-severity `AUD-091`: automatic 256-bit AVX Pippenger had applied a
canonical-scalar, subgroup-only GLV/GLS split under a general MSM contract.
The subgroup-only API recalibration resolves that caller defect.

The cross-language group-domain contract audit now accounts for all 241
relevant public C symbols in 52 disjoint machine-checked families and traces the
native C++ and safe/raw Rust surfaces. It records actual Release validation
separately from Debug assertions and mathematical preconditions. `AUD-092`
through `AUD-100` capture GT/Fp12 naming, complete subgroup predicates,
on-curve-only decoding, the resolved AVX MSM discontinuity, incomplete point-add
preconditions, pairing validation policy, and optimized signing/GT contracts.
`AUD-101` is validated: weighted C and safe Rust pairing reject malformed or
oversized scalar shapes before mutation/native entry, AVX scalar boundaries
reject more than 256 bits, and the shared canonicalizer fails closed in Release
instead of overflowing fixed storage. Focused scalar Release C/C++ and Rust
tests pass. The stale uppercase-coordinate test references were repaired, and
the current native AVX boundary target executes 5/5 on Peach. No sanitizer
sweep, portable-AVX Release build, or full suite ran.

## Main program phases

| Phase | State | Exit criterion |
| --- | --- | --- |
| 12. Audit platform and initial baseline | in progress | CodeQL, Clang/CodeChecker, clang-tidy, Clippy, ASan/UBSan, initial fuzzers, corpora/artifacts, and triage run reproducibly |
| 13. Fundamental scalar trust chain | in progress | Representations, conversions, portable backend, assembly generators/outputs, and CPU selection complete |
| 14. Scalar and extension fields | in progress | Field arithmetic/representation files complete with property, boundary, sanitizer, and oracle evidence |
| 15. Curves and protocols | in progress | G1/G2, serialization, hash-to-curve, pairing, signatures, MSM, and uniqueness files complete |
| 16. AVX-512 | not started | AVX files complete with native hardware and scalar differential evidence |
| 17. FFT and middle product | in progress | Scalar/AVX algorithms and APIs complete with reference, sanitizer, and fuzz evidence |
| 18. C/C++/Rust surfaces | in progress | Public ABI, native C++, Rust FFI/wrappers, ownership, errors, features, and package inputs complete |
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
16. `AUD-070` — inactive-union and whole-object field punning was replaced by
    typed identity storage, raw-limb reciprocal scratch, and an explicit
    little-endian platform contract.
17. `AUD-078` — the isolated generic-modulus C++ shim and its dedicated tests
    were retired; the canonical native C++ implementation remains supported.

18. `cf739f0` and `9c35e86` - allocation-free variable-length Fp decoders and
    their Rust FFI companion were translated and validated in Phase 7A.
19. `f262a6e` - the C2x `bool` ambiguity was translated into the split public
    type header and validated across available C modes in Phase 7B.
20. `b83544b` and `dece82e` - Rust aggregate subtraction and direct dependency
    cleanup were translated and validated in Phase 7C.
21. `07b109d` / `AUD-021` - the upstream bounded-channel optimization was
    superseded by removal of the global Rust pool, unsafe scoped-job lifetime
    erasure, and every internal scheduler path in Phase 8.
22. `2e1948b`, `ed865dc`, and `fe85dc3` - CHERI inverse scratch confinement,
    exact regeneration of 10 ARM artifacts, translator operand hardening, and
    a corrected inverse oracle were validated in Phase 9.
23. `ff66582` and `67fc54b` - pairing-friendly-curves draft-13 generator and
    pairing evidence plus draft-07 basic-signature evidence were translated
    into C++/Rust scalar and AVX tests with accurate release language in Phase
    11.
24. `85808fd` and `33d3f3a` - WASI/Wasmtime CI was rejected as outside the
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

- 681 repository files accounted for;
- 151 in-scope production;
- 113 generated production;
- 234 in-scope support;
- 60 benchmark files explicitly excluded;
- 32 frozen BTx/BTE files pending the required merge;
- 55 retained experimental files, including CUDA;
- 0 planned-removal files;
- 36 documentation/archive files.

Main audit work volume:

| Work group | Phases | Total | Complete | Not complete |
| --- | --- | ---: | ---: | ---: |
| Scalar foundations, fields, curves, and protocols | 13-15 | 192 | 30 | 162 |
| AVX-512 | 16 | 42 | 1 | 41 |
| FFT and middle product | 17 | 10 | 0 | 10 |
| C/C++/Rust APIs and bindings | 18 | 28 | 2 | 26 |
| Support, tests, builds, CI, audit process, and release controls | 12 and 19 | 226 | 17 | 209 |
| **Total in-scope production/generated/support** | **12-19** | **498** | **50** | **448** |

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

Findings: 314 total: 283 validated, 23 fixed pending final evidence, four accepted
risks, one needs a design decision, and three remain open. `AUD-035` is validated with preserved
allocated/caller-table/precomputed modes and 10/10 focused MSVC Release, Linux
ASan, and Linux UBSan tests. `AUD-036` is validated on native
Apple Silicon with AppleClang 17 at commit `0b7ec47`; all 336 C/C++ tests and
the Rust default/full-binding matrix passed. `AUD-037` is validated: pure-C/C++
`cc`, selected-production `vv`, and native intrinsic domains are separated,
and the 9-test direct compatibility target passed on Peach's AVX-512 IFMA CPU.


## Next action

Do not touch BTx/BTE, and do not include experimental CUDA, non-Apple ARM64, or
CHERI in the primary production-readiness gate. Apple Silicon macOS remains in
the supported scalar release matrix. `point.h` and `point.cpp` are complete at
`07d79a48`; `AUD-096` accepts the inherited narrow unsuffixed-add contract with
the complete `add_or_double`/`dadd` alternative documented. Scalar and selected
bulk-addition headers are complete at `bd45afa`; `AUD-109` through `AUD-112`
are validated. Scalar/selected MSM and the safe Rust multipoint wrapper are
complete at `af09f66`; `AUD-113` through `AUD-120` are validated. The
batch-affine MSM review is complete at `913c52f`; `AUD-121` through `AUD-125`
are validated. The `ec_mult.h` hand review and `AUD-126` through `AUD-130`
are complete at `b985529`. Scalar `GT.h`/`GT.cpp` and the selected facade are
complete at `3f72403`: native Miller results are distinct from GT, final
exponentiation is named and typed explicitly, prepared-line factors have their
own type, and the scalar/AVX generator is canonical static data. Scalar and
selected map-to-group headers are complete at `a1ef4b3`: typed hash-to-field
outputs, static host constants, inversion-free AVX bridges, exact RFC vectors,
and public contracts are validated. SHA-256 is complete at `5c69ad6`: the C,
native C++, and raw Rust contracts are accounted for, independent one-shot and
streaming KATs pass, functional zeroing is wide-store optimized, and the public
wrapper securely erases its local context. Key generation is complete at
`f36e35b`: HMAC storage has a valid C++ lifetime, RFC HMAC/HKDF and all retained
version contracts have direct evidence, internal helpers no longer leak, and
complete secret scratch erasure remains in optimized output. Continue Phase 15
with the aggregate packet complete at `091dfb2`: the state machine is
documented, empty export and infinity accumulation are defined, rejected items
are failure-atomic, relocatable DST state has no invented pointer, and safe
Rust uses sealed compile-time G1/G2 orientation instead of `dyn Any`. The
low-level accumulator export is intentionally renamed to the allocation-free
out-parameter `blst_pairing_get_miller_result`; `AUD-152` through `AUD-159` are
validated. `AUD-085` now has a checked, canonical, failure-atomic default GT
decoder plus a canonical `deserialize_unchecked` path that defers only group
membership; focused scalar Debug/Release and selected portable-AVX Debug tests
pass at `0972269`, and `AUD-085` is validated. Coverage
remains a late audit/release metric by user decision; sanitizer, fuzz, and
static-analysis work still accompanies subsystem review. Expand analyzers,
sanitizers, fuzzers, and oracles alongside every subsystem. Apply
`VALIDATION_POLICY.md`: run focused checks in the change loop, batch related
work for subsystem suites, and reserve the full supported matrix for milestone
gates.

The `AUD-094`/`AUD-160` point-domain policy now uses only on-curve and subgroup
preconditions without new validation-state types. Decode and the approved basic
single-point primitive set opt into on-curve use; batch helpers, MSM, raw
Miller/pairing, and protocol algorithms remain subgroup-only regardless of
incidental wider implementation behavior. The scalar/AVX C, selected native
C++, and safe Rust declarations now match the exact cross-language matrix in
`docs/audit/POINT_DOMAIN_MATRIX.md`. Order-11 G1 and order-13 G2 fixtures pass
focused scalar C/native, portable-AVX Debug, and both safe Rust variant flows,
including compressed/uncompressed decoding, Serde, complete validators, and
the approved basic primitive set; `AUD-094` and `AUD-160` are validated.
`AUD-161` preserves the caller-controlled AVX batch predicate as an expert
relation evaluator and adds friendly G1/G2 membership APIs that own 81
domain-separated uniform ternary rows and the subgroup shift from one fresh
unpredictable seed. Portable-AVX Debug low-order and failure-path tests, C99
header compilation, Rust checking, and native MSVC compilation pass. The
finding is fixed pending supported native runtime and performance evidence.
`AUD-162` and `AUD-163` remove the empty pairing facade and unused
scalar fixed-window re-export facade. The selected `blst/curve/point.h` review
is complete at `5d18e42`: `AUD-164` adds true AVX exact-scalar G1/G2 methods,
`AUD-165` makes raw AVX decode failure-atomic, `AUD-166` completes native point
contracts, `AUD-167` records intentional generic/helper visibility, and
`AUD-168` replaces the reserved guard. Portable-AVX Debug contract tests pass
8/8. The owning `blst/avx/curve/point.h` Phase 16 review is now active:
`AUD-169` covers scalar-backed AVX signing, `AUD-170` covers conversion-heavy
point observation, `AUD-171` restores the inherited pre-decorated `to_pk2`
contract, and `AUD-172` records the distinct-object `times_minus_z` invariant.
Commit `5e4ab8d` fixes all four and restores selected `sign2` under `AUD-173`;
both scalar-default and AVX-default portable-AVX Debug contract targets pass
12/12. `AUD-170` through `AUD-173` are validated. `AUD-169` remains fixed
pending native AVX codegen/runtime; sanitizer, fuzz, Apple, and full matrices
remain gates.

For every newly reviewed externally reachable family, update declaration-level
C, selected native C++, and safe Rust contracts in the same packet and add it
to `PUBLIC_API_DOCUMENTATION.csv`. Generated raw Rust FFI continues to inherit
the C contract until its Phase 18 generator review. `AUD-103` is validated:
selected scalar and AVX GT conditional negation are constant-time, while the
scalar variable-time MSM helper retains its public signed-digit branch locally.

The Phase 15 hash-to-field packet is reviewed in the current worktree.
`AUD-174` through `AUD-179` replace limb-punned and partially heap-backed XMD
with byte-safe exact allocation-free expansion, enforce the RFC 9380 8160-byte
limit, erase the safe Rust scalar intermediate, add official XMD and BLS12-381
field vectors, and narrow the selected facade. Focused scalar Debug/Release and
safe Rust tests pass, and optimized MSVC emits vector XOR with fixed stack
storage and no allocator call. The two headers are complete at
`ab212fe80784dc516ca995b6823db6227141c52c`.

The generated scalar Fp addition-chain packet is complete at
`529a79377b6ef848ec90fce35e7580a9f451fc4f`. `AUD-181` through `AUD-184` add exact retained symbolic proof, correct
native zero contracts, restrict the schedules to internal scalar Fp, use typed
assignment, and pin upstream provenance without changing either arithmetic
schedule. Scalar MSVC Debug/Release focused tests pass 3/3. Proposed `AUD-180`
generic scratch erasure was declined; constant-time arithmetic does not claim
comprehensive transient erasure.

Normal GoogleTest audit regressions are now consolidated into `blst_tests`
instead of accumulating one executable per reviewed subsystem. Scalar, AVX,
and sanitizer configurations select their source sets at configure time; the
sanitizer runner uses the same target with its existing narrow BTx/BTE-free
boundary set. Only the special both-enabled native/emulated AVX-512 compatibility
probe remains a separate test executable because it has distinct sources and
compiler flags.

Phase 18 file-by-file review has started with the public compatibility
umbrellas. `AUD-185` records that the FFT header, added after the local binding
split, was omitted from `blst_aux.h` even though the documented umbrella
contract promised every supported split surface. `blst.h` and `blst_aux.h` now
expose scalar, AVX, middle-product, and FFT declarations transitively. Strict
MSVC C11 and C++20 consumers compile, link, and execute an n=1 FFT through only
`blst.h`; no ABI symbol or runtime implementation changed. Both umbrella files
are complete at `98faac3fb001f59f72fd1b131c278c101953f1a0`.

No production file should be marked complete before its related implementation
and validation evidence closes the corresponding upstream rows.

Phase evidence is recorded in `PHASE_01.md`, `PHASE_02.md`, `PHASE_03A.md`,
`PHASE_03B.md`, `PHASE_03C.md`, `PHASE_04.md`, `PHASE_05.md`,
`PHASE_06A.md`, `PHASE_06B.md`, `PHASE_06C.md`, `PHASE_07A.md`,
`PHASE_07B.md`, `PHASE_07C.md`, `PHASE_08.md`, `PHASE_09.md`,
`PHASE_10.md`, `PHASE_11.md`, `PHASE_12.md`, `PHASE_14.md`, `PHASE_16.md`, and
`PHASE_17.md`; the remaining dependency order and exit
gates are in `HARDENING_PLAN.md`.

Phase 16 now includes the isolated MatrixNoK G1 candidate at `7e731d9`.
Native Clang 21 correctness passes 11/11 focused tests, including independent
scalar complete-operation, CT/VT GLV, and FFT round-trip oracles. The candidate
G1 FFT is 9.4% to 9.9% faster than the incumbent and within 1.7% of pinned
VROOM across 256 through 16,384 elements. Production dispatch is unchanged;
formal proof and the remaining G2/Fp12/pairing slices remain open.

The later production integration is now natively exercised at `6ce0d44` on
Peach with Clang 21.1.8, native AVX-512IFMA, and one build worker. The fresh
Release production library and focused executable build; 34/34 MatrixNoK,
C-API-contract, and native-point tests pass, as do all three structural/formal
audit CTests. `AUD-228` fixes redundant-RNS zero observation in public infinity
queries and secret normalization, including G1/G2 identity signing. Identical-
flags code-generation probes have no selected-vector/conversion symbols and
unchanged stack usage; no benchmark ran. `AUD-210`, `AUD-225`, `AUD-226`, and
`AUD-228` are validated, and exact-hash `AVX-P4-002` is accepted. At that
checkpoint, `AUD-227` remained the public curve/GLV/GLS formal gap. Native object inspection opened
`AUD-229` for function-local initialization guards around AVX owners and
public constants.

`AUD-229` is validated at `fe09e26`: selected curve, GLV/GLS, GT, and FFT
owners are constant-initialized; public C literals reference those owners
directly; and the new aligned curve literals have an independent canonical
regenerator and corruption self-test. The focused MSVC portable-Debug target
passes 35/35, both constant validators and the dispatch validator pass, and
the two production C API objects contain no MSVC dynamic-initialization
symbols. Peach Clang 21.1.8 Release also passes 35/35 focused native tests and
all five audit gates. Native guard/atexit relocations fall from 231 to zero,
the selected code-generation probe is byte-identical, hot owner stack frames
are unchanged, and the AVX C API object shrinks by about 19 KiB of text while
eliminating its BSS owners. The five exact-hash proof rows are reaccepted.
`AUD-230` makes all 51 mapped hashes LF-canonical after Peach exposed that
raw-worktree hashing depended on Windows CRLF. No benchmark ran because the
arithmetic probe was byte-identical and hot code generation did not regress.

`AUD-227` is validated at proof commit `f1cfde9`. Accepted `AVX-P4-003` proves
the shared public complete-add/double formulas and curve
preservation, exact G1/G2 scalar decomposition, and the selected constant-time
Booth executors. The variable-time wNAF source recurrence, recomposition, and
exact 129/65-slot exhaustion bounds are also proved. The concrete beta and
minus-psi constants, coordinate maps, curve preservation, and signed-to-
unsigned subgroup-relation bridges are now proved as well. Remaining closure
work is the equivalence between BLST's relation predicates and the
prime-order subgroup definition, followed by final public-dispatch/full
group-law composition. This step changes only Lean and audit infrastructure;
it runs no benchmark or native build.

The subgroup-reduction packet further proves the exact `sigma³ = 1` and
`psi⁶ = -1` coordinate identities, the G1 complete-add endomorphism law, both
BLS cofactor gcd reductions, and conditional G1/G2 equivalences between the
source relation predicates and `[r]P=0`. Full curve cardinality and
prime-subgroup generator span remain explicit theorem premises. The G2
semilinear projective-additivity bridge is now proved by composing a
generic ring-map theorem with a diagonal-scaling theorem; the resulting common
homogeneous scale is proved nonzero from the literal coefficient.
`CurveGroupBridge.lean` now provides a concrete Mathlib projective
`AddCommGroup` carrier for selected G1, proves the nonzero on-curve conversion,
proves CatBLST homogeneous equivalence iff projective equality, and proves zero
and negation correspondence. `Fp2Field.lean` proves the selected quadratic
extension is a field from `p mod 4 = 3`, nonvanishing norm, and an explicit
conjugate-over-norm inverse. The shared projective bridge is now instantiated
for selected G2 too, with exact equality, zero, and negation correspondence.
`CurveGroupAdd.lean` closes the complete-add/group-add equivalence for every
valid selected G1 and G2 homogeneous representative. The proof covers
infinity, doubling, inverse, and distinct-affine cases, proves output
representatives nonzero, and proves the selected fields have no affine
2-torsion by finite-field cube criteria. `CurveBinaryMul.lean` then proves a
logarithmic evaluator built from that exact complete-add owner equal to
Mathlib projective scalar multiplication. `CurveGenerators.lean` validates
the standard BLS12-381 G1/G2 coordinates, proves both generators are killed by
`r`, and checks the selected beta/minus-psi eigenbases. `CurveDispatch.lean`
composes public masking and reduction with the exact selected Booth/wNAF
schedules and proves all four G1/G2 CT/VT dispatches over the conventional
cyclic subgroups. Full ambient-curve cardinality is therefore not required for
the multiplication contract. Proof commit `24de7b4` adds
`CurveMembership.lean`, packages the literal sigma/psi maps as additive maps on
the selected projective groups, and proves
`standard cyclic subgroup ⊆ accepted relation ⊆ r-torsion` for both G1 and G2
without a point-count premise. Commit `1fdc201` then defines the exact accepted
relation kernels as `ZMod r` modules and proves the selected unsigned GLV/GLS
bases plus all four CT/VT multiplication schedules on those domains. The
validator-to-optimized-arithmetic bridge is therefore formal even without the
span theorem. Commit `65b54f6` proves `r` prime by Lucas certificate, proves
both standard generators have exact order `r`, and reduces span equality to
the explicit upper bounds `card accepted ≤ r`. Commit `f6c251b` proves the G1
bound and exact span equality from a generic `2p+1` curve bound and `r² >
2p+1`. Commit `475221a` proves the accepted G2 kernel has cardinality `r` or
`r²`. The corrected BLS12 criterion in El Housni, Guillevic, and Piellard,
ePrint 2022/352, Proposition 5, excludes the rank-two alternative. Commit
`a8075e7` binds that result to the exact CatBLST relation and checks the
BLS12-381 seed, cofactor polynomial, and `gcd(chi(u), c2) = 1` in Lean.
`AUD-231` is validated with the published curve-family theorem and twist order
as an explicit external mathematical dependency; a foundational Lean point
count is not claimed.
Focused module and transitive root Lean
type-checks pass. This
proof-only completion ran no C++ build, Peach job, or benchmark.

The post-selection cleanup at `9573d58` removes GMP/MPIR from the live build
and source surface. Obsolete arbitrary-precision constant validation,
reference tests, benchmarks, and comparison scripts are deleted; selected
MatrixNoK arithmetic remains unchanged. The GMP-free portable-Debug build and
30 focused MatrixNoK/C API tests pass, all 784 retained files are inventoried,
and the exact-source correspondence validator accepts 20 rows and 59 hashes
after rebinding the six affected rows. `AUD-232` is validated. No benchmark,
native run, CUDA work, or BTx audit was performed.

The follow-up reachability cleanup at `b49ae14` removes the unreachable
pre-MatrixNoK `BoundedRing`, baked runtime-constant tables, runtime view, and
default Fp12 policy. Production now names the selected MatrixNoK ring and Fp12
policy explicitly while retaining the reusable arithmetic templates. The
portable-Debug production library and test executable rebuild, 30/30 focused
MatrixNoK/C API tests pass, and the constant and dispatch validators accept the
result. Exact proof-owner hashes are rebound with no selected arithmetic-body
change. `AUD-233` is validated. No benchmark, native run, CUDA work, portable
Release build, or BTx audit was performed.

The subsequent retained-file audit at `b1036f7` removes the obsolete generic
AVX RNS layer. `AUD-234` records padded full-register transfers that were unsafe
for generic non-eight-lane dimensions; `AUD-235` records the template contract
exceeding the quotient and accumulator bounds proved for the fixed engine; and
`AUD-236` records the dead/test-only surface. The selected MatrixNoK dimensions
were exactly eight lanes and were not affected. An independent scalar IFMA
accumulator oracle replaces the same-algorithm generic oracle. The portable-
Debug rebuild and 30/30 focused tests pass, as do the constant and dispatch
validators. Representative optimized cyclotomic-square code generation is
instruction-for-instruction identical. No benchmark or remote job ran.

The fixed MatrixNoK base-field owner audit is complete at `01bcec8`.
`AUD-237` routes production one-shot and batch Fp multiplication through the
performance-selected exact-bound kernels instead of the older delayed path;
the structural validator now pins that dispatch. `AUD-238` adds Debug-only
proof-domain, canonical-import, and raw-object-disjointness checks with zero
NDEBUG control-flow cost. All 31 focused portable-Debug MatrixNoK/C API tests
pass, and the optimized production-spelling probe fully inlines with no
assertion branch. Seven affected correspondence rows are rebound and all 20
AVX correspondence rows validate. No benchmark or remote job ran.

The foundational MSVC wide-integer repair is complete at `33f6a8b`.
`AUD-239` replaces the false `int128_t = int64_t` overflow-safety alias with a
real 16-byte signed carrier. `AUD-240` preserves portable constexpr products
but routes runtime 64-by-64 multiplication through `_umul128`. All 32 focused
portable-Debug tests pass. Optimized production Fp multiplication changes from
eight compatibility calls to eight inline `mul` instructions with no call or
conditional jump. The rest of `constexpr_utils.h` and the foundational bounds,
carrier, and reducer owners were left for the next audit packet; no benchmark
ran.
Follow-up `fbaa3a9` first expressed the same forced-inline policy through the
inherited AVX macro. `AUD-241` then removes that duplicate preprocessor owner at
`2d691ed` and uses CatBLST's canonical `BLST_FORCEINLINE` policy throughout the
AVX templates. A serial portable-Debug rebuild passes 39/39 focused tests and
all structural and exact-hash correspondence checks pass after rebinding. The
optimized native Fp probe retains 244 instructions, a 312-byte frame, eight
inline multiplies, zero calls, and zero branches; no benchmark or remote job
ran.

The remainder of the foundational AVX utility audit is complete at `049a8a1`.
`AUD-242` removes duplicate sizing/log helpers, the dead literal parser, both
obsolete reducer owners, and their unused reference test. The retained bound
and Montgomery regressions now exercise the selected fixed reducer against an
independent scalar oracle. `AUD-243` gives the live modular-product helper one
portable zero-modulus contract. A serial portable-Debug rebuild passes 39/39;
the dispatch and exact correspondence validators pass after rebinding. Native
MSVC Fp code generation remains 244 instructions, a 312-byte frame, eight
inline multiplies, zero calls, and zero branches. No benchmark or remote job
ran.

The paired `bounds.h` and `bounded_element.h` review is complete at `0471f82`.
`AUD-244` and `AUD-245` move the already-proved selected bound and borrow-free
wide-subtraction requirements into zero-cost C++ static assertions. All current
production expressions compile: the maximum selected bound remains 832 below
the 2048 accumulator limit, and every subtraction is offset-dominated.
`AUD-246` through `AUD-248` remove unsafe dead diagnostics and unary negation,
harden bound shifts and radix parameters, make zero construction explicit, and
document overwrite-before-read default storage. The focused MatrixNoK suite
passes 24/24 and the structural validator passes. The optimized native MSVC Fp
probe has a 312-byte frame, eight scalar multiplies, zero calls, and zero
conditional branches. No benchmark, remote job, CUDA, portable-AVX Release, or
BTx audit ran.

The complete `ring_element.h` review closes `AUD-249` through `AUD-252` at
`1e0b6e1`. Unused helpers that could attach unproved physical and cross-base
claims are gone; completion preserves the exact proof-metadata type; RNS
recasts consistently accept 64-bit endpoints; and the owning header now
documents raw observation, initialization, delayed-product ordering, and
producer obligations. The serial portable-AVX Debug target and direct carrier
test pass, as do the structural and exact-correspondence validators. No
benchmark, remote job, native build, CUDA, or portable-AVX Release build ran.

The specialized `matrixnok_fp.h` review closes `AUD-253` through `AUD-257` at
`0562a52`. The inherited `4309` ceiling is confirmed correct: it is the width
of VROOM's reviewed `[-1932,2377]` represented-multiple interval after
negative completion. Its full `40^2` wide encoding is now proved below the
exact generated `M1*M2` capacity. Negative-bound preparation is constrained,
dead runtime field-literal parsing is removed, and raw export, proof-type, and
fixed-batch contracts are explicit. The serial portable-AVX Debug build and
seven focused tests pass; all structural, constant, Lean, and correspondence
checks pass. No benchmark, remote/native run, CUDA, portable-AVX Release, or
BTx/BTE audit ran.

`AUD-258` is closed at `81f21ac`: the zero-caller VROOM test-only
`blst/avx/rns/inversion.h` has been removed. Its 425-operation chain was
independently checked and correctly computes `3p-4`; deletion addresses dead
ownership and interface defects rather than an arithmetic error. Active AVX
inversion remains in `blst/avx/fields/conversion_inversion.h`. The structural
dispatch and audit-manifest gates pass. No build, benchmark, remote/native
run, CUDA, portable-AVX Release, or BTx/BTE work ran.

`AUD-259` through `AUD-263` are closed at `dd4e4a1`. The active MatrixNoK
conversion/inversion owner is now const-correct, selected-backend-only, and
explicit about variable-time nonzero inversion. A confirmed MSVC
secret-dependent carry branch is replaced by the established carry primitive;
optimized inspection shows `adc`/`setb`, fixed-count loops, and no out-of-line
conversion wrapper. The independent validator proves the complete raw-export
machine envelope through REDC. A serial portable-Debug rebuild and 32 focused
tests pass. No benchmark, remote/native execution, CUDA, portable-AVX Release,
sanitizer, broad suite, or BTx/BTE audit ran.

`AUD-264` through `AUD-267` are closed at `13fe98a`. The selected
`FP2Ring<FpRingCandidate>` is now the sole retained MatrixNoK Fp2 arithmetic
owner. It owns the empty base-ring policy by value in one byte rather than
borrowing a potentially dangling reference. Dead prepared conversions, batch
plumbing, the comparison-only microkernel, and its private base-field hooks are
gone; exact proof metadata, provenance, initialization, timing, and raw
observation contracts are explicit. Three focused portable-Debug scalar and
completion differentials pass. The small native MSVC probe has zero calls and
only fixed-control branches. All 20 correspondence rows and 59 hashes validate
after exact rebinding. No benchmark, remote run, broad suite, CUDA,
portable-AVX Release, sanitizer, or BTx/BTE audit ran.

`AUD-268` through `AUD-273` are now closed at source commits `ce22bc2`,
`5042c1b`, and `793b542`. The selected Fp12 owner has one schedule family, one
stateless constants policy, exact conservative split bounds, direct inversion
proof/differentials, and explicit contracts/provenance. MSVC optimized
inspection confirmed and then eliminated the cyclotomic concatenation copy:
the owning wrapper fell from a 3,168-byte frame plus 1,536-byte `memcpy` to a
40-byte frame with exactly two split-helper calls and no copy. Five focused
portable-Debug tests pass, the focused Lean modules type-check, and all 21 AVX
correspondence rows/61 hashes validate. No benchmark or broad/resource-heavy
run was performed.

The complete selected AVX final-exponentiation review is recorded at source
commit `835ffbb`. Its formulas remain supported by the scalar BLST schedule and
accepted Lean evidence, but `AUD-274` and `AUD-275` confirm a cross-toolchain
hot-path ownership defect: MSVC emits at least 305 loop-level 1,536-byte copies
per final exponent and a 38,528-byte hard-part frame, while Clang 18 retains a
35,456-byte frame and the same repeated-copy shape. `AUD-276` records missing
provenance and narrow contracts. The findings are open pending an explicitly
reviewed ping-pong/accumulator repair. Only compile-time code-generation probes
ran; no benchmark, remote job, broad build, CUDA, portable-AVX Release, or
BTx/BTE work ran.

`AUD-274` through `AUD-276` are validated at `35516ff`. Final exponentiation
now owns its distinct-output square locally, uses fixed ping-pong chains, and
retains only a bounded accumulator set. MSVC hard-part stack falls from 38,528
to 12,416 bytes with about 15.6 KiB peak nested stack; Clang 18 falls from
35,456 to 14,720 bytes. Both compilers emit zero repeated-square `memcpy` calls,
reducing nested chain copies from at least 305 to 25 multiply-result boundaries.
Five focused portable-Debug differentials and the three exact final-exponent
Lean modules pass; all 21 correspondence rows and 61 hashes validate after the
source rebind. Public C contracts and provenance are now explicit. A later
serial Peach comparison measures isolated final exponentiation 3.82% faster
(79.176 to 76.148 us median) and complete pairing 1.95% faster (136.713 to
134.045 us). No broad/resource-heavy run was performed.

`AUD-277` is validated at `dd68a0e`. Prepared scalar and AVX Miller evaluation
now preserves the pairing identity for both G1 and G2 infinity, and generic AVX
prepared batches implement identity filtering plus compact continuation. The
already-prefiltered middle-product path remains branch-free in Release. The P4
correspondence row now distinguishes finite affine formula correctness from
public representation dispatch; new Lean theorems prove both observed-
infinity branches and the exact 68-line identity-table schedule. Three focused
identity regressions and four retained finite Miller tests pass. `AUD-278`
through `AUD-282` remain open for the rest of the AVX GT owner review. No
benchmark, remote job, CUDA, portable-AVX Release, or broad runtime suite ran.

`AUD-278`, `AUD-279`, and `AUD-281` are validated at `65206c9`. Native
prepared-line tables now value-construct as one exact zero-fill with no
prepared reductions or per-line initialization loop; public variable-time GT
scalars are no longer wiped; and the selected AVX GT owner has explicit
provenance plus batch, line-table, timing, allocation, alias, conditional-bit,
and Frobenius contracts. Four focused portable-Debug tests and the affected
middle-product object pass. `AUD-282` remains open for reusable scratch during
the middle-product audit, although its current six-buffer allocation contract
is now explicit. No broad build, CUDA, portable-AVX Release, sanitizer, Apple
CI, or BTx/BTE work ran.

`AUD-280` is validated with dense production doubling retained. Correctness at
`4e7c618` showed the cyclotomic candidate is semantically admissible, but one
valid Release process on Peach at `fa26da3` measured Clang 21.1.8 native
AVX-512/IFMA medians of 326 ns dense versus 366 ns cyclotomic and 2.530 us
versus 2.895 us for eight doublings. The candidate regressed 12.3% and 14.4%.
Commit `7295176` removes the losing benchmark/test scaffolding. An earlier
sequential process lacking `-DNDEBUG` was discarded; no process remained after
the measurement.

Phase 17 middle-product safety and contract packet is validated at `5e14afc`.
`AUD-283` rejects zero, oversized, and overflow-shaped lengths before span
construction or allocation, with a `2^30` FFT and `2^29` optimized-middle-
product bound derived from the signed FFT index domain. `AUD-286` documents the
exact diagonal pairing product and subgroup, GT, timing, allocation, exception,
failure-publication, ownership, concurrency, and alias contracts across C,
C++, and Rust. `AUD-290` fixes scalar all-zero infinity crossing into MatrixNoK
without canonicalization; the prior bridge value was accepted by addition but
could fail AVX multiplication inside a group FFT. `AUD-291` makes every
backend-explicit algorithm test call its named scalar or AVX implementation.

MSVC 19.50 scalar-default plus portable-AVX Debug passes seven focused sizing,
C-boundary, scalar-differential, exact-alias, and G1/G2 bridge tests, plus the
corrected explicit AVX FFT comparison. Rust middle-product tests pass 8/8 and
both documentation ledgers pass. Old portable and scalar executables preserve
the AUD-290 before/backend-isolation evidence. `AUD-282`, `AUD-284`, `AUD-285`,
and `AUD-287` through `AUD-289` remain open for reusable final-exp scratch,
precomputation, backend selection, direct GT conversion, and cleanup. No
benchmark, remote run, CUDA, sanitizer, portable-AVX Release, full suite, or
BTx/BTE audit ran.

Phase 17 precomputation and boundary cleanup is validated at `f39063f`.
`AUD-284` replaces one Fr inversion per inverse root with the exact reversed
ordered-power table. `AUD-285` stores one forward and one inverse truncated FFT
plan in the immutable B precompute and reuses them across B construction and A
evaluations; B0 also transforms directly in final storage. `AUD-288` replaces
the 576-byte AVX GT serialize/decode boundary with the direct typed Fp12 bridge.
`AUD-289` removes dead selectors/instrumentation, consolidates the G1 boundary
conversion, fixes formatting residue, and adds Category Labs, Apache-2.0, and
upstream lineage. MSVC 19.50 scalar-default plus portable-AVX Debug builds with
two workers and passes eight focused root, plan, sizing, and C-boundary tests.
`AUD-282` and benchmark-gated `AUD-287` remain open. Frozen BTx rebuilt as a
dependency and emitted its existing nodiscard warning, but was neither changed
nor audited. No benchmark, remote run, CUDA, sanitizer, portable-AVX Release,
or full suite ran.

`AUD-282` is validated at `30c7558`. AVX batch final exponentiation now takes
explicit caller-owned scratch for all six decomposition and batch-inversion
buffers. The reusable multipair middle-product workspace retains that scratch
across evaluations; the one-shot update path owns an explicit local instance.
Once capacity is established, smaller or equal final-exp batches perform no
storage allocation, while growth still occurs before output publication and
retains the existing exception boundary. Portable-AVX Debug passes direct
batch-versus-individual, zero-length, size-oscillation, repeated-workspace, and
existing multipair-versus-FFT cases at sizes 3, 5, 9, and 16. All 21 AVX
correspondence rows and 61 mapped hashes revalidate after the ownership-only
`GT.h` rebind. Only benchmark-gated `AUD-287` remains open in this middle-
product packet. No benchmark, remote run, CUDA, sanitizer, portable-AVX
Release, full suite, or BTx/BTE audit ran.

`AUD-287` is validated at `e3e1a45`. The old `n <= 24` Auto cutoff did not
represent the power-of-two FFT cost jumps and its benchmark did not time the
one-shot production ownership model. An isolated non-shipping target now times
forced scalar/AVX FFT and multipair calls without BTx/BTE linkage. Sequential
Peach runs used Ryzen 9 7950X CPU 2, Clang 21.1.8, separate scalar-default and
AVX-default Release configurations, and five repetitions. Scalar selects
multipair for `n <= 31`, `33..49`, and `65..76`; AVX extends the second range
through 50. The constexpr policies add no allocation, timing probe, or runtime
CPU dispatch. Portable-Debug range and actual dispatch tests pass, and both
native Release configurations rebuild cleanly after `4dcb4fa` and `bdefea3`
mark the inactive scalar/AVX GT bridge overloads `[[maybe_unused]]`. All
findings recorded by the current middle-product pass are now validated. The
thresholds remain deployment-CPU evidence rather than a universal
microarchitecture guarantee. No concurrent benchmark, CUDA, sanitizer,
portable-AVX Release, full suite, or BTx/BTE audit ran.

Phase 17 core FFT Packet 17C is remediated at
`9d1e2bd60e610e88b44c1be36daea7be793898c0`. `AUD-292` removes
mathematically incorrect, uncalled radix-4 and
radix-8 experiments plus stale registrations. `AUD-293` and `AUD-295` define
the retained radix-two root direction, size, alias, zero-tail, destructive
buffer, plan immutability, and output-publication contracts. `AUD-294` rejects
unsupported plan sizes before allocation and replaces overflowable inverse
range addition with subtraction-form validation. `AUD-296` moves composed
timing-key allocation inside the timer's no-throw barrier. Unused reversal-NTT
helpers, the unused default coefficient context, and disabled removed-API test
blocks are gone.

MSVC 19.50 scalar-default plus portable-AVX Debug rebuilt `blst_tests` with two
workers. Eight focused `FftApi.*`, `FftCore.*`, and `Timer.*` tests pass, and
the standalone simple/composed allocation-failure probe exits zero. Frozen
BTx/BTE compiled unchanged and confirmed that it uses only retained ordinary,
polynomial, and truncated transforms. No benchmark, CUDA, full suite,
sanitizer, native AVX-512, remote, or Apple run was performed.
Phase 17 Packet 17F is remediated through implementation commit `1821c58` and
documented interface commit `baf3e45`. `AUD-297` through `AUD-300` and `AUD-302`
are validated by threshold-crossing oracles, fixed-shape reconstruction, flat
tree ownership, remainder-only evaluation, checked reduced scratch, and Debug
invariants. `AUD-301` remains open because the quadratic suffix crossover needs
explicit benchmark approval. Validation follow-up `AUD-303` is fixed at
`2846e0a`. No benchmark or BTx/BTE change is part of this packet.

Phase 17 Packet 17G completes the two small FFT support owners at implementation
commit `0fc26d2`. `AUD-304` adds allocation-before-use validation to the root
builder and zero-Release-cost exact shape assertions to native domain helpers.
`AUD-305` narrows the coefficient context to G1/G2/GT-by-public-Fr and pairing,
documents subgroup and variable-time root contracts, and preserves direct
compile-time dispatch and fully overwritten output temporaries. Scalar and
portable-AVX-selected Debug focused checks pass; `domain.h` compiles standalone.
No benchmark, CUDA, full suite, sanitizer, native AVX-512, remote, Apple, or
BTx/BTE audit/change was performed.

Phase 17 Packet 17H completes `blst/fft/primitive_root.h` and its misspelled
compatibility shim at implementation commit `4a6e951`. `AUD-306` enforces exact
power-of-two root order in Release before output publication. `AUD-307` replaces
guarded runtime table construction with 33 independently regenerated immutable
Montgomery constants; optimized lookup has one indexed 32-byte load and no
initialization path. `AUD-308` removes randomized generator discovery, trial
factorization, dead helpers/macros, and transitive dependencies while retaining
an explicit-generator generic API. Scalar Debug/Release focused tests pass 5/5,
the independent validator and corruption self-test pass, and the header compiles
standalone at MSVC `/W4`. Frozen BTx/BTE compiled unchanged through the shim.
No benchmark, CUDA run, full suite, sanitizer, native AVX-512, remote, or Apple
run was performed.

Phase 17 Packet 17I completes the `blst/fft_c_api.cpp` line-by-line hand review
at `f4063cc`. No arithmetic or boundary-safety defect was found. `AUD-309` and
`AUD-310` are validated: all seven C/safe-Rust FFT and interpolation operations
now have precise representation, size, allocation, failure-publication, alias,
timing, scratch, and concurrency contracts; the four fallible FFT/domain C
declarations are nodiscard. DOC-026 raises the incremental documentation ledger
to 26 audited/focused families. Strict MSVC C11/C++20 headers and five focused
scalar Debug/Release boundary tests pass. No benchmark, CUDA, portable-AVX
Release, full suite, sanitizer, native AVX-512, remote, Apple, or BTx/BTE
audit/change was performed.

Phase 17 Packet 17J completes the `blst/fft/middle_product.h` owner and
`blst/middle_product_c_api.cpp` C-boundary reviews at `9e4ec59`. `AUD-311`
moves all scalar sparse-batch growth before output publication. `AUD-312`
removes selected-backend conversion vectors at compile time while retaining
explicit cross-backend conversions. Contracts now cover immutable expert
precomputes, invalid enum values, live C handles, and exact aliasing with either
combine input. Focused scalar-default portable Debug passes 22/22, scalar-only
Release passes 14/14, AVX-default portable Debug passes 2/2, and strict MSVC
C11/C++20 header checks pass. Arithmetic and dispatch are unchanged, so no
benchmark was run. No CUDA, portable-AVX Release, full suite, sanitizer, native
AVX-512, remote, Apple, or BTx/BTE audit/change was performed.

Phase 17 Packet 17K completes `blst/avx/fft/gt_context.h` at `ec18911` with
proof/evidence rebind `52d71c1`. `AUD-313` through `AUD-315` validate exact
GT-only/public-Fr/lifetime contracts, remove the dead canonical-limb cache and
allocation, and derive inverse scale from the owned roots. The Release hot path
is unchanged. Four focused portable-AVX Debug tests, the MatrixNoK dispatch
gate, and all 21 proof rows/61 hashes pass. No benchmark, CUDA, portable-AVX
Release, full suite, sanitizer, remote, Apple, or BTx/BTE audit/change ran.

Phase 18 expert AVX remediation is implemented in the worktree. `AUD-316`
through `AUD-323` cover final-exponent object identity, generic Fp12 routing,
the safe Rust GT precondition, opaque equality, fallible allocation, complete
77/77 raw declarations, externally used GLV hook retention, and empty GT-MSM
identity. Focused portable-AVX Debug C++ tests pass 2/2; Rust AVX tests pass 7/7
and the non-AVX equality test passes 1/1. Strict C11/C++20 header parsing,
MatrixNoK dispatch, audit manifests, and 21 proof rows/61 hashes pass. Exact
the exact implementation is recorded at `e6163ce`; no benchmark or
broad/platform sweep ran.

Phase 18 Rust generation and binding hardening is complete through
`835b97d802948fabb1e88383a975df07d98cf0ee`. `AUD-324` through `AUD-333` are
validated. The scalar raw binding now regenerates exactly with bindgen 0.65.1,
Clang/libclang 18.1.3, pinned Rust 1.65 rustfmt, explicit repository formatting,
public-name allowlists, distinct opaque C tags, and fail-closed normalization.
Two output locations produced identical bytes; the 253/253 declaration gate,
transform fixture, 13 generated tests, and 13 pairing/uniqueness wrapper tests
pass. A path-scoped CI lane repeats the comparison without a native product
build. No benchmark, CUDA, broad native rebuild, frozen BTx/BTE audit, native
AVX run, push, or hosted CI wait was performed.
