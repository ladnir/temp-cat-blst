# Audit Status

Last updated: 2026-08-06

Waves 1 and 2 of the parallel file audit are integrated through implementation
commit `314c186`. Wave 1 closes the scalar C-boundary lifetime and public
type/contract packets. Wave 2 records `AUD-198` through `AUD-207`, repairs the
approved scalar/extension-field findings, and extends `AUD-060` so portable
conditional negation matches assembly for every nonzero flag. Focused scalar
Release evidence passes 7/7 and optimized small-shift/cneg inspection passes;
`AUD-201` remains open for inverse timing evidence. No benchmark, CUDA, BTx,
sanitizer, Apple job, or broad runtime sweep ran.

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

Findings: 147 total: 120 validated, 16 fixed pending final evidence, four accepted
risks, four need a design decision, and three remain open. `AUD-035` is validated with preserved
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
`AUD-161` separately records that the existing expert AVX batch predicate has
caller-controlled soundness and is not yet a safe general batch-validation
interface. `AUD-162` and `AUD-163` remove the empty pairing facade and unused
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
