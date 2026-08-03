# Phase 12 - Audit Platform and Initial Baseline

## State

Phase 12 is in progress. Automated analysis is a standing workstream through
Phase 19; this phase establishes its reproducible initial platform and triages
the first repository baseline.

## Packet 12A - CodeQL and CI provenance

The first packet is limited to CI, dependency-update, and audit-record changes:

- CodeQL covers retained C/C++ and Rust only;
- C/C++ extraction manually builds portable scalar and portable AVX source
  configurations in Debug mode;
- Rust uses CodeQL's supported no-build extraction mode;
- the C/C++ build target is `blst_static`, so benchmarks, frozen BTx/BTE,
  experimental CUDA, and removed bindings are not part of the audit scan;
- portable-AVX Debug remains in routine CI and automated analysis, while
  portable-AVX Release is excluded because optimized compilation materially
  slows the audit cycle;
- retained GitHub Actions are pinned to immutable full commit SHAs; and
- Dependabot proposes grouped weekly updates for those pins.

Frozen BTx/BTE may continue to compile and execute in the repository's
pre-existing general regression suite. It remains outside the audit and no new
analyzer, sanitizer, or fuzz work targets it.

## Pinned action revisions

| Action | Release | Commit |
| --- | --- | --- |
| `actions/checkout` | `v7.0.1` | `3d3c42e5aac5ba805825da76410c181273ba90b1` |
| `actions/cache` | `v6.1.0` | `55cc8345863c7cc4c66a329aec7e433d2d1c52a9` |
| `github/codeql-action` | `v4.37.4` | `f205ea1c3313d32999d8d6a48b4f6530d4437b38` |
| `ilammy/msvc-dev-cmd` | `v1.13.0` | `0b201ec74fa43914dc39ae48a89fd1d8cb592756` |

These revisions were resolved from their upstream release tags on 2026-07-31.

The first hosted execution completed both language analyses but could not
upload their SARIF because this private repository does not have GitHub Code
Security enabled. The C/C++ job completed in 8:54, including 2:59 for both
scalar and portable-AVX Debug builds. The Rust job reached SARIF generation in
4:23 and reported five files with extraction errors for local triage. These
were service-boundary failures, not analyzer crashes.

Until Code Security is enabled, `.github/workflows/codeql-analysis.yml` is
manual and guarded by the repository variable `CODEQL_ENABLED=true`. The
pinned manual procedure in `CODEQL_LOCAL.md` produces local SARIF without
using an unlicensed CI artifact workaround.

## First local CodeQL baseline

The pinned CodeQL 2.26.2 bundle was installed and its Linux archive verified
against SHA-256
`cb361567fa1bdb9d322da4240f621b36f245e4d7bb97db3c3a2ad7f743c8e8e7`.
Both local suites ran against commit
`e026202ed606ce29ca94f5ceefa55913ebb249de`.

- C/C++ run `e026202ed606-20260801T005352Z` completed in 8:32. It reported
  958 raw results: 921 commented-out-code notices and 37 results across the
  other 11 rules. It scanned 88 of the repository's 318 C/C++ files because
  the trace intentionally builds only the retained `blst_static` production
  target. CodeQL also reported 68 partial extraction diagnostics: 67 in
  `blst/avx/curve/ec.h` and one in `blst/avx/rns/bounded_ring.h`. The portable
  AVX result set is therefore incomplete until this extractor limitation is
  resolved or covered by another analyzer.
- The first Rust attempt exposed ignored Cargo output to build-free
  extraction. The runner was corrected to snapshot only tracked and
  non-ignored working-tree files. Clean run
  `e026202ed606-20260801T011000Z` completed in 1:41, scanned all 14 Rust source
  files, and produced no alerts. Thirteen files extracted without error; one
  file, `bindings/rust/src/lib.rs`, retained `format!` macro-expansion
  diagnostics and is a recorded coverage limitation.
- The clean Rust source snapshot contained 606 files, contained no `out/` or
  `bindings/rust/target/` entries, and recorded archive SHA-256
  `71fdb7aaa2ac887818d4ced644c10d00dded53037eaf30f039781faf72be9ed9`.

The SARIF and databases remain ignored local evidence. Baseline dispositions
and proposed findings `AUD-024` through `AUD-030` are recorded in
`CODEQL_BASELINE_TRIAGE.md`. The user approved the first correctness packet;
`AUD-024` and `AUD-025` are now validated. Review of the FFT boundary also
identified and closed prepared-domain authentication and size-growth findings
`AUD-031` and `AUD-032`.

## First correctness remediation packet

The initial approved packet fixed findings before expanding the audit tooling:

- `Point::sign2` now copies the complete optional affine result, with direct G1
  and G2 tests for combined, affine-only, and serialized-only outputs;
- FFT-domain interpolation rejects sample counts larger than the domain before
  partitioning scratch or writing output;
- scratch and adjacent polynomial size growth use checked arithmetic;
- prepared domains now contain an encoded size tag and require exact storage,
  preventing accidental reuse for a different transform without adding a
  field inversion to the prepared transform path; and
- the radix-two butterfly loop stops before its terminal signed shift.

Focused validation matched `VALIDATION_POLICY.md`:

- Windows scalar builds passed the two native-signature and three FFT-boundary
  regressions;
- Linux GCC ASan passed all five focused C++ regressions;
- Rust FFT module tests passed 11/11 in both default and portable release
  configurations; and
- `git diff --check` passed apart from local line-ending notices.

No benchmark, broad AVX, CUDA, BTx/BTE, or full-suite run was performed for
this focused correctness packet.

## Second CodeQL hardening packet

The second approved packet closes the remaining low-severity actionable
baseline findings without changing algorithms or ownership:

- scalar MSM and bulk-addition template headers now tolerate direct, repeated,
  dependency-repeated, and wrapper-plus-direct inclusion;
- the hot AVX G2 psi helper and three construction-only RNS inputs use const
  references while destination objects retain the same owned values;
- the prepared-root scope rejects temporary roots and all guard copies/moves,
  with explicitly non-throwing scope entry and restoration; and
- prepared-root lookup no longer orders or subtracts unrelated pointers with
  built-in operators. It uses total ordering, constant-time address/alignment
  validation, and exact pointer identity before selecting a compiled exponent.

The prepared-root lookup remains inline, allocation-free, and O(1). It adds no
scan, dynamic dispatch, callback, or synchronization to the FFT butterfly path.
The const-reference changes remove potential copies without adding indirection
beyond the references already required to read each input.

Focused validation:

- dedicated scalar and portable AVX-selected builds passed;
- the repeated scalar-header translation unit compiled and its focused test
  passed;
- compile-time prepared-root construction, noexcept, copy, and move assertions
  passed;
- prepared lookup, unrelated same-value fallback, nested exception restoration,
  full/truncated inverse FFT, and scalar differential checks passed;
- AVX G2 Pippenger/scalar multiplication matched its naive oracle; and
- no benchmark or full-suite sweep was run for these signature and guardrail
  changes.

## Packet 12A review validation

- `scripts/audit/check_audit_manifest.ps1`: passed with 607 files and 63
  upstream commits accounted for.
- `git diff --check`: passed.
- `actionlint` 1.7.12: passed for all retained workflow files.
- WSL Ubuntu with GCC 13.3 and CMake 3.28.3: both proposed CodeQL C/C++
  configurations completed CMake configuration.
- The portable-AVX Debug configuration built `blst_static` successfully with
  two workers. On the Windows-mounted WSL filesystem it took 7:59.79 and
  reported 1,534,072 KiB maximum resident memory. This is not treated as a
  hosted-Linux timing prediction; the first CodeQL run will decide whether the
  AVX extraction remains in the ordinary change loop or moves to a manual or
  subsystem gate.

No tests, benchmarks, sanitizer sweeps, or unrelated builds were run for the
CI-only Packet 12A review itself. The later correctness remediation evidence is
recorded separately above.

## Remaining Phase 12 packets

1. Review and resolve the sanitizer-discovered scalar configuration finding
   `AUD-041` without changing optimized default code generation.
2. Add initial bounded public byte-boundary libFuzzer targets and seed corpora.
3. Complete diagnostic triage and Phase 12 exit-gate evidence.

## AUD-038 commented-code cleanup

The CodeQL baseline's 921 commented-code results were split into two classes.
The 171 results outside generated addition-chain headers represented abandoned
alternate implementations, stale debug fragments, duplicate wrappers, or
interface descriptions written as commented source. They were removed or
rewritten across 30 retained production files without replacing any active
algorithm. The unused, wholly commented-out `test/BasicFp.h` was also deleted.

The 750 results in `recip-addchain.h` and `sqrt-addchain.h` are retained. They
record the generated reciprocal and square-root schedules and therefore serve
as arithmetic provenance; suppressing or deleting them would lose useful
review evidence while no repository-owned generator is available.

Focused validation kept the cycle proportional to the cleanup:

- the scalar Release library and tests built, and 167 field, curve, pairing,
  FFT, middle-product, backend, and header regressions passed;
- the portable-AVX Debug library and tests built successfully;
- a portable-AVX Release rebuild was stopped after its three-minute focused
  budget while still compiling, so it is not claimed as evidence; and
- CodeQL run `01eee6b03950-20260802T223754Z` completed both production
  profiles and reported exactly the retained 750 comments in each, with zero
  commented-code results elsewhere.

No benchmark, full-suite sweep, CUDA, or BTx/BTE audit was performed.

Following this packet, portable-AVX Release was removed from routine CI and
audit validation because optimized compilation exceeded the focused cycle
budget. Portable-AVX Debug remains available in CI, CodeQL, and Clang analyzer
entry points. Native AVX-512 remains the supported production-performance gate,
and the narrow native-versus-emulation compatibility target remains available
when the emulation layer itself changes.

## Split CodeQL and Clang analyzer confirmation

At commit `6bee5ca43f94`, the local CodeQL runner created independent scalar and
portable-AVX databases. Scalar extraction was clean; AVX retained the same 68
known dependent-template diagnostics, now isolated rather than mixed into a
combined database. No actionable `AUD-026` through `AUD-029` alert remains.

Ubuntu Clang 18.1.3 `scan-build -analyze-headers` then compiled and analyzed
all 17 retained scalar and all 18 retained portable-AVX source compilations.
It traversed the AVX group-operation templates CodeQL could not model. Fifteen
raw results reduce to two proposed low-severity cleanups (`AUD-035` and
`AUD-036`) plus documented analyzer false positives. Exact commands, raw
counts, dispositions, and the detached-wrapper logging limitation are in
`CLANG_ANALYZER_BASELINE_TRIAGE.md`.

No benchmarks, full test sweep, CUDA, BTx/BTE, assembly, or unsupported
language bindings were included.

## AUD-037 selected-vector separation

The compatibility layer now has three non-overlapping domains:

- `_mm512_*` and native intrinsic types are the hardware implementation;
- `_cc512_*` and lane-array `cc512i`/`cc256i`, defined only with
  `BLST_AVX512_PORTABLE`, are the pure C/C++ emulation; and
- `_vv512_*` plus `vv512i`/`vv256i` are the compile-time-selected production
  implementation. Native wins when both features are configured.

The former `BLST_AVX_NATIVE` layout switch, ambiguous unprefixed operations,
and uppercase `V*` adapter layer were removed. AVX algorithms and layout
helpers now use `vv` exclusively. Native wrappers call intrinsics directly;
the only representation bridges are explicitly named and used by the direct
compatibility test. Portable definitions are compiled whenever
`BLST_AVX512_PORTABLE` is enabled, including combined builds.

Focused evidence:

- neither, native-only, portable-only, and combined configurations compiled;
- the lightweight combined target compiled only the emulation source and
  direct intrinsic tests, statically selected native `vv`, retained a 64-byte
  standard-layout `cc512i`, and discovered 9 tests;
- portable field/curve integration passed 28/28;
- optimized GCC probes for `madd52hi` and `madd52lo` plus add emitted exactly
  the same instructions through `_vv512_*` as direct `_mm512_*` calls, with no
  calls, conversions, or memory round-trips; and
- structural checks found no remaining `BLST_AVX_NATIVE`, legacy portable
  operation names, ambiguous unprefixed operations, uppercase `V*` adapters,
  production `cc` use, or legacy conversion names.

Peach supplied the missing native execution evidence: its Ryzen 9 7950X
exposes AVX-512 IFMA, and supported GCC 15.2 built the Release combined target
without warnings. CTest discovered and passed all 9 direct native-versus-
emulation tests, including the dedicated IFMA cases. `AUD-037` is validated,
and this supplies the restored direct IFMA evidence for `AUD-036`. GitHub
Actions run `30754724252` supplied that finding's final gate at current commit
`0b7ec47`: native Apple Silicon configured with the experimental ARM override
off, AppleClang 17 built the scalar C/C++ target without the former
conflicting-attribute warning, all 336 C/C++ tests passed, and the Rust default
and full-binding matrix passed. `AUD-036` is validated. No benchmark was run
because native AVX code generation was unchanged.

## clang-tidy baseline

LLVM clang-tidy 18.1.3 ran report-only against separate scalar and
portable-AVX Debug compilation databases. Portable-AVX Release, benchmarks,
BTx/BTE, CUDA, tests, and assembly were excluded. After canonical-location
deduplication, scalar contributed 47 diagnostics and portable AVX added 26,
for 73 unique locations across nine rules.

Seventy-two are bounded conversions, fixed-use macros, ABI-preserving enum
sizes, intentional assignments or signed JSF digits, debug output, or the
required scratch-alignment pointer round-trip. The remaining diagnostic opened
`AUD-039`: an unconstrained fixed-window template combined an inclusive table
bound with an 8-bit loop and booth index. The approved narrow compile-time
contract preserves the hot loop and is validated by boundary assertions,
portable G1/G2 scalar differential tests, and identical native Release code
generation. No clang-tidy fix-it was applied. Exact counts and dispositions are
in `CLANG_TIDY_BASELINE_TRIAGE.md`.

## Clippy baseline

Rust/Cargo/Clippy 1.94.1 ran report-only over the library and build script in
five Debug feature profiles: default, full bindings, portable scalar, AVX
(portable AVX on Windows), and full bindings plus secret serialization.
Benchmarks, tests, examples, frozen BTx/BTE, experimental CUDA, and unsupported
targets were excluded. Portable AVX Release was not compiled.

The default, portable, and AVX profiles each produced 63 canonical broad-lint
diagnostics. Full and full-plus-serde each produced 133 after deduplicating the
two signature-orientation macro expansions. Pedantic/nursery results are
triaged as documentation/style debt, generated ABI fields, deliberate FFI
pointer operations, an intentional mutable aliasing boundary, or transparent
wrapper conversions reserved for the Phase 18 unsafe review.

Three diagnostics belonged to `clippy::all`. The exact remainder test is
retained for compiler compatibility. `AUD-040` resolved the other two: Rust AVX
GT MSM now exposes only zero-length `MaybeUninit` capacity to C++, and its
scratch element count uses overflow-free quotient/remainder arithmetic. Focused
portable-AVX Debug tests passed 4/4; default/full broad Clippy counts fell to
61/131 and `clippy::all` retains only the declined remainder suggestion. Exact
commands, counts, and dispositions are in `CLIPPY_BASELINE_TRIAGE.md`.

## Initial ASan/UBSan baseline

The dedicated `blst_sanitizer_boundary_tests` target keeps focused sanitizer
cycles separate from the full test executable and frozen BTx/BTE. On x86-64
Ubuntu 24.04 with Clang 18.1.3, the scalar portable-dispatch profile selected
123 boundary tests: 113 passed and the 10 AVX-only cases skipped as expected.
ASan plus leak detection completed in 73.970 seconds with no report; UBSan
completed in 333.152 seconds with no report. Benchmarks, CUDA, AVX-512, and
portable-AVX Release were not built.

The reusable runner snapshots the source into Linux-native temporary storage,
pins the compiler version and configuration, retains ignored evidence under
`out/audit/sanitizers`, and requires a nonzero passing test count. Two expensive
tests are deferred to subsystem/fuzz gates rather than every focused cycle.

The baseline exposed proposed `AUD-041`: CMake advertises `BLST_ASM=OFF`, but
native C++ field templates then reference undeclared assembly names; enabling
assembly while disabling ADX also fails unconditionally on x86-64. The clean
baseline therefore retained assembly and ADX while selecting portable scalar
dispatch. Generated assembly remains outside sanitizer instrumentation. Exact
commands, evidence identifiers, exclusions, and limitations are recorded in
`SANITIZER_BASELINE_TRIAGE.md`. No production dispatch fix was applied.

## AUD-035 scalar window contract

The retained `mult_window_n` allocation mode now has an explicit status result
that preserves ordinary call sites which ignore a return value. It rejects
zero sizes, missing first point/scalar
pointers, byte-count overflow, checked table-row multiplication overflow, and
allocation failure before changing the output. It no longer prints or relies
on a debug-only assertion. Caller-populated and precomputed tables plus the
compact later-null pointer convention remain supported.

The default allocator and deterministic failure allocator are selected at
compile time; the production path contains no allocator object or runtime
dispatch. Dedicated G1/G2 tests cover all three table modes, compact pointers,
every rejection, overflow, and allocation failure. The four new tests plus six
existing G1 GLV/G2 GLS multiplication and signing tests passed 10/10 under
MSVC Release, Linux ASan with leak detection, and Linux UBSan. In an optimized
GCC 13 scalar object, both production `mult_glsv` functions remained present
while `mult_window_n`, `malloc`, and `free` had no symbol or reference. No
benchmark was run.
