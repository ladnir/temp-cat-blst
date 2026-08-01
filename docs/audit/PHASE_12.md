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
- portable AVX remains in the initial job on a trial basis and may be deferred
  if its hosted-runner compile cost materially slows the audit cycle;
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

1. Clang Static Analyzer/CodeChecker, clang-tidy, and Clippy baselines.
2. Focused portable-scalar ASan/UBSan audit targets.
3. Initial bounded public byte-boundary libFuzzer targets and seed corpora.
4. Complete diagnostic triage and Phase 12 exit-gate evidence.
