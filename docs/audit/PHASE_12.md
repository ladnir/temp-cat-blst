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
Actual CodeQL results and hosted-runner timings remain pending until the
reviewed workflow is committed, pushed, and executed.

## Packet 12A review validation

- `scripts/audit/check_audit_manifest.ps1`: passed with 603 files and 63
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

No tests, benchmarks, sanitizer sweeps, or unrelated builds were run for this
CI-only packet.

## Remaining Phase 12 packets

1. Clang Static Analyzer/CodeChecker, clang-tidy, and Clippy baselines.
2. Focused portable-scalar ASan/UBSan audit targets.
3. Initial bounded public byte-boundary libFuzzer targets and seed corpora.
4. Complete diagnostic triage and Phase 12 exit-gate evidence.
