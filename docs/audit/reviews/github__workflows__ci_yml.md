# File Review: `.github/workflows/ci.yml`

## Record

| Field | Value |
| --- | --- |
| Path | `.github/workflows/ci.yml` |
| Phase | Phase 19, Packet 19B |
| Reviewed commit | `0c7865082575e15f779debab8186bd0eca0eea7a` |
| Reviewer | Codex audit team |
| Date | 2026-08-12 |
| Disposition | `in_scope_support` |
| Subsystem | Supported CMake and Rust CI |
| Profiles | scalar C++, runtime-dispatched assembly, Apple Silicon, Rust |

## Purpose and trust boundary

This workflow is the ordinary supported-platform correctness gate. It consumes
repository source and dependency registries on pull-request runners, so a green
result must correspond to the checked-out commit and audited build contracts.

## Contracts and invariants

- Pull requests and master pushes run; arbitrary feature pushes do not duplicate
  an open pull request's full matrix.
- CMake uses exactly `linux-scalar`, `linux`, `macos-apple-silicon`, and
  `x64-Release`.
- Apple jobs verify both the GitHub architecture context and `uname` as ARM64.
- Benchmarks, CUDA, portable AVX, and native AVX are absent.
- Rust caches never contain `Cargo.lock`; all Cargo invocations are locked.
- External actions use immutable commits and do not retain credentials.

## Checklist evidence

| # | Audit question | Result |
| --- | --- | --- |
| 1 | Trigger/concurrency | Master/PR/manual with cancellation by ref. |
| 2 | Permissions | Workflow token is `contents: read`. |
| 3 | Matrix correctness | Four distinct audited production presets; duplicate Ubuntu row removed. |
| 4 | Pure scalar coverage | Linux `linux-scalar` disables assembly and vector ISA. |
| 5 | Assembly coverage | Linux/Windows and supported Apple use runtime-dispatched assembly. |
| 6 | AVX/CUDA scope | Native AVX is hardware evidence; portable AVX and CUDA are absent here. |
| 7 | Build/test | Serial workflow jobs build `blst_tests`; CTest reports failures. |
| 8 | Resource bounds | CMake 45 minutes, Rust 60 minutes; four build/test workers per job. |
| 9 | Runner stability | Ubuntu 24.04, Windows 2025, macOS 15 ARM64. |
| 10 | Cache correctness | Registry/Git/target only, keyed by OS, arch, and lock hash. |
| 11 | Dependency resolution | Six Cargo test/smoke commands use `--locked`. |
| 12 | Unsafe expression use | Matrix values are static; no PR title/body interpolation into shells. |
| 13 | Action integrity | Checkout/cache/MSVC commits resolve to stated release tags. |
| 14 | Residual contract | Rust 1.65/MSRV is unresolved as `AUD-351`. |

## Findings and validation

`AUD-347` through `AUD-350` are validated. `AUD-351` remains open. Local
actionlint, preset parsing, action provenance, Cargo metadata, and structural
policy checks pass. Status remains `in_review` until the exact hosted matrix
runs; no local full suite was started for this workflow-only review.

