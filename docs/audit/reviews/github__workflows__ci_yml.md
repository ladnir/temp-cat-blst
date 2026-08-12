# File Review: `.github/workflows/ci.yml`

## Record

| Field | Value |
| --- | --- |
| Path | `.github/workflows/ci.yml` |
| Phase | Phase 19, Packet 19B |
| Reviewed commit | `a71068a47c361c9a1edae4a706b75435023756f6` |
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
- A bounded Rust 1.78 job checks the combined supported feature surface.
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
| 8 | Resource bounds | CMake 45 minutes, Rust 60 minutes, MSRV 30 minutes; MSRV uses two workers. |
| 9 | Runner stability | Ubuntu 24.04, Windows 2025, macOS 15 ARM64. |
| 10 | Cache correctness | Registry/Git/target only, keyed by OS, arch, and lock hash. |
| 11 | Dependency resolution | Six current-toolchain commands and one Rust 1.78 check use `--locked`. |
| 12 | Unsafe expression use | Matrix values are static; no PR title/body interpolation into shells. |
| 13 | Action integrity | Checkout/cache/MSVC commits resolve to stated release tags. |
| 14 | Residual contract | Exact hosted MSRV and platform execution remain. |

## Findings and validation

`AUD-347` through `AUD-351` are validated. Local actionlint, preset parsing,
action provenance, structural policy, and the Rust 1.78 combined-feature check
pass. Status remains `in_review` until the exact hosted matrix runs; no local
full suite was started.
