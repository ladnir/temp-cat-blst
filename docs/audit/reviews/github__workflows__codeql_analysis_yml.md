# File Review: `.github/workflows/codeql-analysis.yml`

## Record

| Field | Value |
| --- | --- |
| Path | `.github/workflows/codeql-analysis.yml` |
| Phase | Phase 19, Packet 19B; standing Phase 12 analysis stream |
| Reviewed commit | `0c7865082575e15f779debab8186bd0eca0eea7a` |
| Reviewer | Codex audit team |
| Date | 2026-08-12 |
| Disposition | `in_scope_support` |
| Subsystem | Static analysis CI |
| Profiles | retained C/C++ scalar, portable-AVX Debug, Rust |

## Purpose and trust boundary

This manual workflow extracts retained C/C++ and Rust production sources for
GitHub CodeQL. SARIF upload crosses the GitHub Code Security service boundary
and is therefore disabled unless the repository service is enabled explicitly.

## Contracts and invariants

- Manual dispatch only; `CODEQL_ENABLED=true` acknowledges service readiness.
- C/C++ uses a manual build of scalar and portable-AVX Debug archives.
- Rust uses CodeQL's supported no-build extraction.
- CUDA, benchmarks, frozen BTx/BTE targets, and removed bindings are excluded.
- Portable AVX Release is forbidden by compile-cost policy.

## Checklist evidence

| # | Audit question | Result |
| --- | --- | --- |
| 1 | Trigger scope | Manual and service-gated. |
| 2 | Permissions | Read-only globally; SARIF job adds only `security-events: write`. |
| 3 | Action integrity | CodeQL init/analyze and checkout are immutable verified commits. |
| 4 | Runner/resource | Ubuntu 24.04, ref cancellation, 45-minute timeout, two workers. |
| 5 | C/C++ extraction | Builds only `blst_static` in two explicit Debug profiles. |
| 6 | Rust extraction | No-build mode covers tracked Rust source. |
| 7 | Experimental scope | CUDA and parallel MSM disabled. |
| 8 | Frozen scope | No BTx/BTE target is built. |
| 9 | Performance | Debug emulation retained; portable Release excluded. |
| 10 | Failure behavior | Build or analysis failure fails the matrix row. |
| 11 | Output | SARIF only through CodeQL's security-events boundary. |
| 12 | Prior evidence | Hosted extraction reached SARIF; local CodeQL baselines are recorded. |
| 13 | Known blind spots | AVX template extraction diagnostics remain covered by Clang analyzer/proofs. |
| 14 | Residual risk | Exact current hosted run waits for Code Security enablement. |

## Findings and validation

No new CodeQL design defect was found. `AUD-349` and `AUD-350` apply to the
shared workflow policy and runner pin. Status remains `in_review` because the
service-gated exact revision has not uploaded SARIF.

