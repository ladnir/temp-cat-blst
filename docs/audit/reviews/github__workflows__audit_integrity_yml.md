# File Review: `.github/workflows/audit-integrity.yml`

## Record

| Field | Value |
| --- | --- |
| Path | `.github/workflows/audit-integrity.yml` |
| Phase | Phase 19, Packet 19B |
| Reviewed commit | `0c7865082575e15f779debab8186bd0eca0eea7a` |
| Reviewer | Codex audit team |
| Date | 2026-08-12 |
| Disposition | `in_scope_support` |
| Subsystem | Audit-policy CI |
| Related files | `scripts/audit/check_audit_manifest.ps1`, `scripts/audit/check_ci_policy.ps1`, Rust policy scripts |

## Purpose and trust boundary

This lightweight workflow rejects drift in audit inventories, retained Rust
declarations/contracts, binding transforms, and CI policy. Pull-request source
is untrusted, so the workflow must run without credentials or write access.

## Contracts and invariants

- Runs for every pull request, master push, and manual request.
- Superseded work for the same ref is cancelled.
- The token is read-only and checkout persistence is disabled.
- Windows owns PowerShell policy checks; Ubuntu owns Bash/Perl transforms.
- Every job has a 15-minute timeout and performs no product build.

## Checklist evidence

| # | Audit question | Result |
| --- | --- | --- |
| 1 | Trigger scope | PR, master push, and manual triggers are intentional. |
| 2 | Untrusted inputs | Checked-out PR files are read/executed without secrets or write token. |
| 3 | Permissions | Workflow-level `contents: read`; persisted credentials disabled. |
| 4 | Action integrity | Checkout uses the verified immutable v7.0.1 commit. |
| 5 | Resource bounds | Ref concurrency and 15-minute per-job timeouts. |
| 6 | Platform pinning | Windows 2025 and Ubuntu 24.04. |
| 7 | Shell semantics | Explicit PowerShell and Bash shells; paths are repository constants. |
| 8 | Failure behavior | Every policy script and transform fails the job directly. |
| 9 | Output/publication | No artifact upload, mutation, cache, or publication. |
| 10 | Product isolation | No benchmark, CUDA, product build, or test executable. |
| 11 | Frozen scope | BTx/BTE is inventoried but not audited or modified. |
| 12 | Performance | Small text/transform checks only. |
| 13 | Validation | actionlint 1.7.12 and every command's local equivalent pass. |
| 14 | Residual risk | Future workflow changes require a delta review and exact hosted execution. |

## Findings and validation

`AUD-349` and `AUD-350` are validated. Commit `5d43ad1` adds the CI policy
gate; `0c78650` tightens its runner and Dependabot invariants. The workflow
blob is unchanged in commit `3b21858`, where hosted run `31665323679` passes
both jobs: manifests, retained CI policy, raw Rust binding integrity, reviewed
Rust contracts, and the binding transform. The file is complete.
