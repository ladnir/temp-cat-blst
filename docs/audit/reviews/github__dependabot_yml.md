# File Review: `.github/dependabot.yml`

## Record

| Field | Value |
| --- | --- |
| Path | `.github/dependabot.yml` |
| Phase | Phase 19, Packet 19B |
| Reviewed commit | `0c7865082575e15f779debab8186bd0eca0eea7a` |
| Reviewer | Codex audit team |
| Date | 2026-08-12 |
| Disposition | `in_scope_support` |
| Subsystem | CI dependency automation |
| Provenance | Local repository policy |

## Purpose and trust boundary

This file asks GitHub's Dependabot service to propose updates to external
GitHub Actions. It must preserve reviewable immutable references rather than
turning an update service into an implicit mutable dependency channel.

## Contracts and invariants

- Only the `github-actions` ecosystem is automated in this packet.
- The repository root is scanned weekly.
- Related action changes are grouped into one reviewable update.
- Workflows themselves require full 40-hex action commit SHAs.
- Dependabot proposals never bypass review or acquire workflow write tokens.

## Checklist evidence

| # | Audit question | Result |
| --- | --- | --- |
| 1 | Purpose and ownership | Weekly action-maintenance policy only. |
| 2 | Inputs and trust | GitHub service metadata and upstream action releases are external inputs. |
| 3 | Output/publication | Produces pull requests; it does not merge or publish. |
| 4 | Failure behavior | Service failure leaves existing immutable pins unchanged. |
| 5 | Security/permissions | No secrets, token grants, or executable steps occur here. |
| 6 | Dependency integrity | Runtime workflows enforce full SHAs independently. |
| 7 | Resource use | One grouped weekly stream bounds update churn. |
| 8 | Platform behavior | GitHub-hosted declarative service; no host-specific script. |
| 9 | Compatibility | Targets the repository default branch and standard Actions ecosystem. |
| 10 | Performance | Not applicable to product code or hot paths. |
| 11 | Concurrency | GitHub manages update scheduling; proposals remain review-gated. |
| 12 | Documentation | Policy is recorded in Phase 12 and Phase 19. |
| 13 | Tests | `check_ci_policy.ps1` validates ecosystem, root, weekly schedule, and group. |
| 14 | Residual risk | Cargo version automation waits for the `AUD-351` MSRV decision. |

## Findings and validation

No defect was found in the retained action-update policy. All current action
SHAs resolve to the human-readable upstream tags recorded in the workflows.
The file is complete at the reviewed commit.

