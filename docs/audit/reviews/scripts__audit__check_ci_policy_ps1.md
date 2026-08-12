# File Review: `scripts/audit/check_ci_policy.ps1`

## Record

| Field | Value |
| --- | --- |
| Path | `scripts/audit/check_ci_policy.ps1` |
| Phase | Phase 19, Packet 19B |
| Reviewed commit | `a71068a47c361c9a1edae4a706b75435023756f6` |
| Reviewer | Codex audit team |
| Date | 2026-08-12 |
| Disposition | `in_scope_support` |
| Subsystem | Structural CI policy gate |
| Provenance | Locally written Apache-2.0/SPDX script |

## Purpose and trust boundary

The script gives routine audit CI a fast offline guard against configuration
drift that a generic YAML parser does not understand: mutable action or runner
references, excess token authority, runaway jobs, incorrect production
presets, unsafe Rust caching, unlocked dependencies, and Release portable AVX.

## Contracts and invariants

- Reads only repository workflow/Dependabot text relative to its own location.
- Performs no network access, process launch, source mutation, or product build.
- Fails closed on absent workflows, jobs, controls, or expected matrix entries.
- Permits local `./` actions but requires every external action to use 40 hex.
- Complements actionlint; it is not represented as a general YAML parser.

## Checklist evidence

| # | Audit question | Result |
| --- | --- | --- |
| 1 | Input ownership | Fixed `.github` files in the checked-out repository. |
| 2 | Path safety | Root derives from `PSScriptRoot`; all reads are literal paths. |
| 3 | Mutation | None. |
| 4 | Network/processes | None. |
| 5 | Error handling | Strict mode and stop-on-error; violations throw descriptive errors. |
| 6 | Workflow discovery | Requires at least one sorted `.yml` workflow. |
| 7 | Job bounds | Counts top-level jobs and requires one timeout each. |
| 8 | Supply chain | Rejects mutable external actions and latest runner aliases. |
| 9 | Permissions/resources | Requires workflow-level permissions and concurrency. |
| 10 | CMake policy | Requires the exact ordered four-preset production matrix. |
| 11 | Rust policy | Requires tracked v4 lock, Rust 1.78, and seven locked test/smoke/check commands. |
| 12 | AVX policy | Requires Debug emulation with assembly/ADX/native AVX disabled. |
| 13 | Dependabot policy | Requires weekly grouped GitHub Actions and Cargo updates. |
| 14 | Blind spot | Regex policy is paired with actionlint semantic YAML validation. |

## Findings and validation

The positive gate reports five workflows, eight immutable external action
uses, four supported presets, and seven locked Cargo commands. actionlint 1.7.12
independently parses all workflow YAML. The script review is complete; hosted
registration remains Phase 19 evidence.
