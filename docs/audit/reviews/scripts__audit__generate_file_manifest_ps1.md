# File Review: `scripts/audit/generate_file_manifest.ps1`

## Record

| Field | Value |
| --- | --- |
| Path | `scripts/audit/generate_file_manifest.ps1` |
| Roadmap phase / queue | Continuous audit-control integrity |
| Reviewed commit | `a9ceb3623cadc794cdf13559db22dde7279ea35a` |
| Reviewer | Codex audit team |
| Date | 2026-08-13 |
| Disposition | `in_scope_support` |
| Subsystem | Authoritative file-manifest regeneration |
| Profiles | PowerShell audit process |
| Generation/provenance | Locally maintained Category-blst audit infrastructure |
| Upstream references | `not_applicable: local audit-control owner` |
| Related files | `docs/audit/FILE_MANIFEST.csv`, `scripts/audit/check_audit_manifest.ps1`, `docs/audit/README.md` |

## Purpose and trust boundary

The script inventories every tracked or unignored repository path and adds a
manifest row for a new path. It must remove deleted paths without changing an
existing human audit decision. The manifest controls production scope and
release-readiness metrics, so silent reclassification is a security-process
failure even when product code is unchanged.

## Contracts and invariants

- Paths are repository-relative, slash-normalized, unique, and sorted.
- New paths receive a deterministic default classification.
- Every nonempty field of an existing row is authoritative and survives
  regeneration, including disposition, subsystem, status, profiles, and notes.
- Deleted paths disappear. The manifest's own path is always present.
- Output is UTF-8 without a byte-order mark and is deterministic on a stable
  repository state.
- The script performs no product build, network operation, or external write.

## Dependency and call-path map

Auditors invoke this script before integrity validation. It reads Git's tracked
and unignored path set plus the existing CSV, applies defaults only to new
paths, and rewrites the CSV. `check_audit_manifest.ps1` then compares the CSV
with Git and validates every disposition/status pair and completed record.

## Checklist evidence

- The path set uses ordinal comparison and explicit slash normalization.
- CSV parsing and emission use PowerShell's structured CSV operations.
- Existing nonempty human fields now override classifier defaults uniformly.
- A regression reproduction retained all 36 experimental ARM/CHERI rows,
  support-file overrides, completed statuses, and review evidence.
- Two consecutive regeneration runs produced identical SHA-256 values.
- The manifest integrity checker validates the resulting exact path set.
- N/A for runtime memory safety, cryptographic timing, randomness, API ABI,
  concurrency, and performance: this is a bounded serial audit tool.

## Findings

`AUD-383`.

## Validation

- Reproduced the prior reclassification of 36 experimental rows and four
  reviewed support overrides.
- Restored the authoritative rows and reran the corrected generator twice.
- Second regeneration was byte-identical.
- Manifest integrity and exact status counts pass after the integration packet.

## Residual risk and decision

Classification defaults can still be imperfect for a new exceptional path;
the generated row remains pending and requires human review. Existing human
decisions are no longer overwritten. The owner is complete for file review.
