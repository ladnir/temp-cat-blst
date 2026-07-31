# CatBLST Production Audit

This directory is the durable record for the production-readiness audit of the
BLST-derived codebase that is expected to become CatBLST.

## Objectives

1. Account for every upstream commit after the recorded merge base and decide
   whether to port, translate, reject, defer, or mark it already equivalent.
2. Deep-review every production and production-support file, one file at a
   time, with explicit evidence.
3. Harden and validate the scalar and AVX-512 implementations.
4. Preserve the FFT and middle-product implementations.
5. Prepare planned removal of CUDA, Go, BTx/BTE, and unsupported language
   bindings.
6. Rename the project and introduce a configurable C ABI prefix only after the
   upstream and code audits are under control.

The current name and ABI remain unchanged during the audit. A rename is not an
audit shortcut, and removal candidates remain inventoried until their deletion
is reviewed and committed.

C++ support means the native C++ implementation and API. The legacy
`bindings/blst.hpp` convenience layer over the C ABI is not supported and is a
planned-removal surface.

The project is pre-production and has no external users. Breaking source or ABI
changes are acceptable when they materially improve safety, correctness,
performance, or API clarity. Do not retain unsafe legacy entry points solely
for compatibility.

## Source of truth

| Artifact | Purpose |
| --- | --- |
| `UPSTREAM_MANIFEST.md` | One decision row for every upstream-only commit |
| `FILE_MANIFEST.csv` | Machine-checkable row for every repository file |
| `CHECKLIST.md` | Required questions and evidence for each file |
| `FILE_REVIEW_TEMPLATE.md` | Template for `reviews/*.md` review records |
| `FINDINGS.md` | Security, correctness, performance, and maintainability findings |
| `STATUS.md` | Current phase, gates, metrics, and next review queue |
| `HARDENING_PLAN.md` | Dependency-aware implementation and validation phases |
| `PHASE_*.md` | Implemented hardening and validation evidence |
| `scripts/audit/generate_file_manifest.ps1` | Regenerates the inventory without erasing review evidence |
| `scripts/audit/check_audit_manifest.ps1` | Enforces file/upstream inventory integrity and completion evidence |

## Scope dispositions

- `in_scope_production`: shipped implementation, ABI, or supported binding.
- `in_scope_support`: tests, benchmarks, build logic, CI, vectors, audit
  controls, and documentation needed to ship safely.
- `generated_production`: generated assembly, constants, or bindings that ship
  or determine shipped code. Generated does not mean unaudited.
- `out_of_scope_planned_removal`: retained temporarily, but intended to be
  removed. The removal itself requires review.
- `out_of_scope_no_audit`: inventoried non-production material that is
  explicitly excluded from file-by-file code audit.
- `documentation_or_archive`: non-shipping historical or explanatory material.

## File review states

- `pending`: not yet reviewed.
- `in_review`: an owner is actively reviewing it.
- `blocked`: review cannot finish until a recorded dependency is resolved.
- `complete`: all required evidence exists and the integrity checker accepts it.
- `accepted_risk`: review is complete but residual risk has been explicitly
  accepted and linked.
- `scope_review_required`: planned-removal classification still needs approval.
- `excluded`: explicitly outside the file-by-file audit.
- `catalogued`: inventoried non-shipping documentation or archive.

## Operating procedure

1. Regenerate `FILE_MANIFEST.csv`.
2. Run the integrity checker before starting a review.
3. Select the next file in dependency order from `STATUS.md`.
4. Copy `FILE_REVIEW_TEMPLATE.md` to `reviews/<stable-name>.md`.
5. Review the file at an exact commit using `CHECKLIST.md`.
6. Record findings in `FINDINGS.md`; do not hide them in the review narrative.
7. Add tests or other validation proportional to the risk.
8. Update the file-manifest row. A file may be marked `complete` only when every
   required column is populated.
9. Run the integrity checker again.

Use one benchmark process at a time. Benchmark evidence must include the build
profile, CPU, instruction-set assumptions, command, sample count, and raw result
location.

## Commands

From the repository root:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File scripts/audit/generate_file_manifest.ps1
powershell.exe -NoProfile -ExecutionPolicy Bypass -File scripts/audit/check_audit_manifest.ps1
```

The generator preserves human-entered review fields by matching the `path`
column. Renames must be reconciled explicitly so that review evidence is not
silently attributed to a different file.
