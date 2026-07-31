# BTx/BTE Freeze Policy

## State

BTx/BTE is frozen pending a required future merge. It is not a removal
candidate, and the current snapshot is not ready for extraction or final audit.

## Rules while frozen

- Do not modify BTx/BTE source, headers, tests, build files, or documentation.
- Do not remove, extract, rename, or otherwise reorganize the subsystem.
- Do not begin its file-by-file audit or use the current snapshot for final
  production-readiness conclusions.
- Do not run BTx/BTE benchmarks. They remain inventoried under the repository's
  general benchmark exclusion.
- Repository-wide correctness suites may continue to execute existing BTx/BTE
  tests as passive regression coverage, provided the subsystem is not changed.
- Cross-cutting work must treat BTx/BTE files as no-touch inputs and must not
  silently alter their generated or public surfaces.

## Unfreeze prerequisite

The required BTx/BTE merge must be identified and completed first. After that
merge, re-inventory the subsystem, record its new baseline, and obtain an
explicit decision before auditing, extracting, or removing it.

Until then, manifest rows use disposition `out_of_scope_frozen`, status
`frozen`, subsystem `btx_bte`, and profile `frozen_pending_merge`.
