# File Review: `.github/workflows/rust-bindings.yml`

## Record

| Field | Value |
| --- | --- |
| Path | `.github/workflows/rust-bindings.yml` |
| Phase | Phase 18 regeneration owner; Phase 19 Packet 19B policy refresh |
| Reviewed commit | `0c7865082575e15f779debab8186bd0eca0eea7a` |
| Reviewer | Codex audit team |
| Date | 2026-08-12 |
| Disposition | `in_scope_support` |
| Subsystem | Rust FFI generation CI |
| Profile | scalar raw binding regeneration only |

## Purpose and trust boundary

This workflow proves that the committed raw Rust binding is exactly generated
from the reviewed C99 scalar header and deterministic transformation policy.
It installs external compiler/generator tools, so versions and output handling
must be explicit and the checkout must remain unmodified.

## Contracts and invariants

- Runs only for binding inputs/outputs, generator policy, and its own workflow.
- Uses Clang/libclang 18, bindgen CLI 0.65.1 with its lock, and Rust 1.65 rustfmt.
- Generates into a temporary file and byte-compares against the tracked output.
- Performs no native product build, publication, or repository write.
- Read-only permissions, cancellation, and a 30-minute timeout apply.

## Checklist evidence

| # | Audit question | Result |
| --- | --- | --- |
| 1 | Trigger scope | Exact regeneration-owner path list plus manual dispatch. |
| 2 | Permissions | Read-only token; checkout credentials disabled. |
| 3 | Action integrity | Immutable verified checkout commit. |
| 4 | Runner | Fixed Ubuntu 24.04. |
| 5 | Toolchain | Clang/libclang major, bindgen version/lock, rustfmt version pinned. |
| 6 | Output ownership | Temporary sibling output; tracked file is comparison-only. |
| 7 | Cleanup | Shell trap removes the temporary generated file. |
| 8 | Failure behavior | Any install, generation, comparison, or policy failure stops. |
| 9 | Transform validation | Bash syntax, Perl syntax, fixture, and PowerShell parity checks. |
| 10 | Product isolation | No CMake/Cargo product compile or benchmark. |
| 11 | Resource bounds | Ref cancellation and 30-minute timeout. |
| 12 | Reproducibility | Two local generation locations previously produced identical bytes. |
| 13 | Platform scope | Linux generator owner only; generated output is platform-neutral ABI text. |
| 14 | Residual risk | Hosted libclang minor/package repository drift remains visible evidence. |

## Findings and validation

`AUD-324` remains validated. Phase 19 adds only shared resource and permission
policy; actionlint and the exact local generator equivalent pass. The file
review remains complete, while first hosted execution of this exact revision is
phase-level evidence rather than an unreviewed file contract.

