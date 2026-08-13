# File Review: `blst/avx/README.md`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/avx/README.md` |
| Roadmap phase / queue | Phase 16 AVX documentation closure |
| Reviewed commit | `5b08edc8f890ef1cbf4c6443162ec4da985879ee` |
| Reviewer | Codex audit team |
| Date | 2026-08-13 |
| Disposition | `in_scope_production` |
| Subsystem | AVX-512 support documentation |
| Profiles | native AVX-512 and portable developer Debug |
| Generation/provenance | handwritten Category Labs documentation |
| Upstream references | VROOM manifest and port map linked from the file |
| Related files | root CMake, AVX source tree, proof ledger, VROOM manifest |

## Purpose and trust boundary

The README tells deployers which backend is production, which CPU and compiler
features it requires, how compile-time selection works, and which proof claims
do and do not cover the implementation. It is not an API tutorial.

## Contracts and invariants

Native AVX requires x86-64 with AVX2, AVX-512F, AVX-512DQ, and AVX-512IFMA.
Selection is compile-time, with no scalar/AVX runtime dispatch. Portable AVX is
developer-only Debug/differential infrastructure and is not deployment or
performance evidence.

## Dependency and call-path map

The layout points to the selected vector, RNS, field, curve, pairing, GT, FFT,
generated-constant, validation, provenance, and formal-evidence owners.

## Upstream/specification comparison

The old import-era text and empty precompute placeholder were obsolete. The
replacement identifies VROOM precisely and links the retained manifest,
adaptation map, proof plan, and exact correspondence ledger.

## Checklist evidence

Sections 1-3, 9, 12, 13, and 14 apply and are covered by the support policy,
layout, provenance, proof-scope, and validation links. Sections 4-8 and 10-11
are N/A because the file contains no executable interface or state.

## Findings

No separate finding. This was documentation remediation discovered during the
four-file Phase 16 review. The obsolete precompute README was deleted.

## Validation

All 15 local links resolve. The documented build policy agrees with audited
CMake and CI records. `git diff --check` passes.

## Automated evidence

N/A for source analyzers and sanitizers. Repository link and structural checks
are the relevant evidence.

## Residual risk and decision

Future backend, compiler, ISA, or proof-boundary changes require a delta review.
The current file is complete at the reviewed commit.
