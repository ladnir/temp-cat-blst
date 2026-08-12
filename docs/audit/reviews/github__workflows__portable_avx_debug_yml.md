# File Review: `.github/workflows/portable-avx-debug.yml`

## Record

| Field | Value |
| --- | --- |
| Path | `.github/workflows/portable-avx-debug.yml` |
| Phase | Phase 19, Packet 19B |
| Reviewed commit | `0c7865082575e15f779debab8186bd0eca0eea7a` |
| Reviewer | Codex audit team |
| Date | 2026-08-12 |
| Disposition | `in_scope_support` |
| Subsystem | Developer AVX emulation CI |
| Profile | portable AVX Debug only |

## Purpose and trust boundary

This path-scoped developer workflow exercises the portable AVX emulation and
its differential/audit tests when relevant sources change. It is not a
production deployment profile and must never be presented as native AVX
release evidence.

## Contracts and invariants

- Debug only; portable-AVX Release is prohibited by compile-cost policy.
- Native AVX, assembly, ADX, CUDA, and benchmarks are disabled.
- The emulation implementation is selected as the default vector provider.
- Relevant AVX, field, curve, FFT, C-boundary, and test-owner paths trigger it.
- Two workers and a 45-minute timeout bound the expensive compile.

## Checklist evidence

| # | Audit question | Result |
| --- | --- | --- |
| 1 | Trigger scope | Path-scoped PR/master push plus manual dispatch. |
| 2 | Permissions | Read-only token; checkout credentials are not persisted. |
| 3 | Concurrency | Superseded work for the same ref is cancelled. |
| 4 | Runner | Fixed Ubuntu 24.04. |
| 5 | Build profile | CMake Debug, tests on, benchmarks off. |
| 6 | Provider isolation | Assembly/default assembly/ADX/native AVX all off. |
| 7 | Emulation selection | Portable AVX and default AVX on. |
| 8 | CUDA/parallel scope | Both off. |
| 9 | Product claim | Developer correctness only; no deployment support claim. |
| 10 | Resource behavior | Build/test use at most two workers and 45 minutes. |
| 11 | Test behavior | Builds `blst_tests` and runs registered tests with failures visible. |
| 12 | Action integrity | Checkout is an immutable verified commit. |
| 13 | Structural gate | CI policy requires Debug and the exact provider settings. |
| 14 | Residual risk | Exact hosted compile/test remains pending. |

## Findings and validation

`AUD-347`, `AUD-349`, and `AUD-350` are validated at the policy level. Status
remains `in_review` until this exact path-scoped job executes remotely. No
portable-AVX Release build was started.

