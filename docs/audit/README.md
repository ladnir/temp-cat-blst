# CatBLST Production Audit

This directory is the durable record for the production-readiness audit of the
BLST-derived codebase that is expected to become CatBLST.

## Objectives

1. Account for every upstream commit after the recorded merge base and decide
   whether to port, translate, reject, defer, or mark it already equivalent.
2. Establish and triage the initial automated-analysis platform before
   advancing the first main audit subsystem, then expand that workstream
   continuously through release qualification.
3. Deep-review every production and production-support file, one file at a
   time, with explicit manual and applicable automated evidence.
4. Harden and validate the scalar and AVX-512 implementations.
5. Preserve the generic FFT and middle-product implementations; keep retained
   interpolation isolated behind a minimal explicit surface.
6. Keep CUDA, non-Apple ARM64, and CHERI retained and explicitly experimental;
   keep Apple Silicon macOS in supported scalar scope. Unsupported language
   bindings were removed in Phases 6B and 6C. Keep BTx/BTE frozen pending its
   required merge.
7. Rename the project and introduce a configurable C ABI prefix only after the
   upstream and code audits are under control.

The current name and ABI remain unchanged during the audit. A rename is not an
audit shortcut, and removal candidates remain inventoried until their deletion
is reviewed and committed.

C++ support means the native C++ implementation and API. The legacy
`bindings/blst.hpp` convenience layer over the C ABI was removed in Phase 6A;
it is not part of the supported C++ surface.

The project is pre-production and has no external users. Breaking source or ABI
changes are acceptable when they materially improve safety, correctness,
performance, or API clarity. Do not retain unsafe legacy entry points solely
for compatibility.

## Source of truth

| Artifact | Purpose |
| --- | --- |
| `UPSTREAM_MANIFEST.md` | One decision row for every upstream-only commit |
| `UPSTREAM_BULK_APPROVAL.md` | Approved translations and validation conditions for the final upstream implementation packets |
| `FILE_MANIFEST.csv` | Machine-checkable row for every repository file |
| `CHECKLIST.md` | Required questions and evidence for each file |
| `FILE_REVIEW_TEMPLATE.md` | Template for `reviews/*.md` review records |
| `FINDINGS.md` | Security, correctness, performance, and maintainability findings |
| `STATUS.md` | Current phase, gates, metrics, and next review queue |
| `HARDENING_PLAN.md` | Dependency-aware implementation and validation phases |
| `PRODUCTION_AUDIT_ROADMAP.md` | Main Phases 12-20, continuous workstreams, subsystem cycles, work volume, and completion gates |
| `../avx/README.md` | Phase 16 AVX-512 modernization, MatrixNoK import boundary, validation profiles, and completion gates |
| `../avx/VROOM_UPSTREAM_MANIFEST.md` | Decision for every VROOM commit after the imported baseline and the selected BLS12-381 translation groups |
| `../avx/PROOF_PLAN.md` | Lean proof layers, C++ correspondence boundary, constant-time evidence, and first theorem milestone |
| `AUTOMATED_ANALYSIS_PLAN.md` | Static-analysis, sanitizer, fuzzing, triage, cadence, and hand-audit integration policy |
| `FUZZING_PLAN.md` | Phased API-complete fuzzing design, domain enforcement, resource policy, and acceptance gates |
| `FUZZ_API_LEDGER.csv` | Symbol-level campaign, domain, factory, invalid-shape, and oracle disposition for supported APIs |
| `CODEQL_LOCAL.md` | Pinned local CodeQL installation, execution, evidence, and license boundary |
| `CODEQL_BASELINE_TRIAGE.md` | Result-by-result baseline dispositions, confirmed findings, and proposed implementation packets |
| `CLANG_ANALYZER_BASELINE_TRIAGE.md` | Split scalar/AVX Clang analyzer evidence, dispositions, and proposed cleanup packet |
| `CLANG_TIDY_BASELINE_TRIAGE.md` | Curated scalar/portable-AVX Debug clang-tidy evidence and dispositions |
| `CLIPPY_BASELINE_TRIAGE.md` | Rust feature-family Clippy evidence, dispositions, and proposed scratch-storage finding |
| `SANITIZER_BASELINE_TRIAGE.md` | Focused scalar ASan/UBSan configuration, evidence, exclusions, and configuration finding |
| `GROUP_API_CONTRACT_AUDIT.md` | Cross-language point/Fp12/GT/pairing domain, validation, timing, and trust-boundary audit |
| `GROUP_API_CONTRACTS.csv` | Machine-checked contract-family inventory for every group-domain public C symbol |
| `PUBLIC_API_DOCUMENTATION.csv` | Incremental C/C++/Rust declaration coverage for already audited API families |
| `../API_CONTRACTS.md` | Public vocabulary for precise preconditions, postconditions, validation, domains, timing, and failure behavior |
| `PHASE_13.md` | Foundational hand-audit findings, C ABI exception policy, remediation sequence, and focused evidence |
| `PHASE_17.md` | FFT native-core redesign, adapter boundary, packaging, and focused validation evidence |
| `VALIDATION_POLICY.md` | Risk-proportional checks and subsystem/milestone validation gates |
| `BTX_FREEZE.md` | No-touch policy and merge prerequisite for frozen BTx/BTE |
| `CUDA_EXPERIMENTAL.md` | Retention and non-production policy for opt-in CUDA |
| `ARM64_CHERI_EXPERIMENTAL.md` | Apple support boundary plus default-deny policy and qualification gates for other ARM64 and CHERI targets |
| `PHASE_*.md` | Implemented hardening and validation evidence |
| `scripts/audit/generate_file_manifest.ps1` | Regenerates the inventory without erasing review evidence |
| `scripts/audit/check_audit_manifest.ps1` | Enforces file/upstream inventory integrity and completion evidence |
| `scripts/audit/check_group_api_contracts.ps1` | Enforces complete and unambiguous public group-domain contract coverage |
| `scripts/audit/check_public_api_documentation.ps1` | Enforces source, anchor, scope, and disposition coverage for audited public API documentation |
| `scripts/audit/check_rust_threading_policy.ps1` | Enforces the caller-threaded Rust runtime contract |
| `scripts/audit/run_sanitizers_local.sh` | Runs the pinned focused scalar ASan/UBSan boundary baseline on x86-64 Linux |
| `scripts/fuzz/api_inventory.py` | Generates and validates the exact C, C++, raw Rust, and safe Rust fuzz ledger |
| `scripts/fuzz/check_structure.py` | Enforces the private adapter boundary and rejects unclassified direct calls |
| `scripts/fuzz/generate_api_inventory.sh` | Regenerates compact Clang AST and rustdoc snapshots with bounded workers |

## Scope dispositions

- `in_scope_production`: shipped implementation, ABI, or supported binding.
- `in_scope_support`: tests, benchmarks, build logic, CI, vectors, audit
  controls, and documentation needed to ship safely.
- `generated_production`: generated assembly, constants, or bindings that ship
  or determine shipped code. Generated does not mean unaudited.
- `out_of_scope_frozen`: inventoried code under an explicit no-touch policy
  until a recorded prerequisite is completed.
- `out_of_scope_experimental`: retained opt-in code outside the primary
  production-readiness profiles and guarantees.
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
- `frozen`: no modification or audit is allowed until the recorded prerequisite
  is completed and the subsystem is explicitly unfrozen.
- `experimental`: retained and classified, but outside the primary production
  readiness claim and file-by-file completion gate.
- `scope_review_required`: planned-removal classification still needs approval.
- `excluded`: explicitly outside the file-by-file audit.
- `catalogued`: inventoried non-shipping documentation or archive.

## Operating procedure

1. Regenerate `FILE_MANIFEST.csv`.
2. Run the integrity checker before starting a review.
3. Establish and triage the initial platform milestone in
   `AUTOMATED_ANALYSIS_PLAN.md` before Phase 13, then keep its automated
   workstream active through the subsystem phases in
   `PRODUCTION_AUDIT_ROADMAP.md`.
4. Select the next file in dependency order from `STATUS.md`.
5. Copy `FILE_REVIEW_TEMPLATE.md` to `reviews/<stable-name>.md`.
6. Review the file at an exact commit using `CHECKLIST.md`.
7. Record findings in `FINDINGS.md`; do not hide them in the review narrative.
8. Add tests or other validation proportional to the risk.
9. Update the file-manifest row. A file may be marked `complete` only when every
   required column is populated.
10. Run the integrity checker again.

Use `VALIDATION_POLICY.md` to select the smallest check set that directly
exercises each change. Full supported-profile matrices are subsystem,
milestone, and release gates rather than per-commit requirements.

Use one benchmark process at a time. Benchmark evidence must include the build
profile, CPU, instruction-set assumptions, command, sample count, and raw result
location.

## Commands

From the repository root:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File scripts/audit/generate_file_manifest.ps1
powershell.exe -NoProfile -ExecutionPolicy Bypass -File scripts/audit/check_audit_manifest.ps1
powershell.exe -NoProfile -ExecutionPolicy Bypass -File scripts/audit/check_group_api_contracts.ps1
powershell.exe -NoProfile -ExecutionPolicy Bypass -File scripts/audit/check_rust_threading_policy.ps1
powershell.exe -NoProfile -ExecutionPolicy Bypass -File scripts/audit/check_win64_asm_symbols.ps1 -BuildDirectory out/build/<profile> -Mode Portable
wsl.exe bash scripts/audit/run_clang_tidy_local.sh all
```

The generator preserves human-entered review fields by matching the `path`
column. Renames must be reconciled explicitly so that review evidence is not
silently attributed to a different file.

Use `-Mode NonPortable` for a build configured with `BLST_PORTABLE=OFF`. The
Win64 symbol checker requires MSVC `dumpbin.exe` on `PATH`, or an explicit
`-Dumpbin` path.
