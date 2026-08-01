# Risk-Proportional Validation Policy

The audit optimizes for short implementation cycles without weakening the
release gate. Validation depth follows the change's risk and affected surface;
the full matrix is not a per-commit requirement.

## Automated baseline and hand-audit sequencing

Establish the initial platform milestone in `AUTOMATED_ANALYSIS_PLAN.md` before
advancing Phase 13. Triage its results first so the hand review can prioritize
concrete suspicious paths. Automated analysis then remains active through the
subsystem phases in `PRODUCTION_AUDIT_ROADMAP.md`; it is not a one-time prelude.
Do not rerun repository-wide CodeQL, Clang Static Analyzer, sanitizer, or fuzz
sweeps for every individual file.

During file review, map baseline results to the file and run only focused
checks that exercise its contracts. Expand fuzz and sanitizer coverage as each
subsystem becomes understood. Rerun the complete applicable automated stack at
the subsystem gate and the full supported stack at milestone/release gates.
Static analysis guides attention but never replaces manual correctness,
constant-time, arithmetic, assembly, or performance review.

## Tier 1 - change loop

Run these checks after each narrow implementation change:

- build only the directly affected target or profile;
- run the new regression and the closest existing tests;
- run formatting, whitespace, or audit-integrity checks that cover the edited
  files; and
- record what was not run and why.

Use the following surface-specific checks:

| Change surface | Minimum focused evidence |
| --- | --- |
| Public C/C++ ABI or header | Standalone C and C++ header compile, focused API regression, and symbol/link check |
| Rust API/FFI | Focused Rust test, formatting, and one native link/build path |
| Scalar arithmetic or protocol logic | Focused oracle/differential regression in each directly affected scalar backend |
| AVX-512 source or dispatch | Rebuild the AVX-512 target and run the closest AVX/differential regression |
| Assembly generator or generated output | Deterministic regeneration/diff plus focused build and symbol/behavior checks for the affected ISA |
| Memory ownership, bounds, or unsafe parsing | Boundary regression plus a focused sanitizer or Miri run when it materially exercises the risk |
| Analyzer or fuzz finding | Finding entry, minimized reproducer where possible, focused regression, and rerun of the originating rule/target |
| Build, packaging, or CI | Parse/configure/package the changed path and build one representative consumer |
| Documentation or audit metadata only | Audit integrity and link/reference checks; no product rebuild |

A failure expands validation one dependency ring at a time. It does not
automatically trigger every supported profile.

## Tier 2 - subsystem gate

After a coherent implementation packet, or several related leaf commits, run
the complete suite for the affected subsystem and its supported consumer
surface once. Examples are a scalar-field packet, the Rust API packet, or an
AVX-512 dispatch packet. Batch nearby changes before paying this cost.

## Tier 3 - milestone and release gates

Run the full supported matrix only at planned checkpoints: completion of a
major audit phase, closure of a production subsystem, rename/ABI integration,
and release candidates. The gate may include direct scalar, portable dispatch,
AVX-512, Rust, sanitizers, cross-platform builds, packaging, and CodeQL as
appropriate to that milestone.

Tool versions, rules/query suites, compilation profiles, fuzz budgets, corpus
identity, exclusions, and false-positive decisions are evidence and must be
recorded. A clean automated run is not evidence of cryptographic correctness,
constant-time behavior, or assembly safety.

## Persistent exclusions

- Never run two benchmark processes concurrently.
- Benchmarks are excluded from the file-by-file audit and run only for an
  explicit performance question or a planned performance gate.
- Frozen BTx/BTE is excluded from new audit, sanitizer, and fuzz work until its
  required merge is ready. Existing repository-wide builds and tests may keep
  compiling and executing it as passive regression coverage, provided the
  frozen subsystem is not changed or used for production-readiness claims.
- Experimental CUDA is outside the primary scalar/AVX-512 gate unless a change
  directly affects CUDA or an explicit experimental checkpoint is scheduled.
- Experimental non-Apple ARM64 and CHERI are blocked by default and remain
  outside the production matrix until their target-specific promotion gates
  are complete. Apple Silicon macOS is supported scalar scope and requires its
  own release evidence.
