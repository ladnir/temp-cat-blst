# Risk-Proportional Validation Policy

The audit optimizes for short implementation cycles without weakening the
release gate. Validation depth follows the change's risk and affected surface;
the full matrix is not a per-commit requirement.

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

## Persistent exclusions

- Never run two benchmark processes concurrently.
- Benchmarks are excluded from the file-by-file audit and run only for an
  explicit performance question or a planned performance gate.
- Frozen BTx/BTE is not built or tested until its required merge is ready.
- Experimental CUDA is outside the primary scalar/AVX-512 gate unless a change
  directly affects CUDA or an explicit experimental checkpoint is scheduled.
- Experimental non-Apple ARM64 and CHERI are blocked by default and remain
  outside the production matrix until their target-specific promotion gates
  are complete. Apple Silicon macOS is supported scalar scope and requires its
  own release evidence.
