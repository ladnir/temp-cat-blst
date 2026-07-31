# Automated Analysis and Fuzzing Plan

## Purpose

Automated analysis is an evidence source for the file-by-file audit, not a
substitute for it. The bootstrap runs before the main hand audit so obvious
defects and high-risk files are identified early. Coverage then expands with
each subsystem after its contracts and independent oracles are understood.

## Scope

- Production profiles: scalar and AVX-512.
- Supported language surfaces: C, native C++, and Rust.
- Supported scalar platform requiring release evidence: Apple Silicon macOS.
- Excluded: benchmarks, frozen BTx/BTE, experimental CUDA, experimental
  non-Apple ARM64/CHERI, removed language bindings, and WASI/WebAssembly.
- No analysis result may justify adding runtime abstraction, allocation,
  dispatch, or generic loop structure to a hot path without a separate
  performance and code-generation review.

## Tool stack

| Tool | Primary purpose | Required scope | Important blind spots |
| --- | --- | --- | --- |
| CodeQL | Repository-wide security and interprocedural data-flow queries | Retained C/C++ and Rust with representative real builds | Handwritten/generated assembly, constant-time behavior, arithmetic proofs, scalar/AVX equivalence |
| Clang Static Analyzer through CodeChecker | Path-sensitive symbolic execution and C/C++ defect discovery | Portable scalar C/C++ compilation database first; expand to buildable shared code | Rust, assembly instructions, unsupported compiler-specific semantics |
| clang-tidy | Compiler-aware correctness, bug-prone, portability, and selected performance checks | Curated C/C++ rule set; generated files analyzed through their generators where appropriate | Whole-program behavior and cryptographic properties |
| Clippy | Rust correctness, suspicious, performance, and FFI-adjacent linting | Default, full bindings, portable, and AVX feature families where supported | Native C/C++ reached through FFI and runtime UB in foreign code |
| ASan/UBSan | Executed memory-safety and undefined-behavior detection | Linux Clang portable scalar baseline and focused subsystem targets | Assembly bodies, unexecuted paths, constant-time behavior |
| Coverage-guided fuzzing | Adversarial state and input exploration | Public byte boundaries first, then protocol and numerical subsystems | Properties not encoded as assertions/oracles and paths outside the corpus budget |

Compiler warnings remain mandatory input to this stack. Analyzer versions,
query/checker sets, build flags, CPU features, and action revisions must be
recorded exactly. GitHub Actions are pinned to immutable full commit SHAs with
the readable release tag in comments.

## Bootstrap sequence

1. Capture clean representative C/C++ and Rust builds and generate the C/C++
   compilation database.
2. Add CodeQL for C/C++ and Rust with manual builds that compile retained
   production sources rather than headers or empty targets.
3. Run Clang Static Analyzer/CodeChecker and the curated clang-tidy rules over
   the portable scalar build.
4. Run Clippy over supported Rust feature families without compiling removed
   or experimental surfaces into the release claim.
5. Run focused Linux Clang ASan/UBSan tests for existing public parsing,
   conversions, group operations, and hardened ownership/bounds paths.
6. Add initial bounded fuzz harnesses for byte conversion, G1/G2
   deserialize/uncompress/validate, and canonical round trips.
7. Triage every diagnostic and crash before the dependency-ordered hand audit
   advances beyond its affected file group.

## Fuzzing expansion order

Fuzz targets are added with the subsystem review that can define their
invariants:

1. Scalar, field, and point byte conversion: arbitrary lengths, modulus
   boundaries, non-canonical encodings, infinity, malformed flags, and
   serialize/deserialize round trips.
2. Hash-to-curve, key generation, signing, and verification: bounded messages,
   DSTs, augmentation, invalid keys/signatures, and cross-language agreement.
3. Pairing, aggregation, MSM, and uniqueness-tree APIs: empty/degenerate
   batches, checked size arithmetic, aliasing contracts, and bounded resource
   use.
4. AVX-512: differential properties against scalar for every exposed
   operation, including conversion boundaries, lazy-reduction limits, and
   inactive lanes. Portable AVX is useful evidence but cannot be reported as
   native AVX-512 execution.
5. FFT and middle product: sizes 0/1, invalid and boundary domains,
   in-place/overlap behavior, round trips, and comparison with simple bounded
   reference implementations.

Corpora begin with retained specification vectors, prior regressions, zero and
all-ones inputs, field moduli and adjacent values, exact buffer boundaries, and
minimized reproductions from every discovered failure. Harnesses must cap input
sizes and work factors so denial-of-service behavior is tested deliberately
rather than making CI nondeterministic.

## Finding and false-positive policy

- Every actionable diagnostic or reproducible failure receives an `AUD-*`
  entry in `FINDINGS.md` before code changes begin.
- Duplicate reports link to the controlling finding rather than being silently
  discarded.
- A false positive records the tool, rule/query, source location, exact
  configuration, technical dismissal, and reviewer. Broad rule suppression is
  not permitted when a narrow annotation or path exclusion is sufficient.
- Analyzer-driven changes receive the same regression, performance, and
  scalar/AVX review as manually discovered changes.
- Critical/high findings block dependent file completion. Lower findings may
  be sequenced or accepted only through the normal residual-risk process.
- A clean run means only that the configured tools emitted no untriaged result;
  it is never evidence of cryptographic correctness or constant-time behavior.

## Cadence

| Gate | Automated evidence |
| --- | --- |
| Pull request/change loop | Formatting and warnings; focused Clippy/tests; short deterministic sanitizer or fuzz regression only when directly relevant |
| Subsystem closure | Complete relevant CodeQL/Clang/Clippy rule sets, subsystem ASan/UBSan suite, bounded fuzz campaign, scalar/AVX differential tests |
| Scheduled CI | Longer fuzz budgets and analyzer configurations too costly for each change; corpus minimization and artifact retention |
| Milestone/release | Full supported static-analysis, sanitizer, fuzz-regression, platform, scalar, AVX-512, C/C++/Rust, and packaging matrix |

This cadence follows `VALIDATION_POLICY.md`: expensive clean sweeps are batched
at subsystem and milestone gates, not repeated after every small edit.

## Per-file integration

Each review record must identify:

- applicable baseline diagnostics and their disposition;
- exact static-analysis configuration and result, or `N/A` with a reason;
- sanitizer tests that execute the file's risky paths;
- fuzz targets, corpus seeds, properties, and maximum input/work bounds;
- assembly, platform, or FFI blind spots requiring manual or differential
  evidence; and
- residual untested behavior.

Files are still reviewed one at a time in the order in `STATUS.md`. A single
subsystem scan may support multiple review records, but each record must map
the evidence to that file's concrete contracts rather than cite a generic
green job.

## Bootstrap exit gate

The bootstrap is complete only when:

- all selected tools run reproducibly from documented commands or immutable CI
  definitions;
- their versions, rules, build profiles, and exclusions are recorded;
- all baseline results are triaged into findings, duplicates, or documented
  false positives;
- initial byte-boundary fuzzers run under ASan/UBSan with bounded deterministic
  CI smoke coverage;
- native AVX-512 and Apple Silicon gaps are stated rather than inferred from
  portable or simulated builds;
- no unresolved critical/high result blocks the first file-review group; and
- `scripts/audit/check_audit_manifest.ps1` passes.
