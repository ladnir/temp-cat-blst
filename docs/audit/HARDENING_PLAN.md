# Production Hardening Implementation Plan

## Ground rules

- Keep BTx/BTE frozen and unmodified until its required future merge is ready;
  see `BTX_FREEZE.md`.
- Do not audit or run benchmarks as part of correctness phases.
- Never run two benchmark processes concurrently.
- Keep changes small enough to review and revert independently.
- Every fix receives a regression before its upstream row is closed.
- Use the risk-proportional tiers in `VALIDATION_POLICY.md`; reserve full
  supported-profile matrices for subsystem and milestone gates.
- Generator repair must reproduce the fork before importing upstream assembly
  changes.

## Phase 0 — checkpoint and baselines

1. Commit the audit infrastructure separately from production fixes.
2. Record the exact local source commit and supported build presets.
3. Capture hashes of all checked-in generated assembly.
4. Confirm existing scalar portable, scalar assembly, AVX-512, and Rust tests
   that are runnable on the available hosts.

Exit gate: cleanly attributable audit baseline; no generated output changed.

## Phase 1 — small, high-confidence correctness fixes

Implement independently:

1. `AUD-004`: guard zero-length salt before the assembly byte-copy primitive.
2. `AUD-007`: call the right-shift primitive from both 256-bit operators.
3. `AUD-008`: place the reduced constant in the assembly-supported operand
   position for Fr Montgomery conversion.
4. `AUD-012`: add correct G1/G2 compressed-infinity regressions.

Validation:

- targeted regression for every defect;
- portable and assembly scalar configurations where applicable;
- scalar/AVX public-path comparison for shared decoding/conversion;
- Rust formatting and targeted tests.

Exit gate: each fix has an isolated implementation commit and passing targeted
validation.

## Phase 2 — Rust ownership and threaded MSM (historical)

The Phase 2 scheduler hardening was completed and validated, then superseded
by Phase 8. Rust no longer owns a scheduler or threaded MSM path; the original
evidence remains relevant to the history of `AUD-001`, `AUD-002`, `AUD-003`,
`AUD-014`, and `AUD-015`.

### 2A. Pairing ownership

- Implement `AUD-001` with one allocation containing native context followed
  by owned DST bytes.
- Reuse allocation capacity on `init` when sufficient.
- On growth, allocate and initialize the replacement before swapping/freeing.
- Correct C/Rust ownership documentation.
- Implement `AUD-002` in the same file-level review, but retain a separable
  commit.

### 2B. Pippenger

- Implement mandatory `AUD-003` oversubscription routing and divisor guards.
- Translate the upstream per-worker result-slot design.
- Audit initialization, one-writer-per-slot proof, synchronization,
  false-sharing, and worker-panic behavior.
- Compose with bounded-channel upstream work; do not reintroduce a removed
  small-path channel.

Exit gate: temporary-DST safety tests, reinitialization tests, Miri/ASan where
available, planner property tests, mocked high-worker tests, and G1/G2
differential MSM.

## Phase 3 — assembly generator trust chain

### 3A. Reproducibility only

- Split assembly, Win64 export, and Rust-bindgen generation.
- Replace stale `src/asm` paths with explicit `blst/asm` inputs.
- Generate into a temporary tree.
- Add deterministic compare/verify mode and an exact output manifest.
- Encode existing intentional fork differences in generators.

Exit gate: two consecutive generations are byte-identical and the first
generation matches all current intended outputs. No upstream assembly
optimization is included yet.

### 3B. Portable-mode and symbol policy

- Standardize on one macro spelling and meaning.
- Pass it explicitly to MASM from CMake.
- Translate `AUD-010`/`41b511a` so `$1` fallback symbols are public/external
  only when the portable-dispatch configuration requires them.
- Verify symbol tables for both modes.

Exit gate: MSVC portable and non-portable builds pass, with expected symbol
allowlists.

### 3C. Accepted upstream assembly chain

- Apply accepted generator-source changes in dependency order.
- Regenerate once through the repaired tool.
- Review only the expected generated-file allowlist.
- Differential-test portable, non-ADX, ADX, and AVX-delegated paths.

Exit gate: generator/source/output provenance is complete and assembly upstream
rows can be closed.

## Phase 4 — allocation-free scalar MSM precomputation

- Implement `AUD-005` using bounded chunking.
- Retain the fork's checked arithmetic and API validation.
- Do not retain any heap fallback.
- Cover both G1/Fp and G2/Fp2, pointer-array and packed layouts, infinity, and
  scratch-threshold crossings.

Exit gate: instrumented tests prove no heap allocation; scalar results match
the reference across both sides of the scratch threshold.

## Phase 5 — uniqueness-tree hardening

State: complete. Implementation and validation evidence are recorded in
`PHASE_05.md`; `AUD-011` is validated.

### Accepted API direction

The project has no external users, so make a deliberate breaking C/Rust API
change rather than preserving an unsafe legacy entry point:

- `blst_uniq_sizeof(capacity)` returns 0 on invalid capacity or checked size
  overflow.
- initialization receives and stores `capacity`;
- insertion returns an explicit result distinguishing inserted, duplicate,
  capacity/depth exhaustion, and invalid input;
- remove the old capacity-blind initialization and boolean-only insertion API;
- migrate the retained Rust helper in the same change.

### Implementation

- Store and enforce node capacity.
- Check allocation-size multiplication/addition.
- Bound traversal before writing fixed path arrays.
- Reject message lengths that cannot be encoded safely in `len_n_colour`.
- Compare unequal lengths without narrowing `size_t` through `long`.
- Distinguish inserted, duplicate, capacity/depth exhaustion, and invalid input.
- Document message-buffer lifetime requirements.

Exit gate: 32/64-bit depth analysis, overflow and over-insertion tests, large
tree tests, and C/Rust wrapper validation.

## Phase 6 — cleanup and release controls

- [x] Remove `bindings/blst.hpp` and all active packaging, build, test, and
  public-documentation references (`PHASE_06A.md`).
- [x] Remove Go bindings, module/lint tooling, CI, CGo shims, and orphaned
  amalgamation (`PHASE_06B.md`).
- [x] Remove SWIG, C#, Emscripten JavaScript, Java, Node.js, and Python
  bindings (`PHASE_06C.md`).
- Keep CUDA opt-in and explicitly experimental; see `CUDA_EXPERIMENTAL.md`.
  Do not modify, audit, extract, or remove frozen BTx/BTE until its required
  merge is complete.
- Begin the standing CodeQL, sanitizer, fuzzing, warnings, and supported-build
  workstream in Phase 12 and expand it through the Phase 19 release gate.
- Begin CatBLST rename and two-prefix ABI work only after the upstream
  implementation rows are closed.

## Phase 8 — caller-threaded Rust runtime

State: complete. Implementation and validation evidence are recorded in
`PHASE_08.md`; `AUD-021` is validated and upstream `07b109d` is closed as
superseded.

- Remove every internal Rust pool, queue, barrier, CPU planner, and scheduler
  dependency.
- Reuse one sequential implementation across host and freestanding targets.
- Preserve optimized native arithmetic kernels and caller-managed composition.
- Enforce the no-internal-scheduler contract with a repository policy check.

Exit gate: default, full-binding, portable, and forced no-std Rust evidence;
deterministic validation errors; no production thread creation or scheduler
dependency.

## Phase 9 — CHERI inverse scratch confinement

State: complete with target-toolchain residuals. Implementation and focused
validation are recorded in `PHASE_09.md`; upstream `2e1948b`, `ed865dc`, and
`fe85dc3` plus `AUD-022` and `AUD-023` are validated.

- Bound both inverse routines' aligned scratch capability to 512 bytes.
- Preserve that bound while alternating between the two 256-byte halves.
- Harden complete-register translation for all affected CHERI synthetic
  instructions.
- Regenerate only the exact 10 affected ARM artifacts and leave `blst.def`
  unchanged.
- Repair the 256-bit inverse test's modulus and independent correctness
  oracle.

Exit gate: deterministic generation, exact changed-output allowlist,
byte-identical ordinary ARM projections, translator regression, and focused
inverse tests. Native AArch64 and Morello/CHERI assembly and execution remain
release evidence to collect when those toolchains are available.

## Phase 10 — Apple support boundary and experimental ARM target gates

State: complete. Implementation and validation evidence are recorded in
`PHASE_10.md`; durable policy is in `ARM64_CHERI_EXPERIMENTAL.md`.

- Retain non-Apple ARM64 and CHERI sources without claiming production support.
- Keep Apple Silicon macOS in supported scalar scope without an override.
- Reject other ARM64 and CHERI targets by default in CMake and Cargo.
- Require separate explicit acknowledgements for ARM64 and CHERI purecap.
- Repeat the restriction in generated configuration headers.
- Keep exact promotion evidence as the only route to removing the gates.

Exit gate: normal x86-64 and Apple Silicon policy paths remain enabled; default
non-Apple ARM64 and CHERI configurations fail with actionable messages;
explicit experimental configurations proceed with prominent warnings; Rust
cannot bypass the policy through a prebuilt native archive.

## Phase 11 - conformance vectors and release language

State: complete. C++/Rust scalar and AVX-portable generator, pairing, and
basic-signature evidence plus narrow draft language are recorded in
`PHASE_11.md`.

## Main production-audit program

Phases 0-11 were upstream catch-up, urgent hardening, scope control, and audit
preparation. Phases 12-20 are the main production-readiness program. Their
durable roadmap and current work volume are in `PRODUCTION_AUDIT_ROADMAP.md`.
Automated analysis is a standing workstream across every phase, not one phase
that is completed and left behind.

## Phase 12 - audit platform and initial baseline

Establish the shared analysis platform described in
`AUTOMATED_ANALYSIS_PLAN.md`: CodeQL, Clang Static Analyzer/CodeChecker,
clang-tidy, Clippy, ASan/UBSan, fuzz infrastructure, corpus/artifact handling,
and exact finding triage. Close upstream `68b58bd` and `f62244e` through the
same CI/provenance work; do not add WASI.

Exit gate: the initial tools execute reproducibly, every baseline result is
triaged, initial byte-boundary fuzzers run under sanitizers, and no unresolved
critical/high result blocks Phase 13. The tools remain active and expand in all
later phases.

## Phase 13 - fundamental representations and scalar backend trust chain

Review common types, limbs, byte/endian helpers, configuration primitives,
portable scalar backend, assembly generators, generated assembly, CPU feature
selection, and boundary conversions. Treat every generator and its outputs as
one trust chain. Add bounded conversion/serialization fuzzing, UB/bounds tests,
and independent integer oracles.

Exit gate: every phase file has a complete review record; generator/output
provenance is exact; scalar portable/assembly behavior agrees at boundaries;
and the phase analyzer, sanitizer, and fuzz gates pass.

## Phase 14 - scalar and extension fields

Review Fr, Fp, Fp2, Fp6/Fp12, reductions, inversion, square roots, hash-to-field
support, representation transitions, and lazy-reduction bounds. Preserve
optimization-friendly structure. Expand property tests and fuzzing around
modulus boundaries, canonical encoding, algebraic identities, and independent
big-integer/reference oracles.

Exit gate: all field files are complete, arithmetic properties and conversion
contracts pass across applicable scalar backends, and no phase finding remains
unresolved without explicit acceptance.

## Phase 15 - curves, pairing, signatures, and scalar MSM

Review G1/G2 point operations, serialization and subgroup checks, map/hash to
curve, GT/pairing, key generation, signature/aggregation, scalar MSM, and the
uniqueness tree. Expand malformed-input, group-law, pairing, signature,
resource-bound, aliasing, and state-machine fuzzing under sanitizers.

Exit gate: every scalar protocol path has manual, oracle/property, sanitizer,
and applicable fuzz evidence; public failure semantics and resource bounds are
documented; findings are resolved or accepted explicitly.

## Phase 16 - AVX-512 backend and scalar differential audit

Review AVX-512 primitives, RNS/reduction, field/curve/GT/pairing/MSM code,
conversion boundaries, instruction gating, generated inputs, and dispatch.
Differential-test every exposed operation against the reviewed scalar oracle,
including edge inputs, inactive lanes, alignment, and lazy-reduction limits.

Exit gate: AVX source files are complete; native AVX-512 hardware execution
covers the supported production path; and scalar/AVX differential campaigns
pass. Portable-AVX Debug remains available for routine compile and analysis
coverage, but portable-AVX Release is not an ordinary audit or release gate
because of its excessive optimized compile cost.

## Phase 17 - FFT and middle product

Review scalar and AVX FFT/middle-product domains, roots, permutations,
normalization, aliasing, sizing, scratch contracts, and representation
transitions. Build bounded simple reference implementations for differential
tests and fuzz invalid, boundary, and in-place cases.

Exit gate: all retained FFT/middle-product files and APIs are complete with
identity/oracle, sanitizer, fuzz, and scalar/AVX evidence.

## Phase 18 - public C, native C++, and Rust surfaces

Review public headers, ABI layouts and symbols, native C++ exposure, C API
contracts, Rust raw FFI and safe wrappers, unsafe invariants, ownership,
features, errors, no-internal-scheduler behavior, and packaging inputs. Add
cross-language differential tests and FFI-focused sanitizer/fuzz coverage.

Exit gate: every supported language surface is complete; declarations and
semantics agree; standalone consumers and feature matrices pass; ABI and
ownership residuals are explicit.

## Phase 19 - support files and release evidence

Review retained tests, vectors, build/configuration, generators and policy
scripts, CI, repository controls, documentation, command-line examples, and
packaging. Close coverage gaps found in Phases 13-18. Run the full supported
CodeQL/Clang/Clippy, sanitizer, fuzz-regression, scalar, AVX-512, Apple Silicon,
C/C++/Rust, and packaging matrix. Benchmarks remain excluded unless an explicit
performance gate is approved, and never run concurrently.

Exit gate: all remaining in-scope support files are complete, the full corpus
and release matrix pass, native Apple/AVX evidence is recorded, and no
unresolved release-blocking finding remains.

## Phase 20 - CatBLST rename and selectable C ABI prefix

Rename the project after the audited implementation and support surfaces are
under control. Implement build-time selection of exactly one public C prefix,
`blst_` or `catblst_`, through generators and headers with no runtime dispatch
or duplicate hot-path abstraction.

Exit gate: both prefix modes build and pass the supported matrix separately;
symbols, packages, documentation, and C/C++/Rust consumers match the selected
identity; no stale BLST-fork claim remains.
