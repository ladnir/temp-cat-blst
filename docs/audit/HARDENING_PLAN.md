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
- Harden CodeQL, sanitizers, warnings, and supported build matrices through
  the Phase 12 automated-analysis bootstrap.
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

## Phase 12 - automated audit bootstrap

State: planned. Detailed scope, cadence, and exit gates are recorded in
`AUTOMATED_ANALYSIS_PLAN.md`.

1. Establish version-pinned CodeQL baselines for retained C/C++ and Rust with
   real representative builds.
2. Establish a C/C++ compilation database and run Clang Static Analyzer through
   CodeChecker plus a curated correctness-focused clang-tidy policy.
3. Run Clippy over each supported Rust feature family and record intentional
   exceptions rather than weakening production code to satisfy style-only
   suggestions.
4. Establish Linux Clang ASan/UBSan scalar coverage, with portable code as the
   instrumented memory/UB oracle for assembly and AVX implementations.
5. Add the first bounded fuzz harnesses at public byte-parsing and conversion
   boundaries, seeded from retained vectors and adversarial edge cases.
6. Triage every result into `FINDINGS.md` or a documented false-positive record.
7. Close upstream `68b58bd` and `f62244e` through the same CI/provenance work;
   do not add WASI.

Exit gate: all configured tools execute reproducibly; versions, rule sets,
profiles, and blind spots are recorded; every baseline result is triaged; no
unresolved critical/high finding permits the dependent file-review chain to
advance; and audit-integrity checks pass.

## Phase 13 - dependency-ordered file audit with expanding automation

Begin only after the Phase 12 baseline is triaged. Review files in the order in
`STATUS.md`. Static-analysis results guide attention but never replace the
manual checklist. Add sanitizer, fuzz, oracle, and scalar/AVX differential
coverage when each subsystem's contracts become understood, then rerun the
complete relevant analyzer/test set once at the subsystem gate.

Exit gate: every in-scope production and support file has a review record tied
to an exact commit, all applicable automated evidence or reasoned `N/A` entries
are recorded, findings are resolved or accepted explicitly, and the manifest
integrity gate passes.
