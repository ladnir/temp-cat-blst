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

## Phase 2 — Rust ownership and threaded MSM

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
- Harden CodeQL, sanitizers, warnings, and supported build matrices.
- Begin CatBLST rename and two-prefix ABI work only after the upstream
  implementation rows are closed.
