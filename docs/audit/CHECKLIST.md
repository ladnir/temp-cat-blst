# Per-File Audit Checklist

Every in-scope file receives a review record. Every item below must be answered
with evidence or marked `N/A` with a reason. “Looks good” is not evidence.

## 1. Identity and purpose

- [ ] Record the exact repository path and reviewed commit.
- [ ] State whether the file is handwritten, generated, vendored, or derived.
- [ ] Identify the subsystem, supported profiles, owners, and public consumers.
- [ ] Describe the file's responsibilities and what must never be delegated to it.
- [ ] Identify generated inputs/outputs and verify reproducibility where applicable.

## 2. Contracts and invariants

- [ ] Document inputs, outputs, preconditions, postconditions, and error behavior.
- [ ] Document representation invariants, valid ranges, alignment, aliasing, and
      ownership expectations.
- [ ] Trace each externally reachable entry point to its callers.
- [ ] Check that debug-only assertions are not relied on for production safety.
- [ ] Check zero, one, maximum, malformed, non-canonical, infinity, and identity
      cases where meaningful.

## 3. Upstream and specification lineage

- [ ] Link applicable rows in `UPSTREAM_MANIFEST.md`.
- [ ] Compare the file with its upstream counterpart at the recorded upstream head.
- [ ] Explain intentional divergences and identify local commits that introduced them.
- [ ] Link the governing standards, papers, test-vector sources, or protocol text.
- [ ] Verify copied/generated code retains required license and provenance.

## 4. Correctness

- [ ] Check arithmetic boundaries, carries, borrows, reductions, overflow, and
      truncation.
- [ ] Check signedness, shifts, casts, endian conversions, and integer-width assumptions.
- [ ] Check loop bounds, termination, indexing, and length-unit consistency.
- [ ] Check state transitions, partial failures, and cleanup paths.
- [ ] Compare independent scalar and AVX-512 results where both exist.
- [ ] Identify unreachable, dead, duplicate, or contradictory logic.

## 5. Memory and lifetime safety

- [ ] Check all allocations, deallocations, stack usage, and allocation failures.
- [ ] Check bounds, alignment, object lifetime, initialization, and use-after-move/free.
- [ ] Check pointer provenance, aliasing, overlap, and restrict assumptions.
- [ ] Check integer calculations used for allocation and buffer sizes.
- [ ] Check caller-owned buffers, DST/message/scalar lifetimes, and retained pointers.
- [ ] Check sensitive temporaries are not copied or retained unnecessarily.

## 6. Undefined and implementation-defined behavior

- [ ] Check shift widths, signed overflow, unaligned access, strict aliasing, and
      invalid pointer arithmetic.
- [ ] Check union/type punning and object-representation assumptions.
- [ ] Check compiler intrinsics and inline assembly constraints.
- [ ] Check language-version and compiler-specific behavior.
- [ ] Check behavior under optimized builds and with assertions disabled.

## 7. Side channels and cryptographic hygiene

- [ ] Classify inputs as public, secret, or mixed.
- [ ] Check secret-dependent branches, addresses, loop counts, and instruction choices.
- [ ] Check variable-time algorithms are unreachable for secret inputs or explicitly named.
- [ ] Check error behavior and early exits do not disclose secrets.
- [ ] Check zeroization requirements and compiler-resistant implementation.
- [ ] Check fault, invalid-curve/subgroup, and non-canonical-input handling where relevant.

## 8. Randomness and secret lifecycle

- [ ] Identify entropy sources and failure behavior.
- [ ] Check domain separation, nonce uniqueness, rejection sampling, and bias.
- [ ] Check secret creation, copying, storage, exposure, and destruction.
- [ ] Check test-only deterministic randomness cannot enter production builds.

## 9. API and ABI

- [ ] Check C, C++, and Rust declarations agree on layout, calling convention, and semantics.
- [ ] Check public structures, enums, return codes, symbol visibility, and versioning.
- [ ] Check nullability, lengths, units, ownership, and thread-safety are documented.
- [ ] Check behavior for unsupported CPU features.
- [ ] Record implications for the future `blst_`/`catblst_` build-time ABI mode.
- [ ] Reject runtime prefix dispatch or other abstraction overhead in hot paths.

## 10. Concurrency

- [ ] Identify global, static, lazy, cached, and thread-local state.
- [ ] Check initialization races, publication, mutation, and reentrancy.
- [ ] Check callback and foreign-function boundaries.
- [ ] Check tests cover parallel use where concurrency is supported.

## 11. Resource exhaustion and denial of service

- [ ] Bound attacker-controlled allocation, stack, recursion, work, and queue growth.
- [ ] Check multiplication/addition used in capacity calculations for overflow.
- [ ] Check large-system, empty-batch, and degenerate workload behavior.
- [ ] Check cancellation and partial-progress behavior where applicable.

## 12. Portability and build behavior

- [ ] Check supported operating systems, compilers, architectures, and language standards.
- [ ] Check feature detection versus feature assumptions.
- [ ] Check warnings at the intended warning levels and warnings-as-errors configurations.
- [ ] Check configuration combinations cannot silently select experimental CUDA or removed bindings.
- [ ] Check headers are standalone and C headers compile as C.
- [ ] Check reproducible generation and clean-tree build behavior.

## 13. Performance-sensitive structure

- [ ] Identify hot paths and record expected complexity.
- [ ] Preserve fixed-width/unrolled loops, batching, inlining, stack layout, and
      predictable access patterns unless evidence supports a change.
- [ ] Check hidden allocation, virtual/type-erased dispatch, unnecessary copies,
      temporaries, and redundant conversions.
- [ ] Inspect scalar and AVX-512 code generation for critical kernels when changed.
- [ ] Record benchmark baseline and regression threshold for performance changes.
- [ ] Run only one benchmark process at a time.

## 14. Tests and evidence

- [ ] Map existing tests to contracts and invariants.
- [ ] Add deterministic regression tests for every fixed defect.
- [ ] Add negative, boundary, differential, and round-trip tests as applicable.
- [ ] Use sanitizers, static analysis, fuzzing, or formal checks where they add coverage.
- [ ] Record exact commands, toolchain versions, CPU features, and results.
- [ ] Record untested configurations and residual risk.

## Specialized review additions

### Assembly and generated assembly

- [ ] Review the generator and generated output as a single trust chain.
- [ ] Verify ABI preservation, unwind metadata, stack alignment, callee-saved
      registers, symbol visibility, and platform syntax.
- [ ] Verify carry chains, flag lifetimes, instruction availability, and fallback selection.
- [ ] Diff regenerated output and explain every change.
- [ ] Test against an independent implementation and inspect representative object code.

### AVX-512

- [ ] Verify CPUID/OSXSAVE/XCR0 gating before any AVX-512 instruction executes.
- [ ] Verify lane ordering, masks, inactive lanes, vector/scalar conversion, and alignment.
- [ ] Verify reduction bounds at every lazy-reduction boundary.
- [ ] Differential-test every exposed operation against scalar across edge and random vectors.
- [ ] Check throttling/down-clock impact and crossover thresholds for dispatch decisions.

### FFT and middle product

- [ ] Verify domain-size, power-of-two, root-order, and inverse-normalization assumptions.
- [ ] Verify bit-reversal/permutation and in-place aliasing behavior.
- [ ] Verify twiddle generation/storage and field-representation transitions.
- [ ] Test sizes 0, 1, minimum supported, maximum supported, and invalid sizes.
- [ ] Differential-test convolution/middle-product identities against a simple oracle.

### Rust

- [ ] Review every `unsafe` block with a stated safety invariant.
- [ ] Verify FFI pointer/length pairs, slice construction, lifetimes, Send/Sync, and ownership.
- [ ] Verify panic behavior across FFI and allocation-failure assumptions.
- [ ] Check `no_std`, feature combinations, dependency versions, and MSRV policy.
- [ ] Run formatting, clippy, tests, and Miri/sanitizers where applicable.

## Completion gate

A file is complete only when:

- all applicable checklist items have evidence;
- all upstream deltas are classified;
- all findings are resolved or explicitly accepted;
- tests and validation are recorded;
- residual risk is stated;
- the review record and manifest cite the same exact commit; and
- `scripts/audit/check_audit_manifest.ps1` passes.
