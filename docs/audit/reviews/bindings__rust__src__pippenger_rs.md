# File Review: `bindings/rust/src/pippenger.rs`

## Record

| Field | Value |
| --- | --- |
| Path | `bindings/rust/src/pippenger.rs` |
| Roadmap phase / queue | Phase 15 scalar MSM supported caller review |
| Reviewed commit | `af09f66ee9097fe35617af33fbcbc7dfd4229237` |
| Reviewer | Codex audit team |
| Date | 2026-08-04 |
| Disposition | in-scope production |
| Subsystem | safe Rust G1/G2 multipoint wrapper |
| Profiles | Rust full bindings over scalar C API |
| Generation/provenance | Hand-maintained safe wrapper |
| Upstream references | not_applicable: local safe Rust wrapper |
| Related files | scalar C declarations/API, scalar MSM implementation, Rust generated bindings/tests |

## Purpose and trust boundary

The file owns affine vectors and safe scratch allocations for G1/G2 batch
conversion, addition, validation, and public-scalar MSM. Unsafe calls are
confined to initialized output publication and the raw C boundary.

## Contracts and invariants

- `MultiPoint::mult` is variable-time and accepts only public scalar data.
- Points must be prevalidated subgroup elements unless `validate` is called.
- Scalars are contiguous little-endian values of `ceil(nbits/8)` bytes each.
- `nbits <= 256`; short or overflowed buffers panic before FFI.
- Empty points or zero bits return identity without allocation or FFI access.
- FFI-written vectors use `MaybeUninit` and become initialized only after the
  C routine writes every element.

## Dependency and call-path map

The wrapper uses generated C declarations for affine conversion, Pippenger,
batch addition, membership, and point helpers. There is no internal scheduler
or global thread pool.

## Upstream/specification comparison

The upstream wrapper structure is retained. Local code removed the former
threading planner, adds checked layout and width rejection, documents narrow
subgroup/public contracts, and removes an unused tile macro parameter.

## Checklist evidence

### 1. Identity and purpose

The full file, macro expansion inputs, safe methods, tests, and C callees were
reviewed.

### 2. Contracts and invariants

Point domain, scalar layout/width, timing, panic, empty, allocation, and
validation behavior are documented.

### 3. Upstream and specification lineage

Local deviations are linked to AUD-003 and the current MSM contract findings.

### 4. Correctness

G1/G2 multipoint multiplication and scalar-layout overflow tests pass.

### 5. Memory and lifetime safety

Checked multiplication prevents allocation-length wrap. `MaybeUninit` vectors
are published only after complete C writes; empty paths do not form FFI heads.

### 6. Undefined and implementation-defined behavior

Unsafe blocks have explicit initialized-storage and buffer-size preconditions.
The C exact-array repair removes the former final-descriptor overread below the
safe wrapper.

### 7. Side channels and cryptographic hygiene

MSM is documented variable-time/public-only. The wrapper does not copy or log
scalar values beyond its required input buffer.

### 8. Randomness and secret lifecycle

N/A: secret scalar use is outside this API contract.

### 9. API and ABI

Safe method names are unchanged. Panics are documented; raw FFI inherits the C
contract. The unused internal macro argument was removed.

### 10. Concurrency

No threads, pool, global scheduler, or shared mutable state.

### 11. Resource exhaustion and denial of service

Point-vector and scratch allocations remain caller-size dependent but all byte
and capacity arithmetic is checked before allocation.

### 12. Portability and build behavior

The full-binding feature compiles on local Windows Rust; scalar C backend tests
pass. Apple remains a later supported-matrix gate.

### 13. Performance-sensitive structure

One scratch vector is allocated for nontrivial MSM and passed directly to C.
No initialization, worker scheduling, callback, or per-point allocation is
added.

### 14. Tests and evidence

Full-binding focused tests pass G1/G2 multiplication and layout overflow.

### Specialized additions

No benchmark or CUDA/AVX expert Rust surface is part of this file review.

## Findings

`AUD-003`, `AUD-040`, `AUD-095`, `AUD-110`, `AUD-111`, and `AUD-114` through
`AUD-116` (through the repaired C boundary).

## Validation

- `cargo test --features full-rust-bindings scalar_layout_rejects_overflow
  --jobs 1`: 1/1 passed.
- `cargo test --features full-rust-bindings multi_point::test_mult --jobs 1`:
  G1/G2 2/2 passed.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | N/A Rust | C boundary covered separately | N/A |
| Clang Static Analyzer / clang-tidy | N/A Rust | C/C++ implementation covered separately | N/A |
| Clippy | Existing multi-profile baseline | Prior unsafe scratch finding AUD-040 resolved | Lint not rerun for documentation-only macro cleanup |
| ASan/UBSan | Not rerun | C boundary milestone evidence pending | Rust ownership not fully modeled |
| Fuzzing | Not yet present | Coverage gap | Safe scalar-layout fuzzing remains planned |

## Residual risk and decision

The file is complete at the recorded commit. Generated raw bindings remain a
Phase 18 generator review; broader Rust and Apple validation remains a
milestone gate.
