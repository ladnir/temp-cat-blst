# File Review: `blst/scalar/curve/bulk_addition.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/scalar/curve/bulk_addition.h` |
| Roadmap phase / queue | Phase 15 scalar curves and protocols, second packet |
| Reviewed commit | `bd45afa58630e7d235d4bbe048e8fc88bc5d5a09` |
| Reviewer | Codex audit team |
| Date | 2026-08-04 |
| Disposition | in-scope production |
| Subsystem | variable-time affine G1/G2 batch addition and fixed-table gather scratch |
| Profiles | scalar, C, selected native C++, safe Rust caller |
| Generation/provenance | C++ translation and hardening of upstream BLST batch addition |
| Upstream references | pinned BLST `54e6e55674722fc2797ebb4bbb71b26d881eb4b8:src/bulk_addition.c`; local scalar refactor history |
| Related files | selected wrapper, `multi_scalar.h`/`.cpp`, scalar C API, C/Rust declarations, focused point tests |

## Purpose and trust boundary

The file batches addition of public affine subgroup points using Montgomery
batch inversion. It also owns the typed staging record used by fixed-table MSM.
Public C callers can supply compact pointer arrays and raw limb scratch; neither
boundary may rely on a nonexistent C++ subobject or allow an exception across
the C ABI.

## Contracts and invariants

- G1/G2 operands are public prime-order subgroup points; execution is
  variable-time and Release performs no membership validation.
- Affine infinity is the all-zero value. Equal, opposite, and infinity inputs
  are handled through the complete mixed-add cases.
- Compact arrays begin with a non-null pointer. Later null entries mean
  contiguous continuation from the prior point.
- Empty or invalid-head batch input writes infinity. Exact output/input storage
  overlap is unsupported.
- Fewer than 16 points use complete serial addition. Larger calls make at most
  one nonthrowing allocation capped at 144 KiB; failure uses the same serial
  result rather than a partial or identity result.
- `BatchAffineNode<Field>` is a genuine affine subobject plus field workspace.
  Its x/y/work footprint exactly matches the old x/y/z scratch footprint, but
  the workspace is not a Jacobian coordinate.
- Raw C fixed-table scratch must be limb-aligned and sized by the public query.
  Node lifetime begins before gather and every later access is to that live
  object.

## Dependency and call-path map

The scalar C API routes `blst_p1s_add` and `blst_p2s_add` here; safe Rust
`MultiPoint::add` reaches those functions. Fixed-table MSM uses the node and
workspace types from `multi_scalar.h`/`.cpp`. The selected native wrapper
exports only `points_add`. CUDA is experimental and outside release evidence.

## Upstream/specification comparison

The field formulas, product-chain order, inversion schedule, 16-point
threshold, and 144 KiB performance cap follow upstream. Upstream C could treat
the prefix structs as C storage; the typed C++ translation could not cast a
live `Point` to an unrelated `AffinePoint`. The remediation restores the
intended physical layout with valid C++ subobjects and keeps the arithmetic
schedule unchanged.

## Checklist evidence

### 1. Identity and purpose

Every declaration, formula, branch, caller, scratch layout, and public wrapper
was inspected. Fixed-table callers were traced only to the boundary needed by
this review; their containing files remain pending.

### 2. Contracts and invariants

Timing, subgroup, infinity, compact-array, overlap, allocation, failure, raw
storage, and output contracts are explicit on native/C/safe-Rust surfaces.

### 3. Upstream and specification lineage

The implementation was compared line by line with pinned upstream
`src/bulk_addition.c` and local refactor history. Formula/product-chain order
is preserved.

### 4. Correctness

Serial oracles cover the threshold, odd counts, chunking, equal/opposite, and
infinity cases for Fp and Fp2. Raw fixed-table MSM matches individual products.

### 5. Memory and lifetime safety

`AUD-109` is closed by a real affine subobject and explicit raw-storage
lifetime start. Index-based workspace access avoids array pointer arithmetic
across separately constructed C-scratch objects. `AUD-111` closes unchecked
allocation and malformed-head pointer arithmetic.

### 6. Undefined and implementation-defined behavior

Static assertions require standard layout, trivial construction/copy/
destruction, exact offsets and size, limb-compatible alignment, and an integral
limb count for Fp/Fp2 nodes.

### 7. Side channels and cryptographic hygiene

This operation is explicitly variable-time and restricted to public points.
No secret schedule or conditional-copy guarantee is claimed.

### 8. Randomness and secret lifecycle

N/A: no entropy or secret input is accepted under the contract.

### 9. API and ABI

Public C signatures and scratch byte counts are unchanged. Safe Rust retains
its method surface. Internal helper names and scratch record types are not a
supported ABI.

### 10. Concurrency

All mutable state is caller-owned or per-call scratch. There is no global
allocator, cache, thread pool, or shared state.

### 11. Resource exhaustion and denial of service

Work remains linear in caller-supplied point count. Temporary heap usage is
capped at 144 KiB and failure preserves correctness through serial execution.

### 12. Portability and build behavior

Fp/Fp2 layout assertions and scalar Debug/Release builds pass on MSVC 19.50.
CUDA annotations remain incidental and experimental; WebAssembly is not a
supported release target.

### 13. Performance-sensitive structure

The batch threshold, product tree, one-inversion schedule, fixed node footprint,
and maximum scratch capacity remain explicit. Optimized code has direct
144/288-byte indexing and no construction store, conversion, or extra affine
copy. No runtime polymorphism or hidden per-point allocation was introduced.

### 14. Tests and evidence

The focused 11-test target passes in scalar Debug and Release. No benchmark was
run without approval; object-code inspection establishes that the representation
repair itself is zero-overhead.

## Findings

`AUD-109`, `AUD-110`, `AUD-111`, and `AUD-112` are validated. Caller tracing
opened `AUD-113` for the separate fixed-table precompute stack budget.

## Validation

- MSVC 19.50 scalar Debug focused execution: 11/11 passed.
- MSVC 19.50 scalar Release focused execution: 11/11 passed.
- Optimized `c_api.cpp.obj`: one checked allocation, direct serial fallback,
  affine-only staging copy, no construction clearing.
- Optimized `multi_scalar.cpp.obj`: workspace start/get arithmetic is inlined;
  accumulator uses direct typed-node strides and no affine conversion.
- Structural scans: no projective-to-affine scratch cast, unprefixed scratch
  macro, generic scratch-size macro, or selected helper re-export.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Existing C/C++ baseline | `AUD-109` came from hand review, not an alert | Object-model intent and algebraic schedule require manual review |
| Clang Static Analyzer / clang-tidy | Existing baseline | No additional file blocker | Does not prove field algebra or performance |
| Clippy | Existing supported profiles | Documentation-only safe wrapper change | Complete Rust containing file remains pending |
| ASan/UBSan | Not rerun for this packet | Deferred focused evidence | Sanitizers cannot prove the corrected C++ object model |
| Fuzzing | Not yet present | Coverage gap | Structured point/batch fuzzing remains planned |

## Residual risk and decision

The file is complete at the recorded commit. `multi_scalar.h`/`.cpp` remain
pending; `AUD-113` and the native fixed-table dispatch/sizing/error paths must
be resolved during that review. Full release-matrix and fuzz evidence remain
milestone work.
