# File Review: `blst/scalar/curve/multi_scalar.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/scalar/curve/multi_scalar.h` |
| Roadmap phase / queue | Phase 15 scalar curves and protocols, third packet |
| Reviewed commit | `af09f66ee9097fe35617af33fbcbc7dfd4229237` |
| Reviewer | Codex audit team |
| Date | 2026-08-04 |
| Disposition | in-scope production |
| Subsystem | scalar fixed-table and Pippenger MSM declarations/kernels |
| Profiles | scalar, C, selected native C++, safe Rust; CUDA excluded |
| Generation/provenance | Hand-maintained C++ translation of upstream BLST MSM |
| Upstream references | pinned BLST `54e6e55674722fc2797ebb4bbb71b26d881eb4b8:src/multi_scalar.c` |
| Related files | implementation, bulk addition, C API, selected facade, Rust multipoint wrapper, focused MSM tests |

## Purpose and trust boundary

The header defines fixed-table sizing/gather contracts, Booth extraction,
Pippenger window selection, bucket operations, and native templates. All MSM
scalars and points are public under this variable-time API. Ordinary operands
are prevalidated prime-order subgroup elements; Release performs no membership
validation.

## Contracts and invariants

- Scalar widths are at most 256 meaningful little-endian bits.
- Fixed-table windows are in `[2,14]`; Pippenger windows are in `[1,24]`.
- Zero points or zero bits produce identity without reading input storage.
- Compact descriptors have a non-null head; later null entries mean contiguous
  continuation. Exact n-entry arrays need no sentinel.
- Scratch is caller-owned, aligned, and sized by the corresponding query.
- All sizing multiplication and shifts fail closed on overflow.
- Internal Debug assertions enforce extractor and nonempty-loop invariants.

The authoritative C declarations are in `bindings/blst_scalar.h`; the selected
native facade is `blst/curve/multi_scalar.h`; safe Rust is
`bindings/rust/src/pippenger.rs`. All three now state subgroup, timing, width,
scratch, descriptor, empty-input, and failure contracts.

## Dependency and call-path map

The C API calls these templates for G1/G2. Fixed-table gather uses reviewed
`BatchAffineWorkspace`; point arithmetic comes from `point.h`. Rust owns its
scratch and calls the C boundary. AVX has a separate Phase 16 implementation.
Experimental CUDA retains only conditionally exposed compatibility helpers.

## Upstream/specification comparison

Booth recoding, bucket integration, fixed-table shape, and heuristic structure
remain upstream-derived. Local hardening caps explicit windows, checks sizes,
defines zero/oversized behavior, removes dead dispatch/debug APIs, and corrects
descriptor traversal without changing the arithmetic schedule.

## Checklist evidence

### 1. Identity and purpose

Every declaration, inline kernel, helper, selected export, and supported caller
was inspected.

### 2. Contracts and invariants

Group domain, public timing, scalar width, windows, compact descriptors,
scratch, aliasing, empty results, and failure behavior are explicit.

### 3. Upstream and specification lineage

The retained recoding and integration loops were compared with pinned upstream;
local surface and safety divergences are recorded in AUD-114 through AUD-120.

### 4. Correctness

Exact descriptor, tile-offset, width, fixed-table, Pippenger, and individual
multiplication oracles pass for G1/G2.

### 5. Memory and lifetime safety

Queries use checked arithmetic. The dead unchecked allocator was removed.
Workspace lifetime repair is detailed in the implementation review.

### 6. Undefined and implementation-defined behavior

Window shifts are bounded, zero-width extraction is unreachable, and pointer
descriptors are never read past the final logical input.

### 7. Side channels and cryptographic hygiene

MSM is explicitly variable-time over public data. Production scalar-bit and
bucket logging was removed; no secret contract is implied.

### 8. Randomness and secret lifecycle

N/A: no entropy or secret input is accepted by this family.

### 9. API and ABI

Public C signatures are unchanged. Internal dead/debug helpers were removed.
The C, native C++, and safe Rust contracts agree.

### 10. Concurrency

All mutable state is caller-owned or automatic. There is no shared pool,
cache, or mutable global.

### 11. Resource exhaustion and denial of service

Queries reject overflow; explicit bucket growth is capped by window 24. Actual
resource provisioning remains the caller's responsibility.

### 12. Portability and build behavior

MSVC scalar Debug/Release and the full-Rust-binding focused build pass. CUDA is
experimental and was not compiled.

### 13. Performance-sensitive structure

Bucket loops, Booth representation, direct prefetch, fixed widths, and typed
scratch layout remain explicit. No runtime abstraction or per-point allocation
was introduced.

### 14. Tests and evidence

Focused C/C++ tests pass 15/15 in Debug and 15/15 in Release across the two
packets; focused safe Rust tests pass 3/3.

### Specialized additions

No assembly generator changed. No benchmark was run because the arithmetic
kernel did not change and benchmark execution requires separate approval.

## Findings

`AUD-026`, `AUD-103`, `AUD-109`, and `AUD-112` remain linked prior findings;
`AUD-114` through `AUD-120` are validated by this review packet.

## Validation

- MSVC 19.50 scalar Debug: 15 focused tests passed.
- MSVC 19.50 scalar Release: 15 focused tests passed.
- Rust full bindings: layout 1/1 and G1/G2 multipoint 2/2 passed.
- Structural scans: no dead dispatcher, production print parameter, scalar
  print call, unbounded window, or post-final descriptor advance.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Existing scalar baseline | No remaining debug-code alert in this file | Algebra and exact descriptor protocol required hand review |
| Clang Static Analyzer / clang-tidy | Existing baselines | No additional blocker | Does not prove field/group arithmetic |
| Clippy | Full Rust bindings focused compile/test | 3/3 focused tests pass | Raw generated FFI remains Phase 18 |
| ASan/UBSan | Not rerun | Milestone evidence pending | Standard object lifetime is addressed structurally, not proved by sanitizer |
| Fuzzing | Not yet present | Coverage gap | MSM boundary fuzzing remains planned |

## Residual risk and decision

The file is complete at the recorded commit. Portable-AVX, sanitizer, fuzz,
and full supported-matrix evidence remain
milestone gates rather than blockers for this scalar file review.
