# File Review: `blst/scalar/curve/multi_scalar_batch_affine.cpp`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/scalar/curve/multi_scalar_batch_affine.cpp` |
| Roadmap phase / queue | Phase 15 scalar curves and protocols, fourth packet |
| Reviewed commit | `913c52f6c226db69ee02f705c7d722f663850b27` |
| Reviewer | Codex audit team |
| Date | 2026-08-05 |
| Disposition | in-scope production |
| Subsystem | variable-time batch-affine subgroup MSM |
| Profiles | scalar, C, native C++, generated raw Rust; CUDA excluded |
| Generation/provenance | Hand-maintained Hypercat/CatBLST-local implementation derived from upstream Pippenger and curve formulas |
| Upstream references | pinned BLST `54e6e55674722fc2797ebb4bbb71b26d881eb4b8:src/multi_scalar.c` for extraction/reduction primitives; no matching upstream batch-affine file |
| Related files | scalar MSM header, C wrapper/declarations, generated Rust declarations, binding tests |

## Purpose and trust boundary

This translation unit evaluates contiguous G1/G2 subgroup MSM by batching
affine bucket additions behind one field reciprocal while routing repeated
bucket collisions through projective accumulation. Points and scalars are
trusted public data. The explicit C entry is an expert, allocation-free API.

## Contracts and invariants

- Points are valid, nonidentity prime-order subgroup elements.
- Scalars are public little-endian values with 1 through 256 meaningful bits.
- Windows lie in `[1,24]`; production selection uses 10 through 15.
- Batch sizes are compile-time 128 or 384.
- Scheduled point and bucket indexes fit `uint32_t`.
- Scratch is disjoint limb-aligned raw storage of at least the queried size.
- Empty input succeeds with infinity; other rejection leaves output unchanged.
- Release does not validate point membership; Debug asserts the narrow domain.

## Dependency and call-path map

Explicit G1/G2 C wrappers canonicalize native Montgomery Fr values into the
scratch prefix, then call the 128- or 384-entry instantiation. The selected
automatic MSM planner chooses this kernel only in its measured range and falls
back to ordinary Pippenger for affine infinity. Generated Rust exposes only
the unsafe raw C declarations; there is no safe expert wrapper.

## Upstream/specification comparison

Window extraction, Booth encoding, bucket integration, and complete XYZZ
operations descend from pinned upstream BLST. The affine batching scheduler is
a local Hypercat/CatBLST extension and has no exact upstream source file.

## Checklist evidence

### 1. Identity and purpose

Every implementation line, template instantiation, C wrapper, declaration,
planner call, generated Rust declaration, test, and local history was traced.

### 2. Contracts and invariants

Subgroup/nonidentity, public timing, scalar width, window, batch size, scratch,
overlap, empty, and failure contracts are explicit.

### 3. Upstream and specification lineage

The shared extractor and point formulas were compared with pinned upstream.
The local scheduler lineage and absence of an exact upstream file are recorded.

### 4. Correctness

The forward/reverse product chain, equal-point doubling slope, opposite-point
cancellation, signed affine addition, projective collision accumulation,
descending bucket sum, and high-to-low window composition were checked.

### 5. Memory and lifetime safety

All sizes use checked arithmetic. Raw scratch is partitioned by aligned byte
offset, every trivial internal object lifetime starts explicitly, and access
uses laundered byte-indexed references rather than fabricated typed arrays.

### 6. Undefined and implementation-defined behavior

Window shifts are limited to the 24-bit extractor domain; scalar widths stop
at 256; point and bucket indexes fit their stored width. Layout and triviality
requirements are statically asserted.

### 7. Side channels and cryptographic hygiene

The algorithm branches and indexes on scalar digits, bucket occupancy, and
point equality. The API is therefore explicitly variable-time and limited to
public data. It logs nothing and retains no scalar state.

### 8. Randomness and secret lifecycle

N/A: no entropy or secret input is accepted by contract.

### 9. API and ABI

No signature or symbol changes. Invalid query parameters return zero; failed
execution leaves output unchanged. Generated raw Rust inherits the C contract.

### 10. Concurrency

All state is caller-owned or automatic. There is no allocation, global state,
thread pool, or internal scheduling.

### 11. Resource exhaustion and denial of service

Scratch grows exponentially with the public window and is caller-provisioned.
The corrected maximum is 24; the production planner uses at most 15 and caps
its measured point ranges.

### 12. Portability and build behavior

The code is ordinary scalar C++23 and does not use architecture intrinsics.
The C scratch contract requires limb alignment. CUDA remains excluded.

### 13. Performance-sensitive structure

Batch sizes remain compile-time constants. The scheduler, product chain,
collision path, and bucket loop are unchanged. Release object inspection finds
no lifetime-start calls or stores and retains direct byte-stride addressing.

### 14. Tests and evidence

Focused differential tests cover G1/G2, batch sizes 128/384, windows 1/2/7/10/14,
scratch reuse, full-width and zero scalars, distinct/repeated points, explicit
equal/opposite cases, projective collisions, invalid windows, short scratch,
empty input, and oversized internal scalar widths.

## Findings

`AUD-121` through `AUD-125`.

## Validation

- MSVC scalar Debug and Release focused packets: 7/7 passed in each profile.
- Strict MSVC `/W4 /WX /O2` translation-unit compilation passed.
- Optimized object: zero scratch-wrapper, lifetime-start, placement-new, or
  allocation calls.
- Public documentation, group-contract, manifest, and structural checks pass.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Existing scalar baseline | AUD-122 lifetime pattern required hand review | Algebra and object lifetime are not proved |
| Clang Static Analyzer / clang-tidy | Existing baselines | No prior batch-affine alert | Does not prove group arithmetic |
| Clippy | N/A implementation | Raw Rust declarations traced | Native C++ implementation not analyzed |
| ASan/UBSan | Not rerun | Milestone evidence pending | Standard lifetime repair is structural |
| Fuzzing | Not yet present | Coverage gap | Differential MSM fuzzing remains planned |

## Residual risk and decision

The implementation review and focused remediation are complete at
`913c52f6c226db69ee02f705c7d722f663850b27`. Sanitizer, fuzz, Apple, and full
supported-matrix evidence remain milestone work.
