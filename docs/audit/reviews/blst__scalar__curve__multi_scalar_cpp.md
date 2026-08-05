# File Review: `blst/scalar/curve/multi_scalar.cpp`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/scalar/curve/multi_scalar.cpp` |
| Roadmap phase / queue | Phase 15 scalar curves and protocols, third packet |
| Reviewed commit | pending remediation commit |
| Reviewer | Codex audit team |
| Date | 2026-08-04 |
| Disposition | in-scope production |
| Subsystem | affine conversion, fixed-table precompute/multiply, Pippenger execution |
| Profiles | scalar, C, selected native C++, safe Rust; CUDA excluded |
| Generation/provenance | Hand-maintained C++ translation of upstream BLST MSM |
| Upstream references | pinned BLST `54e6e55674722fc2797ebb4bbb71b26d881eb4b8:src/multi_scalar.c` |
| Related files | scalar header, bulk addition, point arithmetic, C/Rust wrappers, focused tests |

## Purpose and trust boundary

This translation unit converts projective batches, builds reusable affine
tables, evaluates fixed-table MSM, and integrates Pippenger buckets. Inputs are
public subgroup data. Its C callers require nonthrowing behavior and complete
or fail-closed output publication.

## Contracts and invariants

- Affine table rows contain `2^(wbits-1)` consecutive multiples.
- Only the final packed input array may overlap the precompute output table.
- The fast path uses unused table bytes as typed temporary storage.
- Remaining work attempts one allocation no larger than 144 KiB; failure uses
  an 8 KiB cold stack helper and cannot change the result.
- Raw storage is aligned and each Point/Field lifetime starts before access.
- Pippenger reads exactly n point/scalar descriptors and uses public data only.

## Dependency and call-path map

G1/G2 C wrappers instantiate this file. `point.h` supplies arithmetic,
`bulk_addition.h` supplies lawful fixed-table staging, and the selected native
facade re-exports only reviewed operations. Safe Rust reaches it through C.

## Upstream/specification comparison

Affine product-chain order, row generation, Booth digits, prefetch, and bucket
integration remain equivalent to upstream. Local code replaces C-style prefix
aliasing and dynamic stack allocation with valid C++ object lifetimes and a
bounded resource policy.

## Checklist evidence

### 1. Identity and purpose

All implementation lines, explicit instantiations, callers, scratch layouts,
and history were inspected.

### 2. Contracts and invariants

Table shape, input overlap, allocation fallback, scalar widths, descriptors,
scratch, subgroup, timing, and failure results are documented.

### 3. Upstream and specification lineage

The reverse product chain and chunk recurrence were checked against upstream.
The fallback deliberately uses the general affine converter because later
chunks do not have a unit-z first point.

### 4. Correctness

Independent row conversion and individual-multiplication oracles cover G1/G2,
windows 9 through 14, all supported layouts, exact descriptor arrays, and
nonaligned tiles.

### 5. Memory and lifetime safety

Real affine x subobjects hold general conversion prefixes. Byte-indexed aligned
workspace starts every Point and Field object explicitly. No unrelated array
cast remains. Allocation failure is checked before access.

### 6. Undefined and implementation-defined behavior

All shifts and size products are bounded. Placement/launder access is confined
to trivially destructible aligned Fp/Fp2 Point/Field objects; assertions pin
the relevant layout requirements.

### 7. Side channels and cryptographic hygiene

The implementation is variable-time and restricted to public data. Scalar
printing, pointer printing, and bucket hashing are absent from production.

### 8. Randomness and secret lifecycle

N/A: no randomness or secret inputs under this API contract.

### 9. API and ABI

No C signature changed. Invalid void-boundary inputs fail closed. Internal
debug/dead signatures are not supported ABI and were removed.

### 10. Concurrency

Temporary state is per call; malloc/free are paired locally. There is no global
state or thread pool.

### 11. Resource exhaustion and denial of service

Heap use is capped at 144 KiB and optional. Stack fallback is fixed at 8 KiB.
Work and output table size remain caller-selected and queries reject overflow.

### 12. Portability and build behavior

MSVC 19.50 Debug/Release compile and execute. The ordinary C++ object model no
longer depends on C prefix aliasing. CUDA remains experimental and untested.

### 13. Performance-sensitive structure

In-table reuse returns before allocation/fallback. MSVC Release places the two
8,344-byte probed frames only in non-inlined cold Fp/Fp2 fallback helpers;
primary precompute functions contain no stack probe. Placement lifetime starts
emit no clearing or construction copy.

### 14. Tests and evidence

Debug and Release cover both a heap-sized case and maximum-window cases that
must take the cap-forced fallback. Exact memcmp table equality is checked entry
by entry. Optimized code verifies a null allocation reaches that same helper.

### Specialized additions

No generated assembly or benchmark changed. A benchmark was unnecessary after
object inspection confirmed the fast path avoids the cold frame.

## Findings

`AUD-005`, `AUD-109`, and `AUD-113` through `AUD-120`.

## Validation

- MSVC scalar Debug: focused packet 15/15 passed, including heap-sized and
  cap-forced maximum-table layouts.
- MSVC scalar Release: focused packet 15/15 passed.
- Object symbols: malloc/free expected; no alloca, printf, print helper, or
  normal hash helper. `__chkstk` occurs only in the two cold fallback helpers.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Existing scalar baseline | AUD-120 cleanup completed | Object lifetime and algebra required hand review |
| Clang Static Analyzer / clang-tidy | Existing baselines | No additional blocker | Heap-null routing required source/object inspection |
| Clippy | N/A implementation | Safe caller tested separately | C++ implementation not analyzed by Clippy |
| ASan/UBSan | Not rerun | Milestone evidence pending | Focused standard-lifetime repair validated structurally |
| Fuzzing | Not yet present | Coverage gap | Descriptor/window fuzzing remains planned |

## Residual risk and decision

The file is complete pending insertion of the final remediation commit hash.
Sanitizer/fuzz and non-MSVC release-matrix evidence remain milestone work.
