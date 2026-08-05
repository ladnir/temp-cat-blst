# File Review: `blst/scalar/curve/ec_mult.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/scalar/curve/ec_mult.h` |
| Roadmap phase / queue | Phase 15 scalar curves and protocols, fifth packet |
| Reviewed commit | pending remediation commit |
| Reviewer | Codex audit team |
| Date | 2026-08-05 |
| Disposition | in-scope production |
| Subsystem | fixed-window and public-scalar wNAF multiplication |
| Profiles | scalar host; shared AVX wNAF helper; CUDA excluded from the gate |
| Generation/provenance | C++ translation of upstream BLST fixed-window macros plus CatBLST hardening and wNAF support |
| Upstream references | pinned BLST `54e6e55674722fc2797ebb4bbb71b26d881eb4b8:src/ec_mult.h`; low-order fix `7c535f1` |
| Related files | scalar point implementation, scalar window tests, group multiplication tests, AVX scalar multiplication |

## Purpose and trust boundary

This header implements fixed-window scalar multiplication, constant-time table
gathering, small fixed-window multi-point evaluation, shared window extraction
and Booth encoding, and public-scalar wNAF recoding. It is internal machinery
behind G1/G2 `mult`, `mult_vartime`, and the two inherited `unchecked_mult`
operations.

## Contracts and invariants

- Fixed-window template widths lie in `[2,8]`; production uses four and five.
- `get_wval` extracts at most nine bits and `get_wval_limb` at most 25; callers
  mask loaded high bits and provide every addressed byte.
- Booth input contains at most `w+1` bits and uses a positive width below the
  limb shift limit.
- Fixed-window scalar storage contains `ceil(bits/8)` readable bytes and is
  disjoint from group output; exact point/output aliasing is supported.
- `mult_window_n` accepts explicit point/scalar heads or a null continuation
  after a head, and its output is disjoint from scalar/table storage.
- Fixed-window gather and accumulation are constant-time for fixed public
  widths/counts. wNAF recoding and its executor are variable-time and accept
  only public scalar data.
- Heap allocation is used only when the optional table is absent; overflow or
  allocation failure returns false without changing output.

## Dependency and call-path map

`Point::mult` and `unchecked_mult` use `mult_window` for their generic scalar
paths. Constant-time GLV/GLS constructs stack tables and calls the precomputed
`mult_window_n` form. `Point::mult_vartime` and native AVX public multiplication
share `recode_wnaf_vartime`. Scalar MSM reuses the extraction and Booth helpers.

## Upstream/specification comparison

The extractor, Booth recoding, constant-time gather, fixed table construction,
and complete accumulation structure match pinned upstream. CatBLST replaces
upstream stack allocation in the optional multi-point mode with checked heap
allocation and adds the selected public wNAF recoder. Upstream `7c535f1`
requires complete `dadd` accumulation so low-order unchecked inputs remain
correct; the retained implementation matches that change.

## Checklist evidence

### 1. Identity and purpose

Every helper, template instantiation, scalar/AVX caller, retained test, and
relevant upstream/local history entry was traced.

### 2. Contracts and invariants

Window, bit, source extent, masking, point/table, alias, descriptor, failure,
and timing preconditions are explicit.

### 3. Upstream and specification lineage

The current C++ structure was compared line by line with pinned upstream and
the later low-order repair.

### 4. Correctness

Booth overlap, zero digits, signed table selection, infinity propagation,
low-order accumulation, multi-point window ordering, and wNAF carry handling
were reviewed and tested.

### 5. Memory and lifetime safety

Compile-time table sizes are bounded and allocation multiplication is checked.
Default `malloc` storage is used only for implicit-lifetime point types; caller
tables contain live typed objects. Compact continuation requires a retained
null descriptor after a head when fewer than `npoints` explicit heads exist.

### 6. Undefined and implementation-defined behavior

Typed shifts replace signed literal shifts. Fixed templates reject widths that
would under-size tables, overrun extraction, or truncate selection indexes.
Debug assertions protect runtime extractor and Booth preconditions.

### 7. Side channels and cryptographic hygiene

Fixed-window gather scans every table row and uses conditional copies. Loop
counts and addresses depend only on public sizes. wNAF branches and indexes on
digits and is restricted to public scalars. Secret-capable CT scalar copies
remain securely erased; public VT intermediates are not.

### 8. Randomness and secret lifecycle

N/A for entropy. Secret scalar handling belongs to the constant-time point
caller and retains its wipes.

### 9. API and ABI

No public symbol or signature changes. The internal fallible status is now
non-discardable. The two inherited unchecked C/Rust operations retain their
exact-scalar, on-curve/non-subgroup contract.

### 10. Concurrency

All state is automatic, caller-owned, or a per-call allocation. There is no
global state, thread pool, or internal parallelism.

### 11. Resource exhaustion and denial of service

Optional table allocation uses checked multiplication and returns failure.
Production decomposition supplies bounded stack tables and performs no heap
allocation.

### 12. Portability and build behavior

The implementation is ordinary C++23 scalar code. CUDA remains experimental;
host scalar and shared native-AVX wNAF instantiations are retained.

### 13. Performance-sensitive structure

Production width-four/five table sizes, fixed scans, accumulation order, and
inlining remain unchanged. Public VT no longer performs forced volatile wipes.
The precomputed status invariant adds no Release branch.

### 14. Tests and evidence

Focused tests cover allocated/caller/precomputed modes, invalid inputs,
allocation failure, widths two/eight, partial-byte wNAF, ordinary non-subgroup
points, and upstream order-11/order-13 low-order fixtures.

## Findings

`AUD-035`, `AUD-088`, `AUD-089`, and `AUD-126` through `AUD-130`.

## Validation

- MSVC scalar Release fixed-window packet: 6/6 passed.
- MSVC scalar Release non-subgroup/low-order packet: 2/2 passed.
- Selected CT/VT agreement, exact-alias, full-width boundary, and independent
  crossover-oracle packet: 4/4 passed (focused Release total: 12/12).
- MSVC `/W4 /WX /O2` owning translation-unit compilation passed.
- Optimized 64/128-bit wNAF recoder tails return without erasure loops.
- Optimized G1/G2 precomputed functions contain zero `mult_window_n` calls and
  zero allocation references; the generic optional-allocation COMDAT remains
  available for retained internal compatibility.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Existing scalar baseline | No additional report | Does not prove recoding algebra or constant time |
| Clang Static Analyzer / clang-tidy | Existing baselines | `AUD-035` allocation path previously repaired | Does not prove compiler table scans |
| ASan/UBSan | Prior `AUD-035` focused baseline | Allocated/caller/precomputed modes passed | New bounded template checks are compile-time |
| Fuzzing | Not yet present | Coverage gap | Differential scalar fuzzing remains planned |

## Residual risk and decision

The implementation and focused validation are complete pending the exact
remediation commit hash. Broader sanitizer, fuzz, Apple, and supported-matrix
work remains at milestone gates.
