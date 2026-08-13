# File Review: selected AVX vector facade

## Record

| Field | Value |
| --- | --- |
| Path | `blst/avx/cpu/vector/avx512ifma.hpp`; `blst/avx/cpu/vector/vector_impl.hpp` |
| Roadmap phase / queue | Phase 16 MatrixNoK primitive owner follow-up |
| Reviewed commit | `9e237456a946abe827aa55a8a122e4a040114ffc` |
| Reviewer | Codex audit team |
| Date | 2026-08-13 |
| Disposition | `in_scope_production` |
| Subsystem | Selected eight-lane AVX-512/portable-emulation vector operations |
| Profiles | native AVX-512/IFMA and portable AVX Debug |
| Generation/provenance | Hand-maintained CatBLST adaptation of the MatrixNoK construction from VROOM commit `d6ecd534` |
| Upstream references | `not_applicable: this selected MatrixNoK vector layer is not an upstream BLST owner` |
| Related files | `blst/backend/avx512/intrin_compat.h`, `blst/avx/rns/matrixnok_fixed.h`, `blst/avx/rns/matrixnok_fp.h`, `test/avx/Avx_matrixnok_primitive_tests.cpp`, `proof/avx/CatBLST/AVX/Lane.lean`, `proof/avx/CatBLST/AVX/IFMA.lean` |

## Purpose and trust boundary

The facade represents exactly eight unsigned 64-bit lanes. Native builds map
each value directly to one selected AVX-512 register. Portable Debug builds map
the same operations to the `vv512i` emulation selected by the compatibility
layer. The owner supplies lane arithmetic, IFMA accumulation, permutations,
fixed-width loads and stores, comparisons, and constant-flow selection to the
MatrixNoK field implementation.

The lane values can depend on secrets. Callers supply the arithmetic range
bounds that make delayed reduction and IFMA accumulation correct. The facade
does not validate those higher-level residue bounds.

## Contracts and invariants

- `AVXVector<limbs>` is well-formed only for `limbs == 8`. One value contains
  exactly one 64-byte selected vector.
- Array construction, `load`, `store`, and `to_array` access exactly eight
  consecutive `uint64_t` values. The pointed storage must remain valid for the
  complete operation.
- The default constructors leave storage uninitialized. Internal callers must
  overwrite every used lane before reading it.
- `mullo` and `mulhi` use the AVX-512 IFMA definitions: each multiplicand is
  truncated to 52 bits, and accumulator addition wraps modulo `2^64`.
- Immediate shift counts follow the selected intrinsic layer. `get_scalar_mask`
  requires a public count in `[0, 63]` and asserts that precondition in Debug.
- `extract<index>` requires a compile-time lane index in `[0, 7]`.
- `ctime_select(true, b)` returns the receiver. `ctime_select(false, b)`
  returns `b`. Its selected implementation has no source-level branch or
  scalar-indexed memory access.
- `is_zero` and `equal` test every lane. Their inputs may be secret, and their
  selected implementations have fixed source-level flow.
- The operations allocate no memory, throw no exceptions, and own no shared
  mutable state.

These classes are internal native C++ implementation types. They have no C or
Rust declaration surface and make no stable ABI promise.

## Dependency and call-path map

`vector_impl.hpp` is the stable include point and contains no independent
implementation. MatrixNoK RNS and field owners instantiate `AVXVector` only
with their compile-time ring width. Every retained ring width is eight. Curve,
GT, pairing, FFT, and middle-product owners reach these operations through the
RNS and field layers. Tests and code-generation probes also instantiate only
`AVXVector<8>`.

The selected `vv` layer calls native `_mm512_*` operations directly when
native AVX-512 is enabled. Portable Debug selects the separately implemented
`_cc512_*` operations. This file performs no native/portable runtime dispatch
and contains no native-to-portable conversion.

## Upstream/specification comparison

The original VROOM owner exposed a generic positive limb count although its
loads and stores moved complete eight-lane registers. CatBLST uses only the
eight-lane BLS12-381 MatrixNoK construction. The exact-width assertion makes
that actual construction the compile-time contract and avoids a larger type
redesign during this audit. Dead generic helpers were removed.

The IFMA arithmetic correspondence is recorded by `AVX-P0-001`. The Lean
model defines eight wrapping 64-bit lanes and proves the mapped low-lane and
IFMA statements under their stated no-wrap assumptions. Compiler intrinsic
lowering and processor execution remain external boundaries.

## Checklist evidence

### 1. Identity and purpose

Both headers now carry Apache-2.0, Category Labs, and applicable VROOM lineage.
The implementation and include-facade roles are explicit.

### 2. Contracts and invariants

The exact lane count, storage width, input storage, shift, extraction, IFMA,
selection, initialization, allocation, exception, and timing contracts are
recorded above. Higher-level residue ranges remain owned by the RNS types.

### 3. Upstream and specification lineage

VROOM commit `d6ecd534` is the construction source. CatBLST's selected `vv`
domain and exact eight-lane restriction are intentional local adaptations.

### 4. Correctness

All production and test instantiations resolve to eight lanes. Portable
runtime checks cover construction, zero, equality, both selection outcomes,
and IFMA. Existing MatrixNoK differential tests cover the arithmetic users.
`AVX-P0-001` binds the accepted Lean statements to the exact source hash.

### 5. Memory and lifetime safety

The prior generic contract could transfer eight lanes through an array with
fewer elements. The compile-time equality assertion removes that state. The
remaining loads and stores cover exactly their declared eight-element input or
output. No method returns a reference to temporary storage.

### 6. Undefined and implementation-defined behavior

The Debug assertion prevents an invalid source-language shift in
`get_scalar_mask`. Compile-time extraction rejects an invalid index. Unsigned
lane arithmetic states its modulo-`2^64` behavior. Native intrinsic semantics
remain governed by the selected compatibility owner.

### 7. Side channels and cryptographic hygiene

Selection, equality, and zero testing use vector masks. A focused optimized
MSVC probe contains no conditional branch, call, or scalar-indexed access for
these operations. Public shift counts can control immediate-operation paths in
the selected intrinsic layer.

### 8. Randomness and secret lifecycle

N/A: this deterministic arithmetic owner creates no randomness and does not
own secret buffers. Higher-level owners control secret lifetime.

### 9. API and ABI

N/A: the types are internal C++ templates with no public C or Rust ABI.

### 10. Concurrency

Every operation uses only value arguments and caller-owned output. The owner
has no mutable global state or lazy initialization.

### 11. Resource exhaustion and denial of service

Every value has fixed 64-byte storage. Operations have fixed work and allocate
nothing.

### 12. Portability and build behavior

The native profile requires the selected AVX-512/IFMA toolchain and CPU. The
portable implementation is a developer-only Debug profile. The focused
portable object compile and runtime probe passed under MSVC 19.50.

### 13. Performance-sensitive structure

The exact-width change removes generic loops and helpers without adding a
production abstraction. Optimized native probes lower selection to
`vpblendmq`, IFMA to `vpmadd52luq`, and equality to vector compares plus mask
tests. None contains a wrapper call or memory round trip.

### 14. Tests and evidence

The focused evidence is listed below. Broader MatrixNoK correctness and native
Peach evidence remain attached to the accepted P0-P4 correspondence rows.

### Specialized additions

The exact source SHA-256 is
`608bb008bdd147f8c5d4970a76350b1a6020806425ba3a40a3f40fab19f995ed`.
Changing the implementation reopens `AVX-P0-001`.

## Findings

`AUD-037`, `AUD-102`, `AUD-242`, `AUD-379`, and `AUD-385`.

## Validation

- Structural inventory: every retained `AVXVector` instantiation resolves to
  eight lanes.
- Portable AVX Debug: the focused MatrixNoK primitive test object compiled
  with one Ninja job.
- Portable AVX Debug: a small linked runtime probe passed construction, zero,
  equality, both selection outcomes, and IFMA checks.
- Native MSVC 19.50 `/O2 /Ob3 /DNDEBUG`: focused selection, equality, and IFMA
  probes contain no conditional branch or call. IFMA is a direct
  `vpmadd52luq`. The object SHA-256 is
  `441836bad8faa7548a33f2ed9ff14a017b9b84b21250e5aa45dd845822fd77b8`.
- Native Peach GCC 15.2 `-O3 -DNDEBUG -march=native`: selected equality is one
  `vpcmpq` and selected low IFMA is one `vpmadd52luq`. Both bodies are
  identical to direct-intrinsic probes and contain no call, spill, conversion,
  or memory access.
- No benchmark, broad build, CUDA, BTx/BTE audit, portable AVX Release build,
  sanitizer, remote run, or Apple job was performed.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Existing C/C++ baseline | `AUD-037`; no separate finding in retained exact-width operations | Intrinsic semantics, mathematical bounds, and constant-time lowering require other evidence |
| Clang Static Analyzer / clang-tidy | Existing curated baseline | `AUD-102`; warning ownership is narrow | Native intrinsic lowering and template instantiations are not fully modeled |
| Clippy | N/A | N/A | No Rust source |
| ASan/UBSan | Not rerun for this packet | Deferred milestone evidence | Native AVX and portable-emulation behavior require distinct profiles |
| Fuzzing | Not rerun for this packet | Deferred milestone evidence | Finite vector primitives are better covered by differential and formal tests |

## Residual risk and decision

The source review and focused checks validate the exact eight-lane owner. Lean
does not prove compiler lowering, CPU intrinsic behavior, or constant-time
machine execution. Those common P0 boundaries retain prior Peach execution
and current focused code-generation evidence. A new supported native-toolchain
run remains a release milestone, not a reason to keep this owner pending.
