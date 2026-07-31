# File Review: `blst/common/rb_tree.cpp`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/common/rb_tree.cpp` |
| Reviewed commit | `1e43e0e65c149ea58c999a27f243c7d676560823` |
| Reviewer | Codex |
| Date | 2026-07-30 |
| Disposition | `in_scope_production` |
| Subsystem | common uniqueness utility |
| Profiles | scalar, AVX-512, C, Rust |
| Generation/provenance | Handwritten C++ translation of upstream `src/rb_tree.c` |
| Upstream references | `06a5711`; upstream `src/rb_tree.c` at that commit |
| Related files | `blst/common/rb_tree.h`, `blst/c_api.cpp`, `bindings/blst_scalar.h`, `bindings/rust/src/pairing.rs`, `test/C_api_tests.cpp` |

## Purpose and trust boundary

This file implements an insert-only red-black tree used to decide whether
public message byte strings are distinct. Callers own both the aligned tree
storage and the message buffers. The tree retains pointers rather than copying
message bytes; it must never allocate, free, mutate, or outlive those buffers.

The implementation owns tree balancing, bounds, and the packed
`message_length << 1 | colour` representation. Public ABI translation belongs
in `blst/c_api.cpp`, not here.

## Contracts and invariants

- `blst_uniq_sizeof(capacity)` accepts only nonzero capacities whose exact
  `header + capacity * node_size` fits in `size_t`.
- Storage must meet `blst_uniq_alignof()` and contain at least the reported
  number of bytes. Initialization stores the exact capacity.
- Successful insertion consumes one node. Duplicate and error results do not.
- Duplicate detection precedes the capacity check, so a duplicate remains
  distinguishable when the tree is full.
- A null message pointer is valid only for a zero-length message.
- Lengths greater than `SIZE_MAX >> 1` are rejected before shifting.
- Stored message bytes must remain readable and unchanged until the final tree
  operation completes.
- The root is black; red nodes have black children; every root-to-leaf path has
  equal black height; ordering is lexicographic by bytes and then length.
- Concurrent operations on one tree require external synchronization.

## Dependency and call-path map

`blst/c_api.cpp` is the sole retained native ABI caller. Rust
`pairing::uniq` allocates aligned storage, initializes it with `msgs.len()`,
and maps all explicit native outcomes to its boolean high-level result. Scalar
and AVX-512 builds share this exact implementation; there is no backend
dispatch or duplicated vector implementation.

The tree calls no crypto backend and performs no allocation. It uses only byte
comparison and pointer/index manipulation.

## Upstream/specification comparison

Upstream `src/rb_tree.c` at `06a5711` retains capacity-blind initialization,
unchecked trailing-allocation arithmetic, pointer-width traversal arrays, a
boolean result, and `size_t`-to-`long` length narrowing. Upstream `06a5711`
guards only the separate Go reimplementation's 64-entry walk.

CatBLST deliberately translates that safety intent into the retained shared
C/Rust implementation and goes further: stored capacity, checked sizing,
alignment reporting, explicit outcomes, bounded traversal, safe length
comparison, and documented borrowed-message lifetime. The red-black rotations
and recolouring algorithm remain structurally unchanged.

The Apache-2.0/SPDX provenance header is retained.

## Checklist evidence

### 1. Identity and purpose

The exact path, commit, handwritten provenance, consumers, and responsibility
boundary are recorded above.

### 2. Contracts and invariants

All size, alignment, initialization, result, encoding, ownership, lifetime,
ordering, and synchronization rules are listed above and reflected in the
public C documentation.

### 3. Upstream and specification lineage

Compared directly with `06a5711:src/rb_tree.c`. The only algorithmic lineage is
the conventional red-black insertion routine; the local divergences are the
accepted `AUD-011` hardening contract.

### 4. Correctness

Allocation uses `offsetof(rb_tree, nodes)` plus checked multiplication.
Lexicographic length ordering returns `-1`, `0`, or `1` without narrowing.
The length shift is guarded. Traversal checks the fixed bound before writing
either path array. Existing rotations/recolouring were retained and exercised
by 32,768 ascending insertions followed by duplicate probes.

### 5. Memory and lifetime safety

The implementation performs no heap allocation. Initialization rejects null
and misaligned storage; insertion rejects null and misaligned tree pointers.
Capacity is checked before indexing the trailing node array. Message storage is
borrowed and its lifetime is explicit.

### 6. Undefined and implementation-defined behavior

The only pointer-to-integer conversion is the alignment check. On supported
32- and 64-bit ABIs, `size_t` has pointer width; the x86-32 object compile and
x86-64 MSVC/GCC executions validate those targets. Shift width, array indexing,
and signed narrowing are guarded or removed.

### 7. Side channels and cryptographic hygiene

Messages and equality outcomes are public protocol inputs. Lexicographic
comparison and early duplicate return are intentionally variable-time. The
tree must not be repurposed for secret-string comparison.

### 8. Randomness and secret lifecycle

N/A: this file consumes no randomness and handles no secret key material.

### 9. API and ABI

Result values are explicitly fixed at 0 through 4 in the internal and public C
enums and in generated Rust bindings. The old `blst_uniq_test` symbol is absent;
the archive exports `sizeof`, `alignof`, capacity-aware `init`, and `insert`.
The future prefix generator must rename these symbols with the rest of the C
surface without adding runtime dispatch.

### 10. Concurrency

There is no global, static-mutable, lazy, or thread-local state. Independent
trees are reentrant. A shared tree is mutable and requires external locking.

### 11. Resource exhaustion and denial of service

Capacity arithmetic is checked exactly and insertion cannot exceed it. The
walk uses 64 path entries. On 32-bit ABIs, a 16-byte node and 12-byte header
limit capacity to less than `2^28`, so red-black height is below 56 and the
array covers every allocatable tree. On 64-bit ABIs, exceeding depth 64
requires more than `2^32` nodes in the worst case; it returns
`BLST_UNIQ_DEPTH_EXCEEDED` before writing the arrays. Work remains
`O(message-prefix length * log n)`.

### 12. Portability and build behavior

The implementation compiled and ran with MSVC x64 and GCC 13 x64. GCC 13 also
produced an ELF i386 object with `-m32`. Public headers compiled standalone as
C11 and C++23 with warnings as errors. ARM64 and WebAssembly were not executed.

### 13. Performance-sensitive structure

Insertion remains allocation-free, iterative, and `O(log n)`. The 64-entry
arrays preserve the original 64-bit stack footprint while doubling the unsafe
32-bit limit from 32 to 64 entries. No virtual dispatch, callback, copy, or
message allocation was introduced. No benchmark was run because benchmarks
are excluded from this audit.

### 14. Tests and evidence

`CApi.UniqTree*` covers exact outcomes, duplicate-at-capacity, zero/overflow
capacity, exact overflow boundary, alignment, null/oversized messages, empty
messages, and 32,768 ascending inserts. Rust covers the raw enum contract and
safe helper. Linux ASan/leak detection, full scalar/portable suites, targeted
AVX-512, Rust native/portable/no-std, archive-symbol inspection, and 32-bit
compilation passed. Exact commands and counts are in `PHASE_05.md`.

### Specialized additions

Rust FFI allocation rounds the byte size safely, verifies that native
alignment fits `u64`, checks initialization, and exhaustively matches every
result. Assembly, FFT, and AVX-512 specialized arithmetic checks are N/A
because this backend-independent file contains none of those operations.

## Findings

`AUD-011` — validated.

## Validation

See `docs/audit/PHASE_05.md` for commands, profiles, result counts, symbol
checks, sanitizer coverage, and the 32/64-bit depth analysis.

## Residual risk and decision

Complete for the supported x86-64 scalar/AVX-512 profiles and the compiled
x86-32 arithmetic model. ARM64 and WebAssembly were not executed. Depth
exhaustion above 64 cannot be allocated in practical validation, but the
pre-write guard is direct and the 32/64-bit red-black bounds are documented.
Borrowed-message misuse remains a caller contract and is now explicit.
