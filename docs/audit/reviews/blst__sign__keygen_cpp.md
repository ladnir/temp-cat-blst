# File Review: Key generation, HMAC, HKDF, and EIP-2333

## Record

| Field | Value |
| --- | --- |
| Path | `blst/sign/keygen.cpp` |
| Roadmap phase / queue | Phase 15 scalar curves and protocols, ninth packet |
| Reviewed commit | `f36e35b3831ec7b527e767e8221c00b53ff739a7` |
| Reviewer | Codex audit team |
| Date | 2026-08-05 |
| Disposition | in-scope production |
| Subsystem | HMAC-SHA-256, HKDF-SHA-256, BLS versioned keygen, and EIP-2333 derivation |
| Profiles | scalar implementation; public C/C++ and safe/generated Rust boundaries |
| Generation/provenance | CatBLST C++ adaptation of pinned BLST `src/keygen.c` |
| Related files | `blst/c_api.cpp`, `blst/sign/sha256.h`, `bindings/blst_scalar.h`, `bindings/rust/src/lib.rs`, `test/Keygen_audit_tests.cpp` |

## Purpose and trust boundary

This file transforms caller-controlled key material, salt, and key information
into secret scalars. It also derives EIP-2333 master and child keys through a
fixed 255-block Lamport expansion. Inputs can contain secrets and every local
owning scratch object is therefore security-sensitive. The six retained C
functions are legacy `void` boundaries: they must remain allocation-free and
nonthrowing, publish only complete outputs, and clear output for the explicitly
handled invalid IKM/salt cases.

## Upstream compatibility contract

The pinned upstream comparison covers `blst_keygen`, `blst_keygen_v3`,
`blst_keygen_v4_5`, `blst_keygen_v5`, `blst_derive_master_eip2333`, and
`blst_derive_child_eip2333`. Their C names, signatures, valid-input output, and
version selection are unchanged:

- default keygen uses v4 and the default salt;
- v3 retains its legacy no-retry zero edge;
- v4.5 accepts an explicit salt, while null salt selects the default and
  ignores `salt_len`;
- v5 requires a non-null salt pointer but accepts a non-null zero-length salt;
- null `info` ignores `info_len`;
- IKM shorter than 32 bytes clears output; and
- EIP child output may exactly alias its parent input.

The shared C header is also the supported C++-callable API. The intentionally
removed legacy `blst.hpp` wrapper was not restored. Safe Rust preserves the
same algorithms while rejecting short IKM before FFI; generated raw Rust
inherits the C contract.

## Internal contracts and lifecycle

- `HmacSha256Context::tail` is active limb storage. Hash/copy operations access
  its object representation through `unsigned char`; pad operations use limbs.
- A null HMAC key is an internal reuse sentinel and requires zero key length.
- Raw HKDF expansion accepts at most 255 SHA-256 blocks, is failure-atomic for
  rejected shapes, and accepts a null output only for zero length.
- Versioned keygen produces a canonical nonzero scalar for v4 and later. V3
  produces a canonical scalar and intentionally does not reject zero.
- Public output is written only after derivation completes. Child derivation
  consumes the parent into local scratch before publishing, permitting exact
  parent/output aliasing.
- All internal functions are `noexcept`; the public C boundary performs no
  allocation and cannot receive a C++ exception.
- Debug verifies required output/input pointers and internal HMAC/HKDF shapes.
  Release retains the upstream narrow pointer contract and handled failures.
- Functional zero filling uses ordinary stores. The complete 416-byte keygen
  scratch and 8,576-byte Lamport scratch are securely erased on every path that
  uses them. The latter remains a bounded stack allocation.
- The v4-family loop repeats if reduction produces zero. That specified,
  negligible event remains a timing distinction; no stronger constant-time
  claim is made for the complete keygen operation.

## Correctness and compatibility evidence

1. RFC 4231 HMAC-SHA-256 cases cover ordinary and greater-than-block-size keys,
   including split streaming updates.
2. RFC 5869 HKDF-SHA-256 extract/expand output matches exactly.
3. HKDF zero output, the 8,160-byte maximum, 8,161-byte rejection, null output,
   null info, and failure-atomic storage are tested directly.
4. Frozen independent outputs cover v3, default/v4, v4.5 with nonempty salt,
   v5 with nonempty salt, and v5 with a non-null empty salt.
5. Compatibility regressions cover null-info length ignoring, v4.5 null-salt
   defaulting with ignored `salt_len`, and exact EIP child aliasing.
6. All four published EIP-2333 master/child vectors pass.
7. Scalar MSVC Debug passes 7/7 focused tests, including invalid-shape death
   tests; Release passes the six applicable functional tests.
8. Both full safe-Rust signature variants pass the focused empty-v5-salt test.

## Code-generation and surface evidence

MSVC 19.50 Release emits eight ordinary 64-bit stores for HMAC functional
initialization and only two zero stores for the high limbs not overwritten by
the 48-byte scalar import. No representation copy or volatile functional-zero
loop remains. The final secure erasures compile to fixed volatile byte loops of
416 and 8,576 iterations. Keygen retains a bounded stack frame, and the Lamport
helper retains its existing `0x21c0` frame.

The optimized keygen object has no heap-allocation symbol reference. HMAC,
HKDF, generic keygen, and Lamport helpers are Static object symbols; the C API
object exposes exactly the six bare upstream keygen/derivation symbols.

## Findings and disposition

`AUD-146` through `AUD-151` are validated. The prior allocation and C-boundary
repair `AUD-050` remains validated. No critical/high residual defect is open in
this file.

## Residual risk and deferred evidence

Sanitizer/fuzz coverage, non-MSVC optimized inspection, Apple execution, and
the full supported-language milestone matrix remain release evidence. Compiler-
created secret copies outside explicitly owned objects remain beyond the C++
abstract-machine erasure guarantee. The experimental CUDA path and frozen
BTx/BTE subsystem were not built, edited, or audited. No benchmark was run
because the inspected change removes volatile stores without changing the
cryptographic schedule.

## Completion

Complete Phase 15 file review at implementation commit
`f36e35b3831ec7b527e767e8221c00b53ff739a7`.
