# File Review: hash-to-field headers

## Record

| Field | Value |
| --- | --- |
| Paths | `blst/scalar/fields/hash_to_field.h`; `blst/fields/hash_to_field.h` |
| Roadmap phase / queue | Phase 15 curves and protocols, hash-to-field packet |
| Reviewed commit | `ab212fe80784dc516ca995b6823db6227141c52c` |
| Reviewer | Codex audit team |
| Date | 2026-08-06 |
| Disposition | in-scope production |
| Subsystem | RFC 9380 XMD:SHA-256 expansion and BLS12-381 hash-to-field |
| Profiles | scalar implementation; selected native C++ facade; C and safe/generated Rust boundaries |
| Generation/provenance | CatBLST C++ evolution of upstream BLST hash-to-field code |
| Upstream/specification | pinned BLST lineage; RFC 9380 sections 5.3.1 and 8.8 |
| Related files | `bindings/blst_scalar.h`, `bindings/rust/src/primitives.rs`, `blst/scalar/curve/map_to_G.h`, `test/Hash_to_field_audit_tests.cpp` |

## Purpose and trust boundary

The scalar header implements RFC 9380 `expand_message_xmd` with SHA-256 and
reduces fixed 64-byte chunks into Montgomery-form Fp/Fp2 outputs. Messages,
augmentations, DSTs, output lengths, and all public XMD pointers are
caller-controlled. Hash-to-curve callers use compile-time output sizes; the
public XMD boundary accepts a runtime size through 8160 bytes.

The selected facade intentionally exports only the hardened byte-expansion
wrapper. Field reduction remains scalar-owned because no alternate selected
backend exists and exposing it as backend-selected would be misleading.

## Contracts and invariants

- XMD writes exactly the requested number of bytes, including a final partial
  digest; zero length is a no-op.
- SHA-256 XMD permits at most 255 digest blocks, or 8160 output bytes. Longer
  public requests fail closed by clearing the caller-provided output.
- Nonzero lengths require valid pointers. Message and DST pointers may be null
  only when their lengths are zero. Public input and output ranges do not
  overlap.
- DSTs longer than 255 bytes use the RFC oversize-DST transform.
- The implementation is allocation-free and nonthrowing. It uses fixed local
  SHA-256/chaining storage and a 32-byte final-partial buffer.
- Chaining XOR is defined over bytes. No typed limb access, alignment promise,
  union active-member assumption, or effective-type exception is required.
- `hash_to_field` is compile-time bounded to the XMD domain and writes actual
  Fp/Fp2 subobjects. Outputs are reduced Montgomery-form field elements.
- Rust `blst_scalar::hash_to` reduces 48 XMD bytes modulo r, returns `None` only
  for zero, and erases the expanded intermediate before returning.

## Dependency and call-path map

Hash/encode-to-curve templates call scalar `hash_to_field`, which calls the raw
exact XMD core with compile-time sizes. The C function calls the hardened
`blst_expand_message_xmd` wrapper. The selected C++ namespace re-exports only
that wrapper. Safe Rust scalar hashing calls the C XMD boundary and then the
scalar conversion boundary; generated raw Rust inherits the C declaration.
Frozen BTx/BTE keeps its existing selected-wrapper call and was neither changed
nor audited.

## Upstream/specification comparison

The SHA-256 state schedule and RFC construction remain the inherited design.
CatBLST removes the upstream assumption that hash-to-field always requests a
multiple of 32 bytes, eliminates limb casts over byte storage, rejects the
one-byte block-counter overflow domain, and makes the public wrapper exact and
allocation-free. Official RFC XMD and BLS12-381 field intermediates now provide
direct conformance evidence.

## Checklist evidence

1. Both headers, all declarations/aliases, C and safe Rust boundaries, and
   production callers were traced.
2. Pointer, length, overlap, exact-output, failure-publication, allocation, and
   exception contracts are explicit.
3. The construction and block-count bound were compared with RFC 9380.
4. Official 32/128-byte XMD, oversize-DST, and G1/G2 field-intermediate vectors
   pass.
5. An output pointer offset by one byte proves the operation requires no limb
   alignment; canaries prove exact 65-byte publication.
6. Source inspection finds no limb cast, union type punning, or heap operation.
7. XMD and hash-to-field process public transcript material with fixed SHA-256
   schedules for a public length; no secret-dependent table access was added.
8. The Rust-expanded scalar intermediate is explicitly erased.
9. C, selected native C++, and safe/generated Rust contracts are recorded in
   `PUBLIC_API_DOCUMENTATION.csv`.
10. No dynamic initialization, random state, runtime dispatch, or threading is
    introduced.
11. Work and storage are bounded by the public output length and RFC maximum.
12. Focused scalar MSVC Debug/Release and safe Rust checks pass.
13. Optimized MSVC emits a 624-byte fixed frame and inlined unaligned vector XOR
    with no allocator reference.
14. Broad sanitizers, fuzzing, Apple, and the full release matrix remain later
    milestone evidence.

## Findings

`AUD-174` through `AUD-179`.

## Validation

- MSVC scalar Debug passes 6/6 and Release passes 7/7 focused tests.
- Official evidence includes RFC XMD:SHA-256 32-byte, 128-byte, and oversize-DST
  vectors plus BLS12-381 G1 and G2 hash-to-field intermediates.
- Boundary cases cover zero output, exact 65-byte publication with unaligned
  destination, 8160-byte success, 8161-byte fail-closed output, and Release
  invalid pointer/length pairs.
- Full-binding Rust scalar hashing matches an independent expected reduction;
  the 48-byte intermediate is erased before the result branch.
- Optimized MSVC disassembly contains inlined `vmovdqu`/`vpxor`/`vmovdqu` for
  chaining XOR, a fixed 624-byte frame, and no malloc/calloc/realloc/free call.
- No benchmark, broad suite, sanitizer sweep, CUDA build, BTx/BTE work, Apple
  CI, portable-AVX Release build, or parallel heavy work ran.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Existing local C++ baseline | No packet-specific alert | Exact RFC semantics and C++ byte-object lifetime required direct review |
| Clang Static Analyzer / clang-tidy | Existing repository baselines | No packet-specific blocker | Optimized MSVC behavior inspected directly |
| Rust | Full bindings, focused unit | Pass; `AUD-177` | Does not prove machine-code erasure by itself |
| ASan/UBSan | Deferred milestone evidence | N/A | The invalid typed access was removed structurally; sanitizer coverage remains a gate |
| Fuzzing | Deferred subsystem campaign | N/A | Exact KATs and boundaries do not replace arbitrary-input campaigns |

## Residual risk and decision

The two headers are complete at the recorded implementation commit. The
external XMD contract deliberately forbids overlap instead of adding
overlap-handling overhead. Secure erasure ultimately depends on the separately
audited zeroization primitive and compiler evidence. Broader sanitizer, fuzz,
Apple, and supported-matrix runs remain release evidence rather than change-loop
requirements.
