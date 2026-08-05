# File Review: Safe Rust pairing and uniqueness helpers

## Record

| Field | Value |
| --- | --- |
| Path | `bindings/rust/src/pairing.rs` |
| Roadmap phase / queue | Phase 15 aggregate packet with prior Phase 2/5 repairs |
| Reviewed commit | `091dfb2cae096ee8ed801036e98827b3253b2e84` |
| Reviewer | Codex audit team |
| Date | 2026-08-05 |
| Disposition | in-scope production |
| Subsystem | safe Rust pairing context and message uniqueness wrapper |
| Profiles | safe Rust with full bindings over scalar or selected AVX library |
| Upstream reference | not applicable: local safe Rust ownership and typing layer |

## Ownership and type safety

`Pairing` owns one correctly aligned allocation containing the native context
followed by its DST bytes. Reinitialization reuses capacity or initializes a new
allocation before releasing the old one. Moving the owner preserves the native
relative-DST representation. Exclusive mutation is enforced by `&mut self`;
`Send` is sound because the allocation is stable and the native context has no
thread affinity, while `Sync` and `Clone` are intentionally absent.

Pairing orientation uses sealed compile-time traits. G1 public keys accept only
G2 signatures and G2 public keys accept only G1 signatures. Optional signature
contributions use `Option<&P::Signature>`; there is no runtime type inspection,
downcast, or null produced by a failed cast. Unsafe FFI dispatch remains hidden
behind safe monomorphized methods. Public scalar slices are length-checked before
FFI, and malformed widths return `BLST_BAD_SCALAR`.

The safe API documents validation flags, commit/merge sequencing, trusted raw
inputs, Miller-result semantics, and owned DST lifetime. The uniqueness helper
uses checked allocation sizing, verifies native alignment requirements, and maps
all non-insert outcomes to `false`.

## Evidence and disposition

Thirteen focused full-binding tests pass. They cover owned and empty DSTs,
capacity reuse/growth, moving across threads, committed merge, typed pairing
orientation, explicit missing signature, scalar boundary lengths, empty Miller
identity, and uniqueness results. Earlier native/portable Phase 2 tests validate
the owned-DST and scalar-length repairs; Phase 5 covers the uniqueness ABI.
`AUD-001`, `AUD-002`, `AUD-011`, `AUD-152`, `AUD-155`, and `AUD-159` are
validated. Broader feature combinations, sanitizer/fuzz evidence, and Apple and
native AVX-512 execution remain milestone work.

Complete Phase 15 file review at implementation commit
`091dfb2cae096ee8ed801036e98827b3253b2e84`.
