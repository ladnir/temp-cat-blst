# File Review: `blst/backend/backend_portable.cpp`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/backend/backend_portable.cpp` |
| Roadmap phase / queue | Phase 13, portable-backend trust chain |
| Reviewed commit | `314c186149b9cc95d85a0c54e2dd651f137e4236` |
| Reviewer | Codex audit team |
| Date | 2026-08-02 |
| Disposition | `in_scope_production` |
| Subsystem | Fixed-width portable field arithmetic and inversion/division support |
| Profiles | scalar C++, retained C/Rust consumers through compile-time scalar selection |
| Generation/provenance | Handwritten fork file translated from upstream BLST portable arithmetic |
| Related files | `backend_portable.h`, `backend.h`, scalar `Fp`/`Fr` fields, assembly oracle, `test/Interface_tests.cpp` |

## Purpose and trust boundary

This translation unit implements the fixed-width portable arithmetic used as a
no-assembly backend and as a differential oracle for optimized assembly. It
contains Montgomery arithmetic, modular add/subtract/shift, byte conversion,
constant-time inverse and square testing, signed approximation loops, repeated
squaring, and quotient helpers. Inputs frequently contain secret field or key
material; public loop counts and fixed template widths must not mask unsafe
value-dependent memory behavior.

## Contracts and invariants

- Instantiated widths are the fixed 64-bit-limb BLS12-381 widths declared in
  the header; unsupported 32-bit configuration is rejected before this file.
- Arithmetic is allocation-free, nonthrowing, and uses only caller arrays plus
  fixed automatic arrays.
- Constant-time routines use fixed iteration counts and laundering/masks where
  required. Variable-time/public helpers must be identified by their callers.
- Carry and borrow are propagated through the dedicated wide helpers, never by
  pre-adding into a limb that can wrap.
- `cc_sqr_n_mul_mont_384` requires nonzero count, permits `ret == a`, and
  forbids `ret == b`; the latter is documented and asserted in Debug with no
  Release check.
- The quotient helper consumes an `n+1`-limb dividend/remainder buffer and may
  correct an overestimated quotient by one.

## Review evidence by checklist area

### Identity, lineage, and dependency map

Every definition and explicit instantiation was read against its declaration,
wide primitive, scalar-field caller, and available assembly/differential test.
The file has no generated sections and no runtime backend selection.

### Correctness and undefined behavior

`AUD-053` found three lost-borrow sites: two in `cc_div_3_limbs` and one in
`cc_quot_rem_n`. Each now passes the prior borrow through the helper's carry-in
parameter. The exact maximum-limb quotient regression distinguishes the old
and corrected arithmetic. Existing division and inverse oracles remain intact.
`cc_num_bits` zero reachability through inverse approximation was validated
with the header fix. Incorrect return-value, sign-mask, SHA-copy-size, shift,
and alias documentation was corrected.

### Memory, lifetime, API, and concurrency

All scratch storage has compile-time extent. There is no heap allocation,
exception source, retained pointer, callback, global mutable state, or
thread-local state. Array and alias preconditions remain internal contracts;
public input validation belongs at retained API boundaries.

### Side channels and secret lifecycle

The reviewed constant-time inverse/square loops retain fixed iteration counts,
mask-based selection, and fixed memory access. The borrow repair changes only
carry propagation. The unreferenced general limb printer was removed so this
production object no longer exposes an accidental secret-dump helper.

### Portability and diagnostics

The exact translation unit passes strict MSVC, GCC, and Clang compilation after
removing the inherited broad Clang suppression from the header. Experimental
CUDA consumers are outside this host translation unit's release claim.

### Performance-sensitive structure

Hot arithmetic remains explicitly fixed-width and inlineable with stack arrays
and predictable loops. The fixes introduce no allocation, runtime dispatch,
copy, or new Release branch. The repeated-square alias assertion disappears
under `NDEBUG`. No benchmark was warranted because production instruction
structure did not materially change.

## Findings

- `AUD-036` - validated: dead assignments and broad `optnone` were removed or
  narrowed before this complete review.
- `AUD-053` - validated: all three lost-borrow expressions repaired.
- Adjacent low-risk cleanup: dead secret printer removed and contracts corrected.

## Validation

- Exact maximum-limb quotient regression returns quotient `UINT64_MAX-1` and
  remainder two.
- Portable inverse output matches assembly for zero, one, and
  modulus-minus-one; existing 256/384 inverse oracle tests compile unchanged.
- Directly affected `Interface_tests.cpp` and backend objects compile under the
  scalar assembly configuration.
- Exact post-cleanup translation unit passes MSVC `/W4 /WX` and GCC/Clang
  `-Wall -Wextra -Wpedantic -Werror`.
- `git diff --check` and the audit manifest validator pass at closure.

## Residual risk and decision

Complete for this file at the reviewed commit. Exhaustive arithmetic proof is
not claimed: larger property, sanitizer, UB, and fuzz campaigns remain Phase
12/14 release evidence. The native assembly implementations and field-level
composition retain their own file reviews and differential gates.

## Phase 14 conditional-negation follow-up

The file was reread at `314c186` for `AUD-060`. Portable modular conditional
negation now branchlessly maps zero to false and every nonzero limb value to
true, matching x86-64 and AArch64 assembly. Direct 256/384-bit differential
tests pass for flags zero, one, two, and all-ones, including zero operands.
Optimized MSVC emits `neg/or/shr` normalization and no flag-dependent branch.
