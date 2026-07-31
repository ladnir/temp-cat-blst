# Audit Findings

This is the canonical finding register. Add findings as soon as they are
identified, even if the likely fix is small.

## Severity

- `critical`: practical compromise of secret material, signature validity, or
  production integrity.
- `high`: exploitable correctness, memory-safety, side-channel, or availability issue.
- `medium`: meaningful defect requiring specific conditions or defense-in-depth gap.
- `low`: limited-risk defect, portability issue, or hardening opportunity.
- `info`: maintainability, documentation, or verification improvement.

## Status

`open`, `fix_in_progress`, `needs_decision`, `accepted_risk`, `fixed`,
`validated`, or `closed`.

## Register

| ID | Severity | Status | Area | Summary | Introduced/affected | Owner | Fix/decision | Validation |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| AUD-001 | high | validated | Rust pairing | Safe `Pairing` can retain a DST pointer beyond the slice lifetime | `bindings/rust/src/pairing.rs`, native pairing context | Phase 2 | Own DST in the same aligned allocation after the native context, preserving the non-generic API and reallocating safely on growth | Temporary/empty DST, reuse/growth, merge, threading, and native/portable tests passed; see `PHASE_02.md` |
| AUD-002 | high | validated | Rust pairing | Safe `mul_n_aggregate` can pass an undersized scalar slice to native code | `bindings/rust/src/pairing.rs` | Phase 2 | Translate `96a46e6` with overflow-safe required-byte calculation and return `BLST_BAD_SCALAR` before FFI | Boundary and extreme lengths passed for both group orientations under native and portable builds; see `PHASE_02.md` |
| AUD-003 | medium | validated | Rust MSM | Oversubscribed threaded Pippenger planner can divide or subtract through a zero/underflowing window | `bindings/rust/src/pippenger.rs` | Phase 2 | Translate `6271918` routing, clamp all CPU-derived window widths, and use cache-separated per-worker result slots | Planner matrix and injectable 1,024-CPU G1/G2 tests passed under native and portable builds; see `PHASE_02.md` |
| AUD-004 | high | validated | Key generation | Non-null empty v5 salt reaches nonzero-only assembly byte-copy loop | `blst/sign/keygen.cpp:46` and generated SHA assembly | Phase 1 | Guard the byte-copy call when the salt length is zero; retain null-pointer rejection separately | Exact-vector C and Rust tests passed under portable and assembly/ADX builds; see `PHASE_01.md` |
| AUD-005 | medium | open | Scalar MSM | Precompute stack workaround uses hidden heap allocation and silently zeroes output on failure | `blst/scalar/curve/multi_scalar.cpp:323` | unassigned | Decision accepted: implement bounded allocation-free chunking, retaining checked arithmetic and no heap fallback | G1/G2 threshold, layout, infinity, and instrumented allocation-free tests |
| AUD-006 | info | open | Scope removal | Unsupported C-over-ABI C++ wrapper remains tracked and package-listed | `bindings/blst.hpp`, its commented test, SWIG and Rust package include references | unassigned | Remove the wrapper and all packaging/build references; do not port `e61acc2`/`d4ca0f2` | Clean C/C++/Rust builds and repository-reference check after deletion |
| AUD-007 | medium | validated | Scalar arithmetic | `blst_256_t` right-shift operators invoke left shift | `blst/blst_t.hpp` | Phase 1 | Port `1519e46`, add compile-time portable dispatch, and preserve zero-shift identity | Independent modular-halving oracle passed under portable and assembly/ADX builds; see `PHASE_01.md` |
| AUD-008 | high | validated | Scalar assembly | `Fr::to_mont` places possibly unreduced input in the operand position with tighter assembly bounds | `blst/scalar/fields/Fr.h` conversion path | Phase 1 | Translate the `bf008ec` operand contract by placing the reduced constant first | Boundary, aliasing, round-trip, and 4,096-wide-input tests passed under portable and assembly/ADX builds; see `PHASE_01.md` |
| AUD-009 | high | validated | Generated assembly | Checked-in assembly could not be reproduced from the relocated generators | `build/refresh.sh`, `blst/asm`, generated Win64 macro divergence | Phase 3A | Split generation, canonicalize generator inputs/outputs, encode the existing MASM macro spelling in source, generate only into an empty temporary tree, and enforce an exact output allowlist | Two generations of all 100 outputs were byte-identical and matched the checked-in tree; scalar 310/310 and AVX-default portable 142/142 tests passed; see `PHASE_03A.md` |
| AUD-010 | medium | validated | Windows ABI | Non-portable MSVC objects unnecessarily publish/reference fallback `$1` assembly symbols | x86 generators and generated Win64 objects | Phase 3B | Use `BLST_PORTABLE` as the build option and `__BLST_PORTABLE__` as the assembly guard; conditionally expose both consumer references and provider definitions | Portable mode has the exact 40-entry fallback allowlist and non-portable mode has none across all 13 MASM objects; focused tests passed in both modes; see `PHASE_03B.md` |
| AUD-011 | medium | open | Uniqueness tree | Fixed traversal arrays are indexed without an explicit depth bound; size/failure arithmetic is unaudited | `blst/common/rb_tree.cpp` and C/Rust allocation callers | unassigned | Decision accepted: make a breaking capacity-aware API with explicit inserted/duplicate/error results; remove the unsafe legacy interface | Large-tree, 32-bit, overflow, over-insertion, alignment, duplicate, invalid-length, and resource-failure tests |
| AUD-012 | low | validated | Decoding tests | Retained tests lack direct compressed-infinity coverage with correct G1/G2 input sizes | scalar C/C++, Rust, and AVX public decode wrappers | Phase 1 | Translate the test intent from `17226b2` | G1 48-byte and G2 96-byte tests passed through scalar C, native C++, Rust, and portable AVX-512 public paths; see `PHASE_01.md` |
| AUD-013 | medium | validated | Native C++ | `blst_t.hpp` was not standalone after the fork's type refactor, while instance metadata inflated field-element layout and made size-based methods unsafe | `blst/blst_t.hpp` | Phase 1 review | Restore explicit compatibility with namespaced fork primitives and make compile-time metadata truly static | Direct standalone inclusion, 256/384-bit layout assertions, and portable/assembly right-shift tests passed; see `PHASE_01.md` |
| AUD-014 | high | validated | Rust MSM | Unchecked scalar-size arithmetic and empty indexing can wrap or panic before native MSM calls, potentially admitting an undersized scalar buffer | `bindings/rust/src/pippenger.rs`, `bindings/rust/src/pippenger-no_std.rs` | Phase 2 review | Use overflow-safe ceiling/multiplication, slice-derived pointers, and explicit identity results for empty and zero-bit work | Overflow properties plus empty and zero-bit G1/G2 tests passed in std/no-std and native/portable configurations; see `PHASE_02.md` |
| AUD-015 | high | validated | Rust FFI buffers | Rust MSM constructed initialized `Vec<T>` lengths over uninitialized affine, scratch, and worker-result storage | `bindings/rust/src/pippenger.rs`, `bindings/rust/src/pippenger-no_std.rs`, tests | Phase 2 review | Keep FFI output and scratch storage as `MaybeUninit`, prove full affine initialization before reconstruction, and use uniquely written `UnsafeCell` result slots | Native, portable, no-thread, and no-std compilation plus focused G1/G2 tests and independent soundness review passed; see `PHASE_02.md` |
| AUD-016 | low | open | Windows C ABI packaging | Checked-in `blst.def` is stale relative to the supported C header and its deterministic generator | `build/win64/blst.def`, `bindings/blst.h`, `build/generate-win64-exports.sh` | unassigned | Decide the exported C surface and regenerate during supported C ABI/package hardening; do not mix the 233-to-318 export change into assembly reproducibility | Generator is deterministic; prove shared-library consumers and symbol allowlists after the CatBLST prefix policy is implemented |



## Finding template

### AUD-NNN: `<title>`

- Severity:
- Status:
- Affected files and lines:
- Affected versions/commits:
- Preconditions and threat model:
- Impact:
- Evidence/reproducer:
- Root cause:
- Proposed fix or decision:
- Performance implications:
- Tests and validation:
- Residual risk:
- Owner:
