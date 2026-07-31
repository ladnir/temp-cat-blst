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
| AUD-001 | high | open | Rust pairing | Safe `Pairing` can retain a DST pointer beyond the slice lifetime | `bindings/rust/src/pairing.rs`, native pairing context | unassigned | Decision accepted: own DST in trailing context allocation, preserving the non-generic API and reallocating safely on growth | Temporary-DST, reinit, empty-DST, merge, threading, Miri/ASan tests |
| AUD-002 | high | open | Rust pairing | Safe `mul_n_aggregate` can pass an undersized scalar slice to native code | `bindings/rust/src/pairing.rs:132` | unassigned | Translate `96a46e6` with overflow-safe required-byte calculation | Boundary, extreme-length, G1/G2, and high-level API regressions |
| AUD-003 | medium | open | Rust MSM | Oversubscribed threaded Pippenger planner can divide by zero | `bindings/rust/src/pippenger.rs` | unassigned | Decision accepted: translate `6271918` routing/guards and audit the low-overhead per-worker result-slot design | Direct planner properties and injectable high-CPU G1/G2 tests |
| AUD-004 | high | validated | Key generation | Non-null empty v5 salt reaches nonzero-only assembly byte-copy loop | `blst/sign/keygen.cpp:46` and generated SHA assembly | Phase 1 | Guard the byte-copy call when the salt length is zero; retain null-pointer rejection separately | Exact-vector C and Rust tests passed under portable and assembly/ADX builds; see `PHASE_01.md` |
| AUD-005 | medium | open | Scalar MSM | Precompute stack workaround uses hidden heap allocation and silently zeroes output on failure | `blst/scalar/curve/multi_scalar.cpp:323` | unassigned | Decision accepted: implement bounded allocation-free chunking, retaining checked arithmetic and no heap fallback | G1/G2 threshold, layout, infinity, and instrumented allocation-free tests |
| AUD-006 | info | open | Scope removal | Unsupported C-over-ABI C++ wrapper remains tracked and package-listed | `bindings/blst.hpp`, its commented test, SWIG and Rust package include references | unassigned | Remove the wrapper and all packaging/build references; do not port `e61acc2`/`d4ca0f2` | Clean C/C++/Rust builds and repository-reference check after deletion |
| AUD-007 | medium | validated | Scalar arithmetic | `blst_256_t` right-shift operators invoke left shift | `blst/blst_t.hpp` | Phase 1 | Port `1519e46`, add compile-time portable dispatch, and preserve zero-shift identity | Independent modular-halving oracle passed under portable and assembly/ADX builds; see `PHASE_01.md` |
| AUD-008 | high | validated | Scalar assembly | `Fr::to_mont` places possibly unreduced input in the operand position with tighter assembly bounds | `blst/scalar/fields/Fr.h` conversion path | Phase 1 | Translate the `bf008ec` operand contract by placing the reduced constant first | Boundary, aliasing, round-trip, and 4,096-wide-input tests passed under portable and assembly/ADX builds; see `PHASE_01.md` |
| AUD-009 | high | open | Generated assembly | Checked-in assembly cannot currently be reproduced from the relocated generators | `build/refresh.sh`, `blst/asm`, generated Win64 macro divergence | unassigned | Repair and split deterministic generation before assembly ports | Two byte-identical regenerations and strict output allowlist |
| AUD-010 | medium | open | Windows ABI | Non-portable MSVC objects unnecessarily publish/reference fallback `$1` assembly symbols | x86 generators and generated Win64 objects | unassigned | Decision accepted: use one explicit portable-mode MASM macro policy and conditionally expose fallback symbols after reproducibility baseline | MSVC builds plus `dumpbin /symbols` checks in portable/non-portable modes |
| AUD-011 | medium | open | Uniqueness tree | Fixed traversal arrays are indexed without an explicit depth bound; size/failure arithmetic is unaudited | `blst/common/rb_tree.cpp` and C/Rust allocation callers | unassigned | Decision accepted: make a breaking capacity-aware API with explicit inserted/duplicate/error results; remove the unsafe legacy interface | Large-tree, 32-bit, overflow, over-insertion, alignment, duplicate, invalid-length, and resource-failure tests |
| AUD-012 | low | validated | Decoding tests | Retained tests lack direct compressed-infinity coverage with correct G1/G2 input sizes | scalar C/C++, Rust, and AVX public decode wrappers | Phase 1 | Translate the test intent from `17226b2` | G1 48-byte and G2 96-byte tests passed through scalar C, native C++, Rust, and portable AVX-512 public paths; see `PHASE_01.md` |
| AUD-013 | medium | validated | Native C++ | `blst_t.hpp` was not standalone after the fork's type refactor, while instance metadata inflated field-element layout and made size-based methods unsafe | `blst/blst_t.hpp` | Phase 1 review | Restore explicit compatibility with namespaced fork primitives and make compile-time metadata truly static | Direct standalone inclusion, 256/384-bit layout assertions, and portable/assembly right-shift tests passed; see `PHASE_01.md` |



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
