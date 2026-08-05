# File Review: `blst/scalar/curve/point.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/scalar/curve/point.h` |
| Roadmap phase / queue | Phase 15 scalar curves and protocols, first packet |
| Reviewed commit | `07d79a48cc7d82560f5194f76a137e4ee2527d06` |
| Reviewer | Codex audit team |
| Date | 2026-08-04 |
| Disposition | in-scope production |
| Subsystem | scalar G1/G2 point representations, arithmetic, encoding, validation |
| Profiles | scalar, C, C++, generated Rust; safe Rust complete-add callers traced |
| Generation/provenance | Hand-maintained C++ refactor of upstream BLST point declarations and local extensions |
| Upstream references | pinned BLST `src/point.h`, `src/ec_ops.h`; pre-refactor CatBLST macro implementation before `0f46c46` |
| Related files | `point.cpp`, `bindings/blst_scalar.h`, `bindings/rust/src/bindings.rs`, `bindings/rust/src/lib.rs`, scalar C API and focused point tests |

## Purpose and trust boundary

The header owns the scalar Jacobian, affine, and XYZZ G1/G2 types and their
native operations. Raw field/point representations can cross supported C,
native C++, and Rust boundaries. Arithmetic expects caller-established curve
and subgroup invariants where documented; decoders validate canonical field
range and curve membership but not general subgroup membership.

## Contracts and invariants

- G1/G2 ordinary multiplication requires on-curve prime-order subgroup points.
- The two inherited `unchecked_mult` exceptions require on-curve inputs and
  apply exact scalars without subgroup endomorphisms.
- Unsuffixed `add` uses an incomplete formula. Finite inputs must be distinct
  and nonopposite; `dadd`/`add_or_double` is complete.
- Opposite finite points share affine x and have negated affine y.
- Jacobian infinity uses `z == 0`; affine infinity is the all-zero object.
- Exact whole-object output aliasing is supported by point arithmetic;
  `cswap`/`ccopy` require distinct objects and a 0/1 control bit.
- Conversion, encoding, decoding, validation, and XYZZ operations may be
  variable-time and consume public data.
- High-level native decoding leaves output unchanged on failure.
- `Pointxyzz::from_affine` requires a finite affine input; `equal` compares raw
  representations rather than mathematical points.

The authoritative C contracts are in `bindings/blst_scalar.h`. Generated raw
Rust inherits those contracts. Safe high-level Rust aggregation uses complete
add-or-double operations and therefore does not expose the narrow formula.

## Dependency and call-path map

`point.cpp` implements the declarations using scalar Fp/Fp2 primitives and
fixed-window/GLV helpers. The scalar C API converts exact boundary types and
publishes local results. Pairing/signing/hash-to-curve/MSM consume these types.
AVX has an independent point representation and implementation. CUDA consumers
are experimental and outside release evidence.

## Upstream/specification comparison

Core formulas and legacy C spellings follow upstream BLST. CatBLST exposes
native C++ types and subgroup-optimized CT/VT operations. The narrow add naming
is inherited and retained (`AUD-096`); CatBLST documents the precondition rather
than adding hot-path validation or renaming the compatibility surface.

## Checklist evidence

### 1. Identity and purpose

All declarations, representations, comments, direct supported callers, and
tests were inspected.

### 2. Contracts and invariants

Point domains, scalar domains, infinity encodings, timing, aliasing, decoding,
and XYZZ preconditions are explicit after `AUD-108`.

### 3. Upstream and specification lineage

The old macro/union representation and current typed implementation were
traced. The G2 identity regression is recorded in `AUD-104`.

### 4. Correctness

The undefined `mult_w5` declaration was removed. `add_n_dbl` now documents its
actual accumulator postcondition. Narrow/complete addition paths have focused
tests.

### 5. Memory and lifetime safety

Point and affine are independent live C++ object types. The reviewed header no
longer promises prefix-cast use. Remaining cross-file scratch casting is
isolated as `AUD-109`.

### 6. Undefined and implementation-defined behavior

The reserved double-underscore guard was replaced. CatBLST globally requires
64-bit limbs and supported little-endian targets.

### 7. Side channels and cryptographic hygiene

Constant-time claims are operation-specific. Secret multiplication and
conditionals retain fixed schedules under their stated canonical/domain
preconditions. Public conversion/encoding/validation is not claimed CT.

### 8. Randomness and secret lifecycle

N/A: no entropy. Scalar temporary wiping is implemented in the multiplication
helpers reviewed under `AUD-088`/`AUD-090`.

### 9. API and ABI

Native types remain the selected C++ definitions behind the exact public C
types. C declarations define narrow add and transactional decoder behavior;
generated Rust inherits them.

### 10. Concurrency

Types contain caller-owned state only; generator references are immutable.

### 11. Resource exhaustion and denial of service

Point operations are fixed-size. Public scalar widths are bounded to 256 bits
except separately reviewed internal algorithms.

### 12. Portability and build behavior

The header compiles in MSVC scalar Release/Debug. CUDA annotations are retained
but CUDA is experimental and not release evidence.

### 13. Performance-sensitive structure

Incomplete add, fixed-window multiplication, GLV/GLS, exact aliasing, and
fixed-size stack representations remain explicit. No validation, heap object,
virtual call, or runtime dispatch was added.

### 14. Tests and evidence

Six focused point tests pass under MSVC Release/Debug and Linux ASan. Public
documentation and group-contract scripts pass.

### Specialized additions

Scalar/AVX differential requirements remain attached to multiplication and
membership findings. No benchmark was warranted by the header-only contract
cleanup.

## Findings

`AUD-096`, `AUD-104`, `AUD-105`, `AUD-106`, `AUD-107`, `AUD-108`, and related
cross-file `AUD-109`.

## Validation

- `cmake --build out/build/x64-Release-Scalar --target blst_scalar_point_audit_tests --parallel 1`: passed.
- `cmake --build out/build/x64-Debug-Scalar --target blst_scalar_point_audit_tests --parallel 1`: passed.
- Both focused MSVC executions: 6/6 passed.
- Linux ASan/LSan focused execution: 6/6 passed.
- Public documentation and group-contract integrity scripts: passed.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Existing C/C++ baseline | No header-specific alert | Algebraic domains and timing require manual review |
| Clang Static Analyzer / clang-tidy | Existing baseline | No new blocker | CUDA and generated assembly excluded |
| Clippy | Existing supported profiles | Raw FFI inherits C contract | Safe high-level containing file remains pending |
| ASan/UBSan | GCC 13 ASan/LSan Debug | 6/6 passed; AUD-104 regression covered | UBSan not rerun for this packet |
| Fuzzing | Not yet present | Coverage gap | Point decoding fuzzing remains planned |

## Residual risk and decision

The header is complete at the recorded commit. The inherited narrow `add`
naming is accepted with explicit documentation. `AUD-109` belongs to the next
bulk-addition file review and does not change this header's declared contract.
Decoder subgroup-policy design remains separately tracked by `AUD-094`.
