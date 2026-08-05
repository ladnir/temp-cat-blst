# File Review: `blst/avx/curve/scalar_mult.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/avx/curve/scalar_mult.h` |
| Roadmap phase / queue | Phase 16 AVX-512, reviewed early with the Phase 14 group-multiplication trust chain |
| Reviewed commit | `2744efc5ba9718e0c79c790764ca7f8fdce22d9d` plus the review remediation described below |
| Reviewer | Codex audit team |
| Date | 2026-08-04 |
| Disposition | in-scope production; review complete; cross-file `AUD-091` caller remediated |
| Subsystem | AVX G1/G2 fixed-window and variable-time GLV/GLS multiplication |
| Profiles | portable AVX Debug and native AVX-512/IFMA Release |
| Generation/provenance | Hand-maintained C++ imported with the VROOM backend at `8b73c682`; subsequently adapted to CatBLST contracts |
| Upstream references | No direct supranational BLST counterpart; division primitives retain the pinned BLST backend contract |
| Related files | AVX curve/point/MSM headers, scalar point and backend division code, C/Rust bindings, `test/Group_scalar_mult_tests.cpp`, `test/avx/Avx_msm_tests.cpp` |

## Purpose and trust boundary

This header constructs fixed point tables, performs constant-time Booth scans,
implements public-scalar wNAF execution, splits canonical scalars into GLV/GLS
components, and applies the G1 beta and G2 minus-psi endomorphisms. Single-point
checked multiplication can receive secret scalars. The variable-time and MSM
callers require public scalars. All arithmetic storage is fixed stack or static
owner storage; this file allocates no heap memory and owns no global thread or
mutable cache.

## Contracts and invariants

- Supported single-point scalars have `0 <= nbits <= 256`; nonzero lengths
  require a non-null little-endian scalar.
- GLV/GLS multiplication requires an on-curve prime-order subgroup point and a
  scalar reduced modulo r before the z/z^2 split.
- Table widths are compile-time values 1 through 7. Production uses width four
  tables, five-bit Booth windows, G1 wNAF-5, and G2 wNAF-4.
- Booth indices never exceed the owned table, non-full final windows cannot
  carry, and wNAF magnitudes are odd and within the selected table.
- Constant-time table selection scans every row. Variable-time paths may skip
  zeros and index directly because their scalar is public.
- Constant-time division scratch, split values, schedules, and executor scalar
  arrays are erased after final use with the audited `vec_zero` primitive.
- Exact output aliasing is handled by the owning point boundary, not retained
  pointers in this header.

## Dependency and call-path map

`blst/avx/curve/point.h` canonicalizes checked scalars and calls `GLVMult` or
`GLSMult`; the selected C API converts between boundary and AVX representations.
The scalar backend provides the division and wNAF helpers. `multi_scalar.h`
also calls the split helpers for GLV/GLS Pippenger preparation. That caller
originally exposed `AUD-091`; it now canonicalizes full-width scalars and is
reachable only through the documented subgroup MSM, while explicit unchecked
MSM bypasses the endomorphism split.

## Upstream/specification comparison

The file has no direct upstream BLST source counterpart. It entered the tree
with the VROOM AVX backend and uses BLST's documented division identities:
`div_by_zz` accepts a 255-bit value and produces a 128-bit remainder and
quotient, while `div_by_z` splits a 128-bit component. Constants and signs were
checked against scalar beta/psi code. In particular, the optimized minus-psi
formula relies on `psi_y_i == -psi_y_r (mod p)`. Local CT/VT policy, wiping,
and public API contracts are CatBLST decisions.

## Checklist evidence

### 1. Identity and purpose

All declarations, templates, compile-time validation branches, direct callers,
tests, and local history were inspected.

### 2. Contracts and invariants

Table-width, scalar-width, carry, index, wNAF, subgroup, canonical-scalar, and
secret-lifecycle invariants are recorded above. Debug-only assertions now
encode the zero-overhead internal bounds.

### 3. Upstream and specification lineage

The vague `Modified from blst` marker was replaced with Apache/SPDX provenance
that identifies the VROOM import and CatBLST adaptation. The dead alternate
include layout was removed; the supported tree now has one canonical backend.

### 4. Correctness

The independent test oracle uses simple double-and-add for G1/G2 and basic
cyclotomic square/multiply exponentiation for GT. It covers all order boundary
values, both G1/G2 dispatch crossovers, non-byte-aligned lengths, and four
deterministic full-width random scalars. Single-point CT, VT, and unchecked
results pass under all selected backends. The subsequent caller remediation
adds independent subgroup/reduced and generic/exact-scalar MSM entries and
permanent full-width non-subgroup regressions.

### 5. Memory and lifetime safety

All table and scalar storage has fixed compile-time extent. Direct table lookup
occurs only in the public-scalar path and now has an asserted internal bound.
No pointer escapes and no allocation or exception source is present.

### 6. Undefined and implementation-defined behavior

Widths constrain every shift below the promoted operand width. Byte extraction
uses unsigned values. Backend division retains the established limb/byte union
pattern. Supported CatBLST hosts are 64-bit and little-endian. No dormant
alternate header layout remains.

### 7. Side channels and cryptographic hygiene

CT execution has fixed window counts, fixed table scans, and no scalar-indexed
load. Debug assertions disappear from Release. `AUD-090` records optimized
proof that the secret-derived stack objects are erased. Statistical leakage
testing remains a separate evidence upgrade.

### 8. Randomness and secret lifecycle

N/A for entropy. Caller-owned input scalar memory remains caller-owned. Every
explicit CT scalar copy and schedule owned here is erased; compiler-created
register copies remain outside the C++ object-lifetime guarantee.

### 9. API and ABI

The templates are internal and change no C layout or symbol. Public C/native
C++ documentation now states the 256-bit and null-pointer contracts. Rust's
safe conveniences already enforce both before entering C.

### 10. Concurrency

The file is reentrant. Constants reside in immutable function-local singleton
owners outside this header's call state. There is no global pool or mutation.

### 11. Resource exhaustion and denial of service

Work and stack use are bounded by compile-time widths. Native G1/G2 functions
retain their prior object-code sizes. The allocation-free C MSM boundary uses
caller-owned bucket and split scratch under the remediated `AUD-091` contract.

### 12. Portability and build behavior

MSVC portable-AVX Debug and GCC native AVX-512/IFMA Release compile and execute
the focused target. The compatibility `__has_include` branch for the abandoned
VROOM layout was removed. CUDA is outside this file and the release gate.

### 13. Performance-sensitive structure

Fixed-width tables, direct native field operations, compile-time selection,
full CT scans, and specialized wNAF widths remain intact. VT/MSM receive no new
wipe cost. Native Release G1/G2 multiply symbols remain exactly 0x13e9 and
0x3033 bytes after the review cleanup, so no benchmark was warranted.

### 14. Tests and evidence

The eight-test scalar-multiplication target is permanent and independently
checks the optimized single-point paths. The separate four-test AVX contract
target closes the former MSM blind spot with `2^256-1`, actual on-curve
non-subgroup points, one/two-point execution, and modulo-r subgroup comparisons.

### Specialized additions

Native AVX hardware evidence is required for changes to this file. Portable
AVX Release remains excluded due compile cost. Assembly generators are not
involved.

## Findings

`AUD-027`, `AUD-039`, `AUD-088`, `AUD-089`, `AUD-090`, `AUD-091`.

## Validation

- `cmake --build out/build/x64-Release-Scalar --target blst_group_scalar_mult_tests --config Release --parallel 1`
- `out/build/x64-Release-Scalar/test/blst_group_scalar_mult_tests.exe`: 8/8, 55 ms.
- `cmake --build out/build/x64-Debug-Avx --target blst_group_scalar_mult_tests --config Debug --parallel 1`
- `out/build/x64-Debug-Avx/test/blst_group_scalar_mult_tests.exe`: 8/8, 6.342 s.
- Peach GCC 15.2 native AVX-512/IFMA Release focused build with `-j1` and test:
  8/8, 34 ms.
- Peach `nm -S` records unchanged CT multiply symbol sizes: G1 `0x13e9`, G2
  `0x3033`.
- Structural searches find no TODO/FIXME, dead `__has_include`, unused beta
  inverse/type aliases, ambiguous lineage marker, or unbounded table template.
- No benchmark, broad suite, sanitizer, CUDA, BTx/BTE, or CI run was started.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Existing C/C++ baseline | Prior copy finding `AUD-027` validated | Endomorphism algebra and timing require manual/oracle evidence |
| Clang Static Analyzer / clang-tidy | Phase 12 baseline | Width finding `AUD-039` validated | Not rerun for comment/include/assert cleanup |
| Clippy | Existing supported Rust boundary | N/A for this C++ template | Rust reaches the checked C boundary |
| ASan/UBSan | Not run in this packet | Pending milestone gate | Focused bounds are Debug-asserted; sanitizer evidence remains useful |
| Fuzzing | Not yet present | Coverage gap | Structured scalar/point and MSM differential fuzzing remains planned |

## Residual risk and decision

The single-point code is reviewed and has independent scalar, portable-AVX,
and native AVX evidence. The high-severity `AUD-091` caller defect is validated:
subgroup and generic MSM contracts are explicit and full-width decomposition
canonicalizes before using this file's splitters. The file remains `in_review`
only for the scheduled statistical-leakage, sanitizer, and fuzz evidence
upgrades, not for a known scalar/domain correctness defect.
