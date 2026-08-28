# File Review: `blst/avx/curve/scalar_mult.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/avx/curve/scalar_mult.h` |
| Roadmap phase / queue | Phase 16 AVX-512, reviewed early with the Phase 14 group-multiplication trust chain |
| Reviewed commit | `b6d832c421b8c569955db12941204df0e71d7972` |
| Reviewer | Codex audit team |
| Date | 2026-08-04 |
| Disposition | in-scope production; review complete; cross-file `AUD-091` caller and `AUD-218` stack ownership remediated |
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
reachable only through the documented subgroup MSM. CatBLST exposes no
unchecked MSM family; only upstream's two inherited single-point
`unchecked_mult` exceptions retain the generic fixed-window path.

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
adds independent subgroup/reduced MSM entries. Permanent full-width
non-subgroup regressions are retained only for the two inherited single-point
`unchecked_mult` functions.

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
The inherited CT G2 GLS construction duplicated base and endomorphism tables
across source arrays, wrappers, and the final owner. `AUD-218` now builds all
four tables directly in final storage, reducing measured peak local stack from
about 102 KiB to about 55 KiB without heap allocation or variable-size storage.

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

`AUD-027`, `AUD-039`, `AUD-088`, `AUD-089`, `AUD-090`, `AUD-091`, `AUD-218`.

## Validation

- `cmake --build out/build/x64-Release-Scalar --target blst_group_scalar_mult_tests --config Release --parallel 1`
- `out/build/x64-Release-Scalar/test/blst_group_scalar_mult_tests.exe`: 8/8, 55 ms.
- `cmake --build out/build/x64-Debug-Avx --target blst_group_scalar_mult_tests --config Debug --parallel 1`
- `out/build/x64-Debug-Avx/test/blst_group_scalar_mult_tests.exe`: 8/8, 6.342 s.
- Peach GCC 15.2 native AVX-512/IFMA Release focused build with `-j1` and test:
  8/8, 34 ms.
- Peach `nm -S` records unchanged CT multiply symbol sizes: G1 `0x13e9`, G2
  `0x3033`.
- Peach Clang 21 MatrixNoK G2 evidence at `b6d832c`: 16/16 focused tests pass;
  the CT owner frame is `0xc400`, the unchanged executor frame is `0x1700`,
  and one serial run measures CT GLS at 39.911 us versus 40.969 us before the
  ownership repair.
- Structural searches find no TODO/FIXME, dead `__has_include`, unused beta
  inverse/type aliases, ambiguous lineage marker, or unbounded table template.
- No broad suite, sanitizer, CUDA, BTx/BTE, or CI run was started; only the
  documented serial focused G2 benchmark process ran for `AUD-218`.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Existing C/C++ baseline | Prior copy finding `AUD-027` validated | Endomorphism algebra and timing require manual/oracle evidence |
| Clang Static Analyzer / clang-tidy | Phase 12 baseline | Width finding `AUD-039` validated | Not rerun for comment/include/assert cleanup |
| Clippy | Existing supported Rust boundary | N/A for this C++ template | Rust reaches the checked C boundary |
| ASan/UBSan | Not run in this packet | Pending milestone gate | Focused bounds are Debug-asserted; sanitizer evidence remains useful |
| Fuzzing | Not yet present | Coverage gap | Structured scalar/point and MSM differential fuzzing remains planned |

## Residual risk and decision

Complete for this file at the reviewed source commit. The single-point code has
independent scalar, portable-AVX, and native AVX evidence. The high-severity
`AUD-091` caller defect is validated:
ordinary MSM is subgroup-only and full-width decomposition canonicalizes
before using this file's splitters. Statistical-leakage, sanitizer, and fuzz
evidence remain subsystem release gates, not incomplete source-review work.

## Phase 16 dead-default addendum

Commit `5b08edc8f890ef1cbf4c6443162ec4da985879ee` removes the two unused generic
curve defaults from `GLVMult` and `GLSMult`. Every retained instantiation names
its selected MatrixNoK curve explicitly. The scalar decomposition, Booth/wNAF
executors, tables, and multiplication bodies are unchanged. Focused
portable-Debug compilation and 37/37 related tests pass. `AVX-P4-003` is
rebound to the exact new source hash, and `AUD-387` is validated.

## 2026-08-27 constant-time G1 MSM addendum

BTX partial decryption requires an MSM whose coefficients are secret. The
existing Pippenger family is deliberately variable-time and remains restricted
to public coefficients. This packet adds a separate fixed-window GLV executor.
It decomposes each canonical scalar into two 128-bit components, builds both
public point tables, scans every table row for every five-bit Booth window, and
uses complete point addition. Point count and chunk selection are public;
scalar contents do not affect control flow or memory addresses.

The executor shares one doubling schedule across fixed 32-, 16-, and 8-point
chunks. A bounded no-inline table walk is intentional. The first compile-time
16-component experiment caused the compiler to unroll a large selected-field
expression, consuming several minutes and multiple gigabytes. The bounded
runtime loop retains fixed public bounds, prevents that code-generation
failure, and keeps secret-dependent selection constant-time. Split scalars and
Booth indices are erased before return. Large public tables and split-scalar
staging are supplied by the caller; the template allocates nothing and does
not place the exact 394,240-byte (385 KiB) maximum scratch set on the thread
stack.

Native Peach measurements at batch 128 progressed from 3.131 ms for repeated
AVX multiplication to 1.647 ms with eight-point chunks, 1.548 ms with
sixteen-point chunks, and 1.489 ms with thirty-two-point chunks. Moving the
tables into caller scratch measured 1.486 ms, within noise. The final audited
source measured 1.4851 ms; its scalar reference measured 3.2865 ms under the
same serial Criterion policy.
The permanent Rust differential test covers empty input, every chunk boundary,
mixed chunk/tail sizes, and length/capacity rejection. C boundary tests add the
native-scalar oracle and failure-atomic scratch cases. Both modified fuzz
adapters compile with ASan and UBSan on Peach, and their native AVX and
middle-product corpus replays pass. Statistical leakage and broader sanitizer
evidence remain phase-wide residual gates.

Optimized Clang 21 disassembly of the selected 64-component specialization
(the 32-point GLV batch) was inspected on Peach. Its conditional branches are
loop or shape branches over the public component count, fixed window and limb
counters, and the public bit-offset case used to read a window spanning two
bytes. Booth carry and sign handling lower to `setcc`, `cmov`, shifts, and mask
operations. The existing constant-time table gather remains an all-row scan;
the disassembly contains no scalar-derived branch or scalar-derived table
address. This is a code-generation check for the selected compiler, not a
replacement for statistical leakage testing.

## 2026-08-27 public G1 MSM addendum

Public G1 MSM now has a fixed-batch GLV wNAF-5 executor. Each input scalar is
split into the same two 128-bit components as single-point multiplication. For
each component the executor builds eight odd multiples, records 129 public
wNAF digits, and shares one doubling schedule across the complete batch.
Scalar-derived branches, digit skips, signs, and table indices are intentional:
the entry point requires public scalars.

The wNAF storage is caller-owned and completely overwritten before use. It is
not zero-initialized on construction or erased on return because it contains
only public point tables and public digits. The executor is instantiated only
for fixed 32-, 16-, and 8-point chunks, preventing an unbounded template or
register-pressure expansion. Pippenger remains the large-input algorithm.

On Peach with Clang 21, one pinned core, and sequential Criterion runs, the
eight-point public MSM measures 73.164 us versus 203.54 us for projective
Pippenger and 102.97 us for the constant-time kernel. Shared wNAF remains
faster than Pippenger at 192 points (1.532 ms versus 1.618 ms); the two are tied
at 224 points (1.812 ms versus 1.815 ms). The production crossover therefore
selects wNAF through 192 points and Pippenger above it.

`TODO(AVX-formal: G1-VARTIME-MSM-01)`: after the dispatch and window policy
stabilize, add Lean validation of the GLV split, signed wNAF reconstruction,
endomorphism sign convention, and shared-doubling accumulation. The current
packet has differential, sanitizer, and fuzz evidence but does not claim a
formal proof of this new executor.
