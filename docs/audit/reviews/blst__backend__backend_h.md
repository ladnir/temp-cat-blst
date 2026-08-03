# File Review: `blst/backend/backend.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/backend/backend.h` |
| Roadmap phase / queue | Phase 13, scalar selector after portable backend |
| Reviewed commit | `dc07be802c39a3c04a17347a5cec64f74d95ecd7` |
| Reviewer | Codex audit team |
| Date | 2026-08-02 |
| Disposition | `in_scope_production` |
| Subsystem | Compile-time scalar backend selection, SHA-256 selection, fixed scalar division helpers |
| Profiles | Scalar C++, C/Rust consumers, host assembly and pure portable scalar |
| Generation/provenance | Locally composed header with upstream BLST `src/vect.h` and `src/vect.c` arithmetic lineage |
| Upstream references | Current upstream snapshot `54e6e55674722fc2797ebb4bbb71b26d881eb4b8`; local CUDA namespace refactor `4e8775d1`; scalar selector repair `33914f2` |
| Related files | `backend_portable.h/.cpp`, `backend_asm.h`, `blst_t.hpp`, scalar fields, SHA/keygen/hash-to-field callers, generated scalar assembly |

## Purpose and trust boundary

This header is the single compile-time boundary between pure C++ scalar
operations and supported host assembly. It also selects SHA-256 support
functions and provides fixed-size division helpers used by scalar decomposition.
Its operands can contain secret scalars, field elements, and key-derivation
state. Callers depend on fixed-width stack-only behavior, correct in-place
operation, compile-time-only dispatch, and no exception or allocation boundary.

## Contracts and invariants

- CatBLST uses 64-bit limbs exclusively.
- Default host assembly is selected only when `BLST_ASM` and
  `BLST_DEFAULT_ASM` are defined and device compilation is inactive.
- Otherwise the `cc_*` pure C++ implementation is selected. There is no
  runtime selector in this header.
- A zero-bit 256- or 384-bit modular left/right shift is identity. Exact
  in-place aliasing (`ret == a`) and disjoint output are supported; arbitrary
  partial overlap is not a contract.
- Nonzero shifts retain the underlying backend precondition and fixed modulus.
- SHA compression updates caller-owned state and receives complete 64-byte
  blocks. Supported callers pass a positive block count.
- Division helpers operate on their documented fixed-width arrays and baked
  BLS12-381 constants.
- No function allocates, throws, retains a pointer, performs I/O, or owns
  caller storage.

## Dependency and call-path map

- Scalar fields and `blst_t.hpp` call the `_v256` and `_v384` selectors.
- C and Rust consumers reach the same selectors through retained field, curve,
  signature, and pairing boundaries.
- SHA-256 update/finalization, key generation, and hash-to-field call the
  selected SHA wrappers.
- Scalar point decomposition calls `div_by_z` and `div_by_zz`; those call the
  selected fixed-width division primitives.
- MSM code calls `vec_prefetch` with a nonzero compile-time object size.
- The header calls declarations from `backend_asm.h` or implementations from
  `backend_portable.h/.cpp`. The experimental CUDA bridge is included only
  when `BLST_CUDA` is configured and is not release evidence.

## Upstream/specification comparison

The division loops and constants match upstream BLST `src/vect.c` at audited
head `54e6e55674722fc2797ebb4bbb71b26d881eb4b8`, with the unreachable 32-bit
limb expression simplified for CatBLST's 64-bit-only contract. Assembly
declarations correspond to upstream `src/vect.h`. The selected wrappers are a
local architectural layer introduced to preserve pure C++ and assembly builds
without runtime dispatch.

Commit `4e8775d1` accidentally changed the SHA compression selector from the
host/default-assembly rule to an unconditional portable call. `AUD-055`
restores the prior/upstream intent while retaining the device fallback.

## Checklist evidence

### 1. Identity and purpose

The entire file, its history, manifest row, direct includes, and every caller
of its shift, SHA, division, and prefetch functions were inspected.

### 2. Contracts and invariants

Fixed widths, selection macros, zero-count behavior, aliasing, block-count,
array-size, modulus, and ownership contracts are recorded above. Every
supported shift entry now enforces the zero identity at this boundary.

### 3. Upstream and specification lineage

The fixed division routines were compared with upstream `src/vect.c`; SHA
selection was compared with upstream `src/sha256.h` and local pre-CUDA history.
The file now carries explicit Apache-2.0 and local-composition provenance.

### 4. Correctness

All four selected shifts route zero to an exact fixed-width copy before either
backend. Native C++ 256-bit left/right and 384-bit left/right operators pass
in-place and out-of-place identity tests. Deterministic two-block SHA state
updates match between assembly, portable, and the selected wrapper.

### 5. Memory and lifetime safety

The file uses caller storage and fixed automatic values only. Zero-shift copies
use direct limb assignments so exact in-place operation is defined without a
`restrict`-qualified copy helper. No storage escapes.

### 6. Undefined and implementation-defined behavior

`AUD-051` removes every supported path that could pass zero to assembly loops
or portable functions with explicit nonzero assumptions. Compile-time width
and architecture selection does not depend on runtime data. The obsolete
32-bit conditional array expression in `div_by_z` is gone.

### 7. Side channels and cryptographic hygiene

Backend selection is preprocessing-only. Arithmetic wrappers introduce no
secret-dependent dispatch. Shift count and SHA block count are control inputs,
not secret field values; the zero-count check replaces an unsafe backend
precondition. Arithmetic constant-time properties remain owned by the selected
backend files and their separate reviews.

### 8. Randomness and secret lifecycle

N/A: this file creates no randomness. It neither persists nor logs secret
material and adds no temporary secret buffers.

### 9. API and ABI

These are internal C++ inline selectors, but supported C, C++, and Rust paths
depend on them. Public signatures and representations are unchanged. The
header guard changed only from a misleading internal name to
`BLST_BACKEND_H`.

### 10. Concurrency

The file owns no global, static mutable, cached, lazy, thread-local, or worker
state. All operations use caller-owned values.

### 11. Resource exhaustion and denial of service

`AUD-051` closes the zero-count assembly-loop denial-of-service path for every
selected shift caller. Loop bounds elsewhere are fixed widths, positive SHA
block counts supplied by checked callers, or public input lengths.

### 12. Portability and build behavior

The exact header compiles under MSVC 19.50 `/W4 /WX` in default-assembly and
no-assembly configurations. Every host assembly selection condition now also
requires `BLST_ASM`. Experimental CUDA is conditionally isolated and was not
built. Apple/AArch64 qualification remains covered by the release CI gate, not
this local x64 packet.

### 13. Performance-sensitive structure

All arithmetic selection remains forced-inline and compile-time. Moving the
zero checks from field wrappers to this shared boundary adds no second branch
to field calls. MSVC `/O2 /Ob3` constant-count probes for all four shifts are
direct tail calls to their assembly symbols with no zero branch or copy. The
restored SHA wrapper is a direct tail call to assembly in the production
profile and to `cc_*` in the portable profile. No benchmark was warranted.

### 14. Tests and evidence

Focused native tests passed 6/6: all 256/384 shift cases and deterministic SHA
selection. Focused no-assembly shift tests passed 5/5. Exact test objects were
compiled individually and linked into focused executables without rebuilding
the repository. Strict selector probes compile in both configurations.

### Specialized additions

The selected SHA test uses deterministic bytes rather than `random_device`, so
failures are reproducible. Direct assembly/portable equality and selected-path
code generation cover the dispatch regression separately.

## Findings

- `AUD-041` - validated: scalar assembly/portable configuration selection.
- `AUD-051` - validated: zero-count shifts at every selected boundary.
- `AUD-055` - validated: production SHA compression selection restored.

## Validation

- MSVC 19.50 exact native/default-assembly focused tests: 6/6 passed.
- MSVC 19.50 exact no-assembly focused shift tests: 5/5 passed.
- MSVC 19.50 `/O2 /Ob3 /W4 /WX` selector probes compile for native and
  portable configurations.
- Native constant-count shift probes tail-call `lshift_mod_256`,
  `rshift_mod_256`, `lshift_mod_384`, and `rshift_mod_384` directly.
- Native SHA probe tail-calls `blst_sha256_block_data_order`; portable SHA
  probe tail-calls `cc_blst_sha256_block_data_order`.
- Strict probe object SHA-256 hashes:
  - native: `937454395CCEE25F9551FDA154D27FD00960CE7D65450D2107425DE2478539C4`
  - portable: `6C09655A5BC0C3DE7964F5B3B8B035E41F9B85CA5B45EBC8FE926BA38C053710`
- Structural searches found no duplicate `BLST_HD`, stale interface guard,
  old 32-bit divisor expression, or default-assembly condition lacking the
  assembly-presence requirement.
- `git diff --check` and the audit manifest validator pass at closure.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Phase 12 scalar database and security queries | No additional selector-header finding | Does not model handwritten assembly timing or generated loop contracts |
| Clang Static Analyzer / clang-tidy | Phase 12 compile database and triaged baseline | No selector-header diagnostic | Compile-time configuration coverage required the separate probes |
| Clippy | Retained Rust feature families | N/A for C++ header | Rust reaches this file only through native linkage |
| ASan/UBSan | Phase 12 pure-C scalar baseline | `AUD-041` evidence retained | Not rerun for inline dispatch-only fixes; zero-count behavior has direct tests |
| Fuzzing | No direct selector target | Deferred to field/protocol consumers | Assembly selection and code generation require structural evidence |

## Residual risk and decision

Complete for supported host scalar contracts at reviewed commit
`dc07be802c39a3c04a17347a5cec64f74d95ecd7`. Assembly implementations and
CPU capability initialization retain their own pending/generated-file reviews.
Apple/AArch64 must continue to pass CI before release. Experimental CUDA is
outside this decision and was not used as evidence. Broader field, protocol,
and public API properties remain assigned to later file reviews.
