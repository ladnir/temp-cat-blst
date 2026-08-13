# File Review: `test/avx512_compat_test.cpp`

## Record

| Field | Value |
| --- | --- |
| Path | `test/avx512_compat_test.cpp` |
| Roadmap phase / queue | Phase 16 selected-vector boundary follow-up |
| Reviewed commit | `9e237456a946abe827aa55a8a122e4a040114ffc` |
| Reviewer | Codex audit team |
| Date | 2026-08-13 |
| Disposition | `in_scope_support` |
| Subsystem | Native-versus-portable AVX-512 differential tests |
| Profiles | combined native AVX-512/IFMA plus portable AVX |
| Generation/provenance | Hand-maintained Category Labs test |
| Upstream references | `not_applicable: CatBLST-specific compatibility test` |
| Related files | `blst/backend/avx512/intrin_compat.h`, `blst/backend/avx512/intrin_compat.cpp`, `test/CMakeLists.txt` |

## Purpose and trust boundary

This test compares each retained pure C/C++ AVX emulation operation with the
corresponding native Intel intrinsic. It must run only when both implementations
are enabled and on AVX-512/IFMA hardware. Inputs are deterministic test values,
not secrets.

## Contracts and invariants

- Configuration fails at preprocessing unless native and portable support are
  both enabled.
- Static assertions require `vv` to select native intrinsic types while `cc`
  remains an independent lane-array representation.
- Explicit bridges are confined to this differential boundary.
- The test executable must discover a nonzero test count.
- Random cases use a fixed seed and are reproducible.

The file exposes no C, native C++, or Rust API.

## Dependency and call-path map

GoogleTest invokes ten cases. They exercise construction, arithmetic, exact
comparison, shifts, unaligned transfers, mask blending, variable permutation,
extraction, low-128-bit operations, broadcast, and both IFMA halves. The file
does not link or exercise the complete product library.

## Upstream/specification comparison

N/A: this is a CatBLST regression and differential test. Native Intel
intrinsics are the executable reference for the emulated operation semantics.

## Checklist evidence

### 1. Identity and purpose

The file's isolated compatibility-test role is explicit.

### 2. Contracts and invariants

Configuration, type separation, deterministic input, and direct-comparison
requirements are enforced in source.

### 3. Upstream and specification lineage

N/A beyond Intel intrinsic semantics; no upstream BLST counterpart exists.

### 4. Correctness

Boundary and deterministic random inputs compare every retained operation.
IFMA tests include 52-bit truncation, carries, accumulator wrap, and random
full-width lanes.

### 5. Memory and lifetime safety

Aligned arrays, deliberately unaligned byte buffers, and exact 64-byte
transfers cover both load and store boundaries.

### 6. Undefined and implementation-defined behavior

High-bit unsigned literals are compared through `bit_cast`. Shift cases include
0, 63, 64, 65, and 255. Extraction indices are valid compile-time constants.

### 7. Side channels and cryptographic hygiene

N/A: this file tests functional correspondence. It does not establish
constant-time behavior; optimized production probes provide that evidence.

### 8. Randomness and secret lifecycle

The pseudorandom generator has a fixed public seed and creates no secret.

### 9. API and ABI

N/A: test-only translation unit.

### 10. Concurrency

Each fixture owns its generator. Tests use no shared mutable state.

### 11. Resource exhaustion and denial of service

Every loop has a fixed small bound and every object has fixed storage.

### 12. Portability and build behavior

The target is intentionally hardware-specific and combined-profile-only. It
passed on Peach GCC 15.2 with native AVX-512/IFMA.

### 13. Performance-sensitive structure

N/A: no benchmark or production implementation. Its limited target avoids a
portable-AVX product rebuild.

### 14. Tests and evidence

All 10 discovered GoogleTest cases passed. The CMake target rejects an empty
test discovery result.

### Specialized additions

The reviewed source SHA-256 is
`aaaff8a8a6778d58772e52a02cd66680da559eb8bb83d201055f93d09131a985`.

## Findings

`AUD-036`, `AUD-037`, `AUD-385`, and `AUD-386`.

## Validation

Peach GCC 15.2 Release built this focused target with one job and ran its ten
tests successfully. Local MSVC 19.50 compiled the combined-profile object.
No broad suite, benchmark, CUDA, BTx/BTE audit, sanitizer, Apple CI, or
portable AVX Release build ran.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Existing C/C++ baseline | No test-source finding | Runtime ISA availability and intrinsic lowering need hardware evidence |
| Clang Static Analyzer / clang-tidy | Existing curated baseline | No additional finding | Differential semantic coverage is supplied by execution |
| Clippy | N/A | N/A | No Rust source |
| ASan/UBSan | Not rerun | Deferred milestone evidence | Focused fixed-buffer cases passed without sanitizer instrumentation |
| Fuzzing | N/A | N/A | Deterministic random and boundary vectors cover this finite adapter surface |

## Residual risk and decision

The test covers the complete retained emulation surface on one supported native
GCC host. Other supported compiler profiles remain release-milestone evidence.
That limitation does not leave this support-file review open.
