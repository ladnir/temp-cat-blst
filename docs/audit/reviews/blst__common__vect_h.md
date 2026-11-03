# File Review: `blst/common/vect.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/common/vect.h` |
| Roadmap phase / queue | Phase 13, common vector boundary after assembly interface |
| Reviewed commit | `80ae0f9a23da1c30415498241e204297e15157af` |
| Reviewer | Codex audit team |
| Date | 2026-08-03 |
| Disposition | `in_scope_production` |
| Subsystem | Fixed-width copy, zero, select, equality, and assembly predicate boundary |
| Profiles | scalar C++, scalar assembly, C/Rust compiled consumers |
| Generation/provenance | Hand-maintained C++ header with upstream BLST `src/vect.h` lineage |
| Upstream references | Upstream head `54e6e55674722fc2797ebb4bbb71b26d881eb4b8`, `src/vect.h` and `src/vect.c` |
| Related files | `blst/common/defines.h`, `blst/asm/add_mod_384-x86_64.pl`, generated vector assembly, `test/Interface_tests.cpp` |

## Purpose and trust boundary

This header is the lowest common byte/vector boundary used by field and curve
code. Values may be secret. Callers require fixed work for a fixed public
length, exact self-copy support, bit-valued low-level selection, compiler-
resistant zeroing, and correct empty predicate behavior.

## Contracts and invariants

- `vec_copy` accepts exact self-alias or disjoint ranges; partial overlap is
  outside the contract and asserted/assumed.
- Known self-copy is a no-op. The distinct host path uses `memcpy`; the device
  fallback retains an explicit byte loop.
- `vec_select`/`vec_czero`/`vec_cswap` receive a bit. Public native C++ callers
  normalize arbitrary `int` truth values before this boundary.
- Zero-length zero/equality/copy operations do not dereference inputs.
- Fixed-size assembly declarations match the generated ABI.

## Dependency and call-path map

Scalar fields, extension fields, points, pairing, signature code, and the
native compatibility templates include this header. Host assembly supplies
fixed select and predicate functions when enabled; portable code supplies the
fallback behavior. C and Rust reach these operations through compiled APIs.

## Upstream/specification comparison

The byte-loop and fixed assembly boundary follow upstream vector helpers.
CatBLST intentionally adds namespace/type integration, defined exact self-copy,
Debug bit-contract assertions, empty predicate guards, and 64-bit generated
predicate counts.

## Checklist evidence

### 1. Identity and purpose

All declarations, inline operations, assembly generators, and supported
callers were inventoried in the common-foundation and assembly-interface waves.

### 2. Contracts and invariants

Aliasing, length, truth-bit, output, and partial-overlap contracts are stated
above and encoded at the boundary.

### 3. Upstream and specification lineage

Compared with upstream `src/vect.h/.c`; local divergence closes actual C++
aliasing and cross-provider truth-value inconsistencies.

### 4. Correctness

Exact self-copy returns before the restricted distinct path. Public `int`
selection maps `0` to false and every other value to true. Generated predicates
now preserve the full `size_t` count on x86-64.

### 5. Memory and lifetime safety

No storage is retained or allocated. Exact alias and disjoint ranges are safe.
Partial overlap is deliberately unsupported rather than silently assigned a
direction-dependent meaning.

### 6. Undefined and implementation-defined behavior

The reachable restrict/self-alias violation is removed. Restricted access is
established only after equality and overlap preconditions. Compiler-specific
barriers remain isolated in `defines.h` evidence.

### 7. Side channels and cryptographic hygiene

Fixed-size select and copy operations have no secret-dependent address. The
selection bit may be secret; Debug assertions disappear in Release. Provider
assembly remains subject to architecture-specific constant-time evidence.

### 8. Randomness and secret lifecycle

No randomness. `vec_zero` uses volatile stores plus the supported compiler
barrier; lifecycle correctness of each caller is reviewed separately.

### 9. API and ABI

The fixed assembly declarations use C linkage and exact byte widths. This is an
internal header, but its behavior underpins every public language surface.

### 10. Concurrency

All operations are stateless and reentrant; callers own every byte accessed.

### 11. Resource exhaustion and denial of service

Work is linear in a public `size_t` length. No allocation occurs. Higher-level
callers remain responsible for bounding attacker-controlled work.

### 12. Portability and build behavior

MSVC host and portable-AVX Debug compilation pass. The CUDA device loop remains
experimental and outside release evidence. Native Apple provider execution is
still tracked by the assembly review.

### 13. Performance-sensitive structure

MSVC `/O2 /Ob3` emits a single `ret` for known self-copy. A disjoint 48-byte
copy emits one equality branch and three vector loads/stores, with no overlap
versioning, helper call, allocation, or runtime dispatch.

### 14. Tests and evidence

Exact self-copy and non-bit public selection regressions pass. Empty predicates
and fixed native calls were validated in the prior assembly-interface packet.

### Specialized additions

The x86 generator/output pair uses 64-bit predicate counts and reproduces
deterministically. AVX-512 `cc`/`vv` operations are reviewed in their own files.

## Findings

`AUD-056`, `AUD-059`, `AUD-060`, and `AUD-061` - validated.

## Validation

- Focused scalar tests: 2/2 passed.
- MSVC 19.50 `/O2 /Ob3 /W4 /WX` copy probes: self=`ret`; distinct=one
  equality branch plus three vector moves.
- Generated Win64/ELF predicates use 64-bit count operations.
- `./scripts/build/refresh.sh verify-assembly`: 100/100 passed.
- Single-worker MSVC Debug target build: passed.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Phase 12 scalar/portable databases | No remaining mapped issue | Restrict legality and constant-time need source/codegen review |
| Clang Static Analyzer / clang-tidy | Retained profiles | No additional mapped diagnostic | Cannot prove generated provider equivalence |
| Clippy | N/A | N/A | Rust consumes compiled ABI |
| ASan/UBSan | Prior scalar baselines | No mapped report | Exact self-alias UB required source-contract repair |
| Fuzzing | Deferred to field/curve packets | N/A | Finite helper cases covered directly; caller composition remains |

## Residual risk and decision

Complete for the common vector contract at the reviewed commit. Partial overlap
remains intentionally unsupported. Apple assembly execution and experimental
CUDA device compilation are separate evidence gates and are not implied by
this closure.
