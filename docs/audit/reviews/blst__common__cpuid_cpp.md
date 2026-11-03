# File Review: `blst/common/cpuid.cpp`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/common/cpuid.cpp` |
| Roadmap phase / queue | Phase 13, common platform-selection foundation |
| Reviewed commit | `2f59883caf67191aae95abb905598f405f8139e8` |
| Reviewer | Codex audit team |
| Date | 2026-08-03 |
| Disposition | `in_scope_production` |
| Subsystem | Host CPU capability discovery for generated scalar/SHA assembly |
| Profiles | scalar assembly, C, C++, Rust; x86-64 and AArch64 host builds |
| Generation/provenance | Hand-maintained C++ adaptation of upstream BLST `src/cpuid.c` |
| Upstream references | Upstream head `54e6e55674722fc2797ebb4bbb71b26d881eb4b8`, `src/cpuid.c` |
| Related files | `blst/CMakeLists.txt`, x86-64 and AArch64 assembly generators and generated outputs |

## Purpose and trust boundary

This translation unit initializes the process-wide capability word consumed by
generated assembly. It reads public CPU/OS feature state only; no secret value
crosses the boundary. Callers rely on exact C symbol linkage, initialization
before normal API calls, and conservative zero bits when a feature cannot be
proven.

## Contracts and invariants

- The only externally visible object is `blst_platform_cap`, with C linkage.
- Project-owned identifiers do not use reserved double-underscore names.
- Query and initialization helpers have translation-unit-local linkage.
- x86-64 bits are enabled only after the corresponding CPUID leaf and feature
  bit are present.
- AArch64 Linux reads `AT_HWCAP`; Apple uses the platform feature policy in
  the retained upstream-derived branch.
- Initialization allocates nothing, throws nothing, and performs no dispatch
  based on secret data.

## Dependency and call-path map

`blst/CMakeLists.txt` compiles this file whenever assembly is enabled. Generated
SHA and arithmetic providers reference `blst_platform_cap`, which causes the
object and its static initializer to be retained from a static library. C,
C++, and Rust reach the selected providers through the same compiled library;
the capability word is not a public API.

## Upstream/specification comparison

The feature masks and OS checks were compared with upstream `src/cpuid.c`.
CatBLST's intentional changes are C++ compilation, explicit `extern "C"` on
the data boundary, internal linkage on helpers, and the non-reserved symbol
name propagated through every generator and checked-in assembly output.

## Checklist evidence

### 1. Identity and purpose

The file, build inclusion, generated consumers, and exact upstream source are
recorded above.

### 2. Contracts and invariants

The C data symbol, internal helpers, conservative feature-bit policy, and
startup lifetime are explicit in source and object symbols.

### 3. Upstream and specification lineage

CPUID/XGETBV and AArch64 HWCAP logic retain upstream semantics. The symbol and
linkage changes are local integration repairs rather than feature changes.

### 4. Correctness

ADX, SHA-NI, and SSSE3 bits follow their corresponding leaves. Unsupported
cases preserve a zero capability word.

### 5. Memory and lifetime safety

The translation unit owns one static integer and retains no external pointer.
CPUID output uses fixed local arrays. There is no allocation or variable-sized
memory access.

### 6. Undefined and implementation-defined behavior

Compiler intrinsics and OS APIs are selected by explicit platform macros.
The former C-to-C++ linkage ambiguity and reserved project identifiers are
removed.

### 7. Side channels and cryptographic hygiene

N/A for secret-dependent behavior: all inputs are public machine features.
The selected providers still require their own constant-time review.

### 8. Randomness and secret lifecycle

N/A: no entropy or secret material.

### 9. API and ABI

The assembly ABI is the exact unmangled object symbol `blst_platform_cap`.
It is internal to the compiled library and unaffected by the future public C
prefix configuration.

### 10. Concurrency

Initialization occurs during static startup, before supported API use. The
word is read-only thereafter, so normal concurrent calls require no lock.

### 11. Resource exhaustion and denial of service

N/A: fixed work, no allocation, and no input-controlled loop.

### 12. Portability and build behavior

MSVC x64 strict compilation and the assembly-enabled CMake target pass.
Supported run `31665323681` builds and tests the Apple Silicon assembly-backed
library and retained platform initializer.

### 13. Performance-sensitive structure

Feature detection runs once. Hot assembly reads one capability word and no
runtime abstraction, heap object, or thread pool is introduced.

### 14. Tests and evidence

The MSVC object exposes the exact C symbol and local helpers. The assembly
target links successfully, and generator verification reproduces all 100
checked-in outputs after the rename.

### Specialized additions

Assembly generators and their outputs use the same symbol. This file has no
FFT, AVX-512 vector-emulation, or Rust-unsafe implementation of its own.

## Findings

`AUD-058` - validated.

## Validation

- MSVC 19.50 `/W4 /WX /O2` compile: passed.
- `dumpbin /symbols` showed external `blst_platform_cap`; helpers were local.
- `cmake --build out/build/x64-Debug --target blst_tests --config Debug
  --parallel 1`: passed with assembly enabled.
- `./scripts/build/refresh.sh verify-assembly`: all 100 outputs passed.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Phase 12 retained scalar database | `AUD-058` found by hand dependency audit | Static analysis does not prove static-library retention or native CPU behavior |
| Clang Static Analyzer / clang-tidy | Prior retained profiles | No additional mapped diagnostic | Native runtime feature-state exhaustiveness is platform-owned |
| Clippy | Retained Rust profiles | N/A for C++ initializer | Rust consumes compiled providers |
| ASan/UBSan | Prior scalar baselines | No mapped report | Does not validate capability-bit selection |
| Fuzzing | N/A | N/A | Fixed public machine state, no input grammar |

## Residual risk and decision

The file review is complete at the recorded source commit. Hosted run
`31665323681` supplies the remaining Apple Silicon assembly-backed build and
execution evidence, so `AUD-058` is validated and this file is complete.
Future supported OS/compiler additions must add a conservative feature probe
rather than assuming availability.

## Batch E SHA capability follow-up

The x86 initializer was reread at `2f59883`. CPUID leaf 1 ECX bit 9 now maps
to capability bit 2 for SSSE3. Existing capability bits 0 and 1 retain ADX and
SHA-NI. Unsupported CPUID leaves leave their capability bits clear.

Forced provider tests demonstrate that the three resulting states select C,
SSSE3, and SHA-NI correctly. Capability discovery remains one-time public
initialization and adds no secret-dependent behavior. `AUD-375` is validated.
