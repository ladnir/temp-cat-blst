# File Review: `blst/CMakeLists.txt`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/CMakeLists.txt` |
| Roadmap phase / queue | Phase 19, Packet 19A product target owner |
| Reviewed commit | `20c075a5000b8f7c29bb06a714112b1c8f0f9742` |
| Reviewer | Codex audit team |
| Date | 2026-08-12 |
| Disposition | `in_scope_production` |
| Subsystem | `core` |
| Profiles | scalar, assembly, AVX-512, experimental CUDA/parallel MSM |
| Generation/provenance | Handwritten local CMake target definition |
| Upstream references | not_applicable: local C++ build architecture |
| Related files | root CMake, `config.h.in`, AVX CMake, production translation units |

## Purpose and trust boundary

This file selects every production translation unit and assembly object,
generates `config.h`, and defines the compile/link usage requirements of
`blst_static`. Its public properties affect embedded consumers.

## Contracts and invariants

- Every production source appears once.
- Implementation optimization, warning, ADX, and sanitizer compilation flags
  do not become ordinary consumer compile options.
- Ordinary native C++ requires C++20; experimental CUDA/parallel MSM requires
  C++23.
- The library is position-independent and uses the caller/toolchain's normal
  optimization and MSVC runtime policy.
- Native selected AVX C++ consumers receive the exact AVX usage target.

## Dependency and call-path map

Root CMake supplies options; this file owns `blst_static`, checked-in assembly,
generated `config.h`, `blst_avx_headers`, sanitizer propagation, and optional
experimental CUDA sources. C/Rust link the archive; native C++ also consumes
its public include and language requirements.

## Upstream/specification comparison

N/A as a file. Source membership is compared through the individual production
file reviews and upstream manifest.

## Checklist evidence

### 1. Identity and purpose

All source lists, conditions, options, target properties, and links inspected.

### 2. Contracts and invariants

Compile usage is separated into private implementation and legitimate public
language/selected-AVX requirements.

### 3. Upstream and specification lineage

Covered by owning source records.

### 4. Correctness

Scalar and assembly Release archives link; the prior duplicate portable source
entry is removed.

### 5. Memory and lifetime safety

N/A at configuration time.

### 6. Undefined and implementation-defined behavior

No implicit host CPU tuning remains; `-fwrapv` is private where retained.

### 7. Side channels and cryptographic hygiene

No runtime dispatch or abstraction was added. Exact native AVX flags are owned
by the AVX usage target.

### 8. Randomness and secret lifecycle

N/A.

### 9. API and ABI

The public C++ language feature is C++20. Install/export remains deferred.

### 10. Concurrency

No build definition adds an internal thread pool.

### 11. Resource exhaustion and denial of service

N/A beyond root target defaults.

### 12. Portability and build behavior

MSVC scalar/assembly and hosted Linux/Apple execution pass. Native AVX remains
residual.

### 13. Performance-sensitive structure

Release optimization is returned to toolchain defaults; explicit arithmetic
semantics and ISA flags are private/selected as required.

### 14. Tests and evidence

Serial static builds and a linked downstream consumer passed.

### Specialized additions

Public transitive compile properties were inspected through the consumer's
actual `compile_commands.json` entry.

## Findings

`AUD-340`, `AUD-343`; related root findings `AUD-341`, `AUD-342`, `AUD-346`.

## Validation

See Phase 19 Packet 19A. The downstream executable linked the archive and did
not inherit `/arch`, warning, optimization, host-native, or `-fwrapv` flags.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CMake/Ninja | MSVC 19.50 plus hosted run `31665323681` | scalar/assembly Linux, Windows, and Apple pass | Native AVX remains; experimental CUDA is outside the production claim |
| Static/sanitizer/fuzz tools | deferred | N/A | Product sources retain their subsystem evidence |

## Residual risk and decision

Status remains `in_review` pending exact native AVX-512/IFMA release-target
execution. Supported Linux scalar/assembly, Windows assembly, and Apple Silicon
are green in run `31665323681`. CUDA is experimental and its C++23 build
remains outside the production claim.
