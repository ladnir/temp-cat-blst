# File Review: `CMakeLists.txt`

## Record

| Field | Value |
| --- | --- |
| Path | `CMakeLists.txt` |
| Roadmap phase / queue | Phase 19, Packet 19A root configuration owner |
| Reviewed commit | `20c075a5000b8f7c29bb06a714112b1c8f0f9742` |
| Reviewer | Codex audit team |
| Date | 2026-08-12 |
| Disposition | `in_scope_support` |
| Subsystem | `build_configuration` |
| Profiles | scalar, assembly, AVX-512, Apple policy, experimental CUDA policy |
| Generation/provenance | Handwritten local CMake project |
| Upstream references | not_applicable: upstream BLST has a different build system |
| Related files | `CMakePresets.json`, `blst/CMakeLists.txt`, `blst/config.h.in`, `test/CMakeLists.txt`, `bindings/rust/build.rs` |

## Purpose and trust boundary

This file admits targets, compilers, architectures, backend combinations,
sanitizers, tests, and experimental features before any product source is
compiled. Its cache is an input from users, CI, Rust, and parent projects; it
must reject configurations whose generated `config.h` or instruction set would
misrepresent the resulting binary.

## Contracts and invariants

- Production targets are 64-bit x86-64 and Apple Silicon; other ARM64 and
  CHERI require their explicit experimental acknowledgements.
- ADX/BMI2 and native AVX-512 are x86-64-only deployment requirements.
- Assembly/default-assembly and CUDA/CUDA-MSM choices must be consistent.
- Shared and removed VROOM-vendor modes fail rather than do nothing.
- `BLST_FAST_BUILD` is valid only with experimental CUDA.
- Tests default on only at the top level; benchmarks default off.
- CUDA-off configuration must not mutate CUDA toolchain cache variables.

## Dependency and call-path map

The root creates the sanitizer usage target, configures backend policy, then
adds `blst/`, frozen `BTE/`, optional `test/`, and excluded `benchmark/`.
`bindings/rust/build.rs` invokes it with tests, benchmarks, CUDA, and frozen
experimental paths disabled. `blst/config.h.in` receives the selected feature
values.

## Upstream/specification comparison

No one-to-one upstream CMake owner exists. Supported platform policy comes
from Phases 10 and 12-18; compiler/ISA policy comes from the selected scalar
and MatrixNoK backend reviews.

## Checklist evidence

### 1. Identity and purpose

All cache options, target additions, and status output were reviewed.

### 2. Contracts and invariants

Invalid option relationships now fail before target generation.

### 3. Upstream and specification lineage

N/A beyond the local policy sources above.

### 4. Correctness

Fresh scalar and assembly configurations generate; six invalid configurations
fail as intended.

### 5. Memory and lifetime safety

N/A: configuration code owns no runtime object.

### 6. Undefined and implementation-defined behavior

Architecture and pointer-width gates prevent unsupported generated layouts and
assembly selection.

### 7. Side channels and cryptographic hygiene

ISA selection is explicit; runtime arithmetic review remains with its owners.

### 8. Randomness and secret lifecycle

N/A.

### 9. API and ABI

The file generates the feature header. Shared/export packaging remains open.

### 10. Concurrency

No scheduler or thread pool is introduced.

### 11. Resource exhaustion and denial of service

Benchmarks and embedded tests no longer expand an ordinary build by default.

### 12. Portability and build behavior

The complete architecture, compiler, sanitizer, CUDA-off, and option paths
were inspected. Native cross-platform execution remains residual evidence.

### 13. Performance-sensitive structure

Implicit host-native tuning is removed; explicit native configurations retain
their required ISA without runtime abstraction.

### 14. Tests and evidence

See `docs/audit/PHASE_19.md`, Packet 19A.

### Specialized additions

CMake cache mutation, parent-project defaults, and try-compile propagation were
reviewed explicitly.

## Findings

`AUD-340` through `AUD-343`, `AUD-346`; `AUD-016` remains open for exports.

## Validation

MSVC 19.50 scalar and assembly Release configuration/build, invalid-option
configurations, embedded consumer configuration/link, preset parsing, and
compile-database inspection passed as recorded in Phase 19.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CMake | 4.2.0-compatible parser, project minimum 3.27 | passed; `AUD-340`-`AUD-346` | Linux/Apple/native AVX hosted execution pending |
| Static/sanitizer/fuzz tools | N/A | N/A | Build-policy file has no runtime input parser |

## Residual risk and decision

Status remains `in_review` until Linux, Apple, native AVX, and negative
cross-toolchain evidence execute. CUDA remains experimental, and final
install/shared/export behavior is deferred to `AUD-016` and Phase 20.

