# File Review: `CMakePresets.json`

## Record

| Field | Value |
| --- | --- |
| Path | `CMakePresets.json` |
| Roadmap phase / queue | Phase 19, Packet 19A preset policy |
| Reviewed commit | `20c075a5000b8f7c29bb06a714112b1c8f0f9742` |
| Reviewer | Codex audit team |
| Date | 2026-08-12 |
| Disposition | `in_scope_support` |
| Subsystem | `build_configuration` |
| Profiles | scalar, assembly, ADX/BMI2, AVX-512/IFMA, Apple Silicon |
| Generation/provenance | Handwritten local CMake preset document |
| Upstream references | not_applicable: local supported matrix |
| Related files | root and `blst/` CMake files, CI workflows, README |

## Purpose and trust boundary

The presets are the obvious user and developer entry points. A visible preset
must not imply unsupported deployment, hidden benchmark work, or inactive
machine-local CUDA configuration.

## Contracts and invariants

Visible presets name their deployment requirements. Benchmarks, CUDA, and
portable AVX emulation are absent. Scalar C++, runtime-dispatched assembly,
ADX/BMI2-native, AVX-512/IFMA, and Apple Silicon are distinct choices.

## Dependency and call-path map

Each visible preset inherits a hidden common, scalar, assembly, or Windows-x64
policy and passes only cache values consumed by root CMake.

## Upstream/specification comparison

N/A: this is Category Labs' supported build matrix.

## Checklist evidence

### 1. Identity and purpose

Every preset and inherited cache value was expanded and reviewed.

### 2. Contracts and invariants

CPU-specific presets state required deployment features.

### 3. Upstream and specification lineage

N/A.

### 4. Correctness

`cmake --list-presets` parses and lists exactly ten visible presets.

### 5. Memory and lifetime safety

N/A.

### 6. Undefined and implementation-defined behavior

Explicit target/ISA choices replace implicit host-native behavior.

### 7. Side channels and cryptographic hygiene

N/A beyond selecting the reviewed backend.

### 8. Randomness and secret lifecycle

N/A.

### 9. API and ABI

All presets currently build the compatibility-named static library.

### 10. Concurrency

No parallel build count or internal scheduler is imposed.

### 11. Resource exhaustion and denial of service

Benchmarks default off in the shared hidden base.

### 12. Portability and build behavior

Machine-local CUDA paths and developer portable-AVX presets are removed.

### 13. Performance-sensitive structure

Native ADX/AVX presets remain explicit; general presets avoid host tuning.

### 14. Tests and evidence

Parsing passed locally; hosted Linux/Apple/AVX execution remains pending.

### Specialized additions

Inheritance conflicts and visible/hidden classification were checked through
CMake's parser.

## Findings

`AUD-345`, with related `AUD-340`, `AUD-342`, and `AUD-346`.

## Validation

Local CMake lists all ten intended visible names and no portable-AVX or CUDA
preset. The Windows scalar/assembly policies also configured through equivalent
fresh build directories.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CMake preset parser | local CMake | passed | Host conditions and remote compilers not executed |
| Static/sanitizer/fuzz tools | N/A | N/A | Declarative non-runtime file |

## Residual risk and decision

Status remains `in_review` until each supported hosted preset has configured
and its applicable release target has executed.

