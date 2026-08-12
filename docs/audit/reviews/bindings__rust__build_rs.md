# File Review: `bindings/rust/build.rs`

## Record

| Field | Value |
| --- | --- |
| Path | `bindings/rust/build.rs` |
| Roadmap phase / queue | Phase 18 Rust build and packaging bridge |
| Reviewed commit | `906fc2c6964fcff3f8d188b5d373fd25334aa828` |
| Reviewer | Codex audit team |
| Date | 2026-08-12 |
| Disposition | in-scope production |
| Subsystem | Rust-to-native build routing and static linkage |
| Profiles | Rust; scalar; native AVX-512; Windows portable AVX developer mode |
| Generation/provenance | Handwritten Rust build script with upstream BLST lineage |
| Upstream references | Common base `f48500c1fdbefa7c0bf9800bccd65d28236799c1`; local fail-closed routing intentionally diverges |
| Related files | `Cargo.toml`, `README.md`, root `CMakeLists.txt`, `blst/config.h.in`, CI workflows |

## Purpose and trust boundary

The script translates Cargo target, profile, and feature metadata into one
native CMake build and publishes the resulting archive/configuration metadata
to Rust and dependent build scripts. Environment variables selecting a
toolchain or prebuilt archive cross the trust boundary: they are explicit
caller assertions about target compatibility, not autodetected facts.

## Contracts and invariants

- `OPT_LEVEL=0` selects native Debug; every optimized Cargo level selects
  native Release independently of the profile's name.
- `HOST != TARGET` requires an existing explicit `CMAKE_TOOLCHAIN_FILE` before
  CMake starts. Native Apple Silicon (`HOST == TARGET`) remains ordinary.
- `avx` is x86-64 only. On Windows it means portable AVX developer emulation;
  elsewhere it means native AVX. Optimized Windows emulation is explicit but is
  not deployment or native-performance evidence.
- `portable`, `avx`, and `force-adx` combinations fail closed where their
  semantics conflict. Host CPUID is never used to infer ADX for a cross build.
- No archive is discovered implicitly. A prebuilt requires both
  `BLST_PREBUILT_LIB` and the directory directly containing its `config.h`;
  common policy, runtime linkage, rebuild watches, and `DEP_BLST_*` metadata
  still apply.
- Native outputs use profile-specific deterministic directories and exact
  candidate paths. The script does not recursively delete build trees.
- CUDA, benchmarks, BTx/BTE, and the internal parallel pipeline are disabled.
  The script creates no threads; CMake and Cargo retain process scheduling.

## Checklist evidence

1. Identity/purpose: one Cargo build-script entry point and pure policy helpers.
2. Contracts: target/profile/feature/prebuilt behavior is documented above and
   in the package README.
3. Lineage: upstream convenience behavior was narrowed where it could silently
   build or select the wrong native artifact.
4. Correctness: policy helpers cover optimized levels, conflicts, architecture,
   Windows portable AVX, and cross-toolchain presence.
5-6. Memory/UB: ordinary owned `String`/`PathBuf`/`Command` values; no unsafe
   Rust or FFI.
7-8. Side channels/randomness: N/A; build-time public metadata only.
9. API/ABI: emits one static archive, platform C++ runtime linkage, include
   metadata, and checked cfg declarations.
10. Concurrency: no background work or scheduler.
11. Resources: deterministic directories replace recursive search/removal;
    native compilation remains the caller/CI resource boundary.
12. Portability: target architecture/OS/environment are separated from host;
    unsupported ARM64/CHERI policy remains fail closed.
13. Performance: optimized Rust profiles cannot accidentally link Debug C++.
14. Tests: four warning-clean build-policy unit tests plus focused explicit-
    prebuilt and no-std package checks.

## Findings

`AUD-326`, `AUD-327`, `AUD-328`, `AUD-329`, `AUD-330`, `AUD-333`.

## Validation

- `rustc --test build.rs -D warnings`: 4/4 passed.
- Explicit `BLST_PREBUILT_LIB`/`BLST_PREBUILT_CONFIG_INCLUDE` full-binding
  focused Rust link/test passed.
- Portable AVX is absent from ordinary release presets/CI and retained in the
  path-scoped Debug developer lane; explicit optimized selection remains legal.
- `BLST_TEST_NO_STD=1 cargo check --features serde-secret` passed with warnings
  denied.

## Automated evidence

| Tool/profile | Version and configuration | Result | Blind spots |
| --- | --- | --- | --- |
| rustc policy tests | local Rust toolchain, `-D warnings` | 4/4 passed | Does not execute every external CMake generator |
| Cargo focused link | MSVC cached scalar archive, full bindings | passed | Explicit prebuilt path, not a fresh native build |
| Clippy/CodeQL | Phase 12 Rust baselines | no controlling diagnostic beyond listed findings | Build-script process behavior still needs platform CI |
| Sanitizers/fuzzing | N/A | N/A | Build routing is covered by deterministic policy tests rather than data fuzzing |

## Residual risk and decision

Fresh source builds on Linux, Windows, and Apple remain ordinary CI evidence at
the candidate commit; a user-supplied cross toolchain or prebuilt archive can
still be wrong because the caller explicitly owns that assertion. The script
now refuses every previously silent path and is complete at the reviewed
commit.
