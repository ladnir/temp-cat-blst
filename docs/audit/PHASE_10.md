# Phase 10 Experimental ARM64 and CHERI Build Gates

## Record

- Date: 2026-07-31
- Initial implementation commit: `0a94602`
- Apple Silicon correction: `bc93888`
- State: complete
- Scope: CMake, generated configuration, Rust native build orchestration, and
  platform-support documentation
- Benchmarks: not run or audited
- CUDA: experimental and untouched
- BTx/BTE: frozen and untouched

## Decision

Retaining non-Apple ARM64 and Morello/CHERI code must not imply production
support before those targets have compiled and executed under their real
toolchains. They are therefore blocked by default and described as
experimental targets that may work. Their source and generated artifacts
remain available for future qualification.

The initial gate also blocked Apple Silicon macOS. That was broader than the
project's intended platform scope and was corrected immediately in `bc93888`.
Apple Silicon is supported without an override, while its real macOS release
matrix remains outstanding evidence rather than an assumed pass.

An explicit override is preferable to deletion because it permits controlled
bring-up without allowing ordinary CMake, Cargo, package, or CI builds to ship
the targets silently.

## Implementation

- CMake detects ARM64 through the system processor, Apple architecture list,
  generator platform, or compiler architecture ID.
- CMake detects CHERI pure-capability mode with the compiler-defined
  `__CHERI_PURE_CAPABILITY__` macro.
- Apple Silicon macOS proceeds without an experimental override.
- Other ARM64 targets fail configuration unless
  `BLST_ALLOW_UNTESTED_ARM64=ON`.
- CHERI requires both the ARM64 acknowledgement and
  `BLST_ALLOW_UNTESTED_CHERI=ON`.
- Explicitly permitted configurations emit a prominent unsupported-target
  warning.
- Generated `config.h` enforces the same policy at compilation.
- Rust checks the target before its prebuilt-library shortcut, preventing a
  supplied archive from bypassing the policy. It forwards explicitly parsed
  environment acknowledgements to CMake.
- Root and Rust documentation state that the targets may work but are not
  fully supported, release-qualified, or ABI-guaranteed.

## Focused validation

| Surface | Evidence | Result |
| --- | --- | --- |
| Normal x86-64 CMake | Clean configure with both acknowledgements absent | Passed with no experimental-target warning |
| Apple Silicon policy | Simulated Darwin ARM64 portable configuration with no override | Passed the policy gate and generated the Apple-support definition; no Apple object was compiled |
| Apple Rust policy | `cargo check --target aarch64-apple-darwin --lib` without opt-in | Reached the post-policy path; the unavailable target toolchain then prevented compilation |
| Non-Apple ARM64 default | Simulated Linux AArch64 target | Configuration rejected with the required policy message |
| Non-Apple ARM64 override | Same target with `BLST_ALLOW_UNTESTED_ARM64=ON` | Configuration passed and emitted the unsupported-target warning |
| CHERI default | Simulated AArch64 compiler defining `__CHERI_PURE_CAPABILITY__`, ARM64 only acknowledged | Configuration rejected and required both acknowledgements |
| CHERI override | Same compiler with both acknowledgements | Configuration passed and emitted the CHERI warning |
| Generated header | Preprocessor probes for ARM64, CHERI, and both explicit overrides | Both default compilations rejected; double opt-in passed |
| Rust default ARM64 | `cargo check --target aarch64-pc-windows-msvc --lib` without opt-in | Rejected by `build.rs` with the ARM64 policy error |
| Normal Rust | `cargo check --lib` | Passed without warnings |
| Rust formatting | `cargo fmt --check` | Passed |
| Repository formatting | `git diff --check` | Passed |

No Apple, other ARM64, or CHERI object was compiled or executed. This phase
controls claims and prevents accidental unsupported builds; it does not close
the target-validation gaps recorded in `PHASE_09.md` or the supported Apple
release matrix. No full suite or benchmark was warranted because production
code and runtime behavior did not change.
