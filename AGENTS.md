# Repository Guidelines

## Project Structure & Module Organization
Core library code lives in `blst/` with arithmetic back ends under `backend/`, `curve/`, and `fields/`; shared utilities are in `common/`. Entry points and demos stay in `blst/mains/`. Platform-tuned assembly lives in `build/`, while tests mirror the module layout inside `test/` (e.g., `G1_tests.cpp`, `Fp_tests.cpp`). Language bindings for Go, Rust, Python, Java, Node.js, C#, and C++ sit under `bindings/` with their own README files. Benchmarks reside in `benchmark/`, and generated artifacts from CMake presets land under `out/`.

## Build, Test, and Development Commands
Configure with the preset matching your toolchain, e.g. `cmake --preset x64-Debug` for MSVC, `cmake --preset ninja-clang` for LLVM, or `cmake --preset linux` inside WSL. Build with `cmake --build out/build/<preset>` to produce `blst_static` and any examples. Run tests using `ctest --test-dir out/build/<preset>` or by calling `out\build\<preset>\blst_tests.exe`. Add `-D BLST_CUDA=ON` when configuring if a CUDA toolkit is installed, and mirror that setting when developing GPU kernels.

## Coding Style & Naming Conventions
Follow the existing C/C++ layout: four-space indentation, braces on their own line for functions, and compact blocks for control flow. Public symbols keep the `blst_` prefix in lower_snake_case; macros, build flags, and constants use uppercase (e.g., `BLST_PORTABLE`). Prefer fixed-width types (`uint64_t`, `limb_t`) and keep headers self-contained. Bindings should mimic the structure of established language folders while respecting idiomatic naming.

## Testing Guidelines
GoogleTest powers unit coverage. Place new suites beside the closest peer file and name tests with intent-revealing `TEST(Module, behavior)` cases. Guard GPU-specific checks behind `#ifdef BLST_CUDA` and ensure a CPU analogue exists. Use deterministic vectors from `bindings/vectors/` when expanding coverage, and document any new fixtures in the relevant README.

## Commit & Pull Request Guidelines
Commits follow concise, imperative subjects ("bench fix", "linux cuda compile"); elaborate in the body when behavior or portability changes. Pull requests should link issues, summarize platform matrices tested, and attach logs for benchmarks or binding updates. Note any new cache, ADX, or CUDA requirements and confirm the relevant `ctest` runs.

## Security & Configuration Tips
Minimize long-lived secrets: keep private material on the stack or securely zeroized buffers. When disabling `BLST_PORTABLE`, state the target CPU features so reviewers can gauge risk. For CUDA work, capture the tested driver/runtime versions and include fallbacks or guards for non-GPU builds.

