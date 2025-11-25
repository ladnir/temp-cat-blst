# Repository Guidelines

## Project Structure & Module Organization
Core math and protocol logic sits in `blst/`. Numeric back ends are split into `backend/`, `curve/`, and `fields/`, while shared helpers live under `common/`. Command-line entry points and demos stay inside `blst/mains/`. Architecture-specific assembly and generated tables go into `build/`. Tests shadow their source peers inside `test/` (e.g., `test/G1_tests.cpp`, `test/Fp_tests.cpp`). Bindings for Go, Rust, Python, Java, Node.js, C#, and C++ are versioned under `bindings/` with language-specific READMEs. Benchmarks belong in `benchmark/`, and CMake presets emit artifacts beneath `out/`.

## Build, Test, and Development Commands
Configure with the preset that matches your toolchain, such as `cmake --preset x64-Debug` (MSVC), `cmake --preset ninja-clang` (LLVM), or `cmake --preset linux` (WSL). Build outputs via `cmake --build out/build/<preset>` to produce `blst_static` and sample targets. Run the suite using `ctest --test-dir out/build/<preset>` or execute `outuild\<preset>lst_tests.exe` directly. Pass `-D BLST_CUDA=ON` during configuration when compiling GPU kernels and mirror that flag in subsequent builds.

## Coding Style & Naming Conventions
Use four-space indentation; keep function braces on their own line and prefer tight control-flow blocks. Public C/C++ APIs retain the `blst_` prefix in lower_snake_case. Macros, compile flags, and constants are uppercase (e.g., `BLST_PORTABLE`). Favor fixed-width types (`uint64_t`, `limb_t`) and ensure headers compile standalone. Follow each binding folder's idioms without breaking the shared structure.

## Testing Guidelines
GoogleTest provides unit coverage. Add suites next to the closest implementation file and name cases `TEST(Module, behavior)`. GPU checks must be wrapped in `#ifdef BLST_CUDA` with CPU equivalents. Reuse deterministic vectors from `bindings/vectors/` and document new fixtures in the relevant README. Always run the preset-appropriate `ctest` job before submitting.

## Commit & Pull Request Guidelines
Write concise, imperative commit subjects (`bench fix`, `linux cuda compile`) and expand in the body when behavior or portability shifts. Pull requests should link issues, call out toolchains tested, and attach logs for benchmarks or updated bindings. Highlight any new ADX/cache/CUDA requirements and confirm the matching `ctest` run.

## Security & Configuration Tips
Avoid long-lived secrets; keep private material on the stack or zeroize buffers promptly. When opting out of `BLST_PORTABLE`, list the instruction set you target so reviewers can assess portability. CUDA contributions must note driver/runtime versions and include fallbacks for non-GPU builds.
