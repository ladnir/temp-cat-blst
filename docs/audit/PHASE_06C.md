# Phase 6C Unsupported Language Binding Removal Evidence

## Record

- Date: 2026-07-31
- Implementation commit: `fd0e698c08e3d9ec0a8d46bc4788c147a3bd3fb9`
- State: validated and committed
- Scope: SWIG, C#, Emscripten JavaScript, Java, Node.js, Python, and `AUD-020`
- Benchmarks: excluded and not run
- CUDA: disabled and not exercised

## Decision

CatBLST supports C, native C++, and Rust. The SWIG, C#, Emscripten JavaScript,
Java, Node.js, and Python surfaces were proof-of-concept or legacy bindings
outside that support policy. The project has no external users or compatibility
requirement for them, so Phase 6C removes them rather than auditing, repairing,
packaging, or publishing them.

The binding-only upstream rows remain `not_applicable`. The removed wrappers
did not define retained native C++ behavior, and their language/runtime
compatibility changes have no supported consumer.

## Removed surface

The implementation commit deletes 21 tracked files:

- `bindings/blst.swg`;
- four C# files;
- four Emscripten JavaScript files;
- four Java files;
- seven Node.js/TypeScript files; and
- the Python proof-of-concept driver.

It also removes their manifest-classification rules and public repository
documentation. Shared C headers, Rust bindgen/header-trimming tools, and shared
test vectors remain because they are independent retained inputs.

The production key-generation comment that named a historical SWIG Python
caller was rewritten to state the actual invariant: validate a nullable `info`
pointer before trusting its caller-supplied length. The defensive behavior is
unchanged.

## Validation

| Surface | Configuration | Result |
| --- | --- | --- |
| C and native C++ | Windows MSVC Debug, direct ADX scalar, post-rebuild | 282 registered: 273 passed, nine expected AVX-only skips, zero failed |
| AVX-512 | Windows MSVC Debug, portable AVX-512 default | `blst_tests` rebuilt successfully |
| Rust | `cargo test --features full-rust-bindings` | 71 unit tests and one compressed-infinity integration test passed |
| Rust packaging | `cargo package --list --allow-dirty` | No removed binding or SWIG path is packaged |
| Repository | Filename and content reference searches | No active non-audit SWIG or deleted-binding reference remains |
| Formatting/integrity | Cargo formatting and Git whitespace | Passed |

## Residual validation

- Removed language runtimes and SWIG were not executed because they are outside
  the supported production surface.
- At Phase 6C close, the Rust README still mentioned Emscripten as an optional
  compiler for a Rust WebAssembly target. The later Phase 12 scope decision
  removed WASI/WebAssembly from the supported matrix and removed that guidance.
- No ARM64 or WebAssembly build was run for this removal-only batch.
- `blst/sign/keygen.cpp` received only a comment clarification and remains
  pending its normal full file review.
