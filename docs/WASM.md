# WebAssembly scalar backend

The initial browser target is `wasm32-unknown-emscripten`, using the existing
scalar C++ backend. Pointers and `size_t` are 32 bits; arithmetic limbs remain
64 bits. Canonical field, point, and pairing encodings match native builds.
Native assembly and AVX implementations are unchanged.

This restores backend portability. It does not yet provide a high-level
JavaScript wallet or Confidential Monad prover package.

## Build and test

Install the official [Emscripten SDK](https://emscripten.org/docs/getting_started/downloads.html)
and activate version **6.0.13**. The new CI workflow pins the SDK driver source
to `d7b4130688bfa348392b921ab04015aa15c4f154`.
Run these commands from an SDK-activated shell:

```sh
emcmake cmake -S . -B out/build/wasm -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DBLST_TEST=ON -DBLST_BENCHMARK=OFF
cmake --build out/build/wasm \
  --target blst_tests blst_portability_vectors middle_product_workspace_allocation_test \
  --parallel 4
ctest --test-dir out/build/wasm --output-on-failure --parallel 2 -E '^audit\.'
```

CMake selects the portable scalar backend automatically. It rejects assembly,
ADX, AVX, CUDA, parallel MSM, and wasm64 configurations for this target.
Other unsupported 32-bit native targets remain rejected.
BTE's BLAKE3 dependency also selects portable C rather than inferred x86 SIMD.

The static library is `out/build/wasm/blst/libblst.a`. Applications link it with
`em++ -fwasm-exceptions`. Native WebAssembly exception handling preserves the
C API's exception-to-error handling without JavaScript-based exception
instrumentation. The browser must support WebAssembly exception handling.
Do not mix C++ exception ABIs between translation units.

The unit-test executable uses a 4 MiB stack and growable memory. These are
test settings, not library-wide requirements or application memory estimates.
Applications choose and measure their own limits.

This profile is single-threaded. It requires neither pthreads nor shared Wasm
memory. Run expensive client operations in a Web Worker to keep the UI
responsive. Threaded execution and Wasm SIMD are not qualified by this patch.

## Browser execution

The build copies a browser test harness next to the test executable:

```sh
python3 -m http.server 8765 --bind 127.0.0.1 --directory out/build/wasm/test
```

Open `http://127.0.0.1:8765/browser-smoke.html`. The page downloads the generated
JavaScript and Wasm, runs arithmetic/FFT/conformance tests locally, and displays
`PASS` or a failure. It is a test harness, not a wallet; it contains no secrets
and makes no chain RPC requests. The CI job executes this page in headless Chrome.

## Rust consumers

Install the Rust target and configure the C++ linker, explicit CMake toolchain,
and matching exception flags. In an SDK-activated Bash shell:

```sh
rustup target add wasm32-unknown-emscripten
export CMAKE_TOOLCHAIN_FILE="$EMSDK/upstream/emscripten/cmake/Modules/Platform/Emscripten.cmake"
export CARGO_TARGET_WASM32_UNKNOWN_EMSCRIPTEN_LINKER=em++
export CARGO_TARGET_WASM32_UNKNOWN_EMSCRIPTEN_RUNNER=node
export RUSTFLAGS='-C link-arg=-fwasm-exceptions -C link-arg=-sSTACK_SIZE=4194304 -C link-arg=-sALLOW_MEMORY_GROWTH=1'
cargo run --manifest-path bindings/rust/Cargo.toml --locked --release \
  --target wasm32-unknown-emscripten --example smoke
cargo test --manifest-path bindings/rust/Cargo.toml --locked --release \
  --target wasm32-unknown-emscripten --lib
```

On Windows with SDK 6, set the linker to the full path of
`$EMSDK/upstream/emscripten/em++.exe`; the SDK no longer supplies the
`emcc.bat` name assumed by Rust's default Windows target configuration.
Use `em++`, not `emcc`, so libc++'s matching exception runtime is linked.
The binding links libc++, not the native Linux `libstdc++` runtime.

To reuse the CMake-built archive instead of rebuilding it inside Cargo:

```sh
export BLST_PREBUILT_LIB="$PWD/out/build/wasm/blst/libblst.a"
export BLST_PREBUILT_CONFIG_INCLUDE="$PWD/out/build/wasm/blst"
```

The ordinary explicit-toolchain check still applies. The archive and generated
configuration must match the Cargo target. This profile does not currently
support `wasm32-unknown-unknown`, wasm-bindgen, or WASI directly.

## Validation and limits

The local validation used Emscripten 6.0.13, Rust 1.94.1, Node 24, and headless
Chrome on Windows. Native regression tests used MSVC with x86-64 assembly.
The CI workflow additionally specifies an Ubuntu run; that hosted run is not
claimed until it executes.

The tests cover field/group arithmetic, hash-to-curve, signatures, pairings,
prepared FFTs, MSMs, decoding, invalid sizes, overflow rejection, and allocation
failure. `blst_portability_vectors` emits deterministic canonical results for
SHA-256 boundaries, hash-to-curve, pairings, field operations, MSM, and FFT.
Compare native and Wasm executables with:

```sh
node test/wasm/compare-vectors.cjs native-vector-executable out/build/wasm/test/blst_portability_vectors.js
```

The matching vector output has SHA-256:

```text
a90e527c7cab3aeb6cfbac4aca9da38a3ba88f54b3c04f1903db7cd5d1d1a9c9
```

AVX-only tests report skips. The single-threaded profile explicitly skips the
C++ concurrency test and the Rust concurrency tests. Four optional BTx tests
also skip unless their external shared-vector fixture is configured, as on
native builds. Ignored benchmark tests are not executed.

The complete CTest invocation also runs source-bound formal audit checks.
At the base revision, `audit.avx_proof_correspondence` reports a stale hash for
the unchanged `blst/avx/curve/scalar_mult.h`. This port does not refresh or weaken
that ledger. The new runtime CI job excludes `audit.*`; existing audit jobs
retain their own checks.

These tests establish cross-target functional evidence, not a constant-time
claim about browser JITs or a new cryptographic security audit. No benchmark
or browser proving-time claim is included.
