<div align="left">
  <img src="blst_logo_small.png" alt="category-blst">
</div>

`category-blst` is a substantially modified fork of
[BLST](https://github.com/supranational/blst) for high-performance BLS12-381
cryptography. It retains BLST's scalar and assembly foundation, incorporates an
AVX-512/IFMA implementation derived largely from
[VROOM](https://github.com/SimonLangowski/VROOM), and includes an experimental
CUDA backend.

The library provides low-level field, curve, pairing, serialization, and
signature operations, together with FFT, multi-scalar multiplication, and
middle-product facilities. C, native C++, and Rust are the supported language
surfaces.

## Backends and platforms

| Backend | Status | Intended use |
| --- | --- | --- |
| Scalar and assembly | Supported | General-purpose x86_64 and Apple Silicon deployments |
| Native AVX-512/IFMA | Supported | High-throughput x86_64 deployments with the required ISA |
| CUDA | Experimental | GPU experimentation; excluded from the production support gate |

The supported production platforms are little-endian 64-bit x86_64 and Apple
Silicon macOS. Other ARM64 and Morello/CHERI sources are retained but are not
currently supported platforms. Native AVX-512 binaries require AVX-512IFMA and
the other instruction-set features selected by the build; they must not be
deployed to incompatible CPUs.

The library does not create a global thread pool. Applications control
parallelism, scheduling, resource limits, and randomness.

## BLST compatibility

The project preserves substantial BLST compatibility so existing C and Rust
consumers can be ported without an unnecessary rewrite. Public C symbols,
library artifacts, and the Rust crate currently retain their established
`blst` names. Category Labs extensions and native C++ APIs are not necessarily
available in upstream BLST, and compatibility should not be interpreted as
identical behavior for undocumented or internal interfaces.

## APIs and contracts

- [C umbrella header](bindings/blst.h)
- [Scalar C API](bindings/blst_scalar.h)
- [AVX C API](bindings/blst_avx.h)
- [FFT C API](bindings/blst_fft.h)
- [Middle-product C API](bindings/blst_middle_product.h)
- [Native C++ implementation and headers](blst)
- [Rust binding](bindings/rust)

This is a low-level cryptographic library. Some fast operations deliberately
have narrow preconditions and do not repeat validation already expected from
the caller. The shared terminology for validation, mathematical domains,
aliasing, and timing guarantees is defined in the
[API contract vocabulary](docs/API_CONTRACTS.md). Function declarations and
their associated tests are the authoritative source for individual contracts.

## Examples

Tests and examples:

- [C API tests](test/C_api_tests.cpp)
- [Signature conformance tests](test/Conformance_tests.cpp)
- [Aggregation tests](test/Aggregate_audit_tests.cpp)
- [Serialization and validation tests](test/CompressedInfinity_tests.cpp)
- [Pairing tests](test/Pairing_tests.cpp)
- [G1](test/G1_tests.cpp), [G2](test/G2_tests.cpp), and [GT](test/GT_tests.cpp) tests
- [FFT tests](test/Fft_tests.cpp) and [FFT C API tests](test/Fft_api_tests.cpp)
- [MSM tests](test/Msm_tests.cpp)
- [Middle-product tests](test/MiddleProduct_c_api_tests.cpp)
- [AVX tests](test/avx)
- [Rust conformance tests](bindings/rust/tests/conformance.rs)
- [Rust smoke example](bindings/rust/examples/smoke.rs)

## Building

The native library requires CMake and a C++20 compiler. A basic scalar build is:

```sh
cmake -S . -B out/build/release -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DBLST_BENCHMARK=OFF
cmake --build out/build/release --target blst_static
```

To build and run the main test executable:

```sh
cmake --build out/build/release --target blst_tests
ctest --test-dir out/build/release --output-on-failure
```

Important backend options include:

- `BLST_AVX512=ON` — compile the native AVX-512 backend.
- `BLST_DEFAULT_AVX=ON` — select the native AVX backend as the default
  implementation.
- `BLST_CUDA=ON` — enable the experimental CUDA backend.

Repository presets in [CMakePresets.json](CMakePresets.json) provide common
Windows and Linux configurations. Build artifacts currently retain the
`blst`/`libblst` compatibility names.

For Rust build features, supported targets, and examples, see the
[Rust binding README](bindings/rust/README.md).

## Security

Inputs crossing a trust boundary must be deserialized and validated according
to the selected API contract before they are passed to subgroup-only or other
narrow-domain operations. Upstream BLST security reviews do not automatically
cover Category Labs modifications. Please report suspected vulnerabilities as
described in [SECURITY.md](SECURITY.md).

## Provenance

category-blst builds on work from:

- [BLST](https://github.com/supranational/blst), which provides the original
  library, scalar arithmetic, and architecture-specific assembly foundation.
- [VROOM](https://github.com/SimonLangowski/VROOM), from which much of the
  AVX-512 design and implementation lineage is derived.
- Category Labs contributors, who maintain the combined library and its
  extensions, hardening, tests, and platform integration.

Detailed retained AVX source pins and adaptation boundaries are recorded in
the [AVX documentation](docs/avx/README.md).

## License

category-blst is licensed under the
[Apache License 2.0](LICENSE). Retained upstream and third-party material
remains subject to its recorded notices and attribution requirements.
