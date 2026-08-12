# blst [![Crates.io](https://img.shields.io/crates/v/blst.svg)](https://crates.io/crates/blst)

The `blst` crate provides a rust interface to the blst BLS12-381 signature library.

## Build
[bindgen](https://github.com/rust-lang/rust-bindgen) is used to generate FFI bindings from the C headers under `bindings/`. Then [build.rs](https://github.com/supranational/blst/blob/master/bindings/rust/build.rs) invokes CMake to build the native blst library inside Cargo's target area and links the resulting static library into the crate. On Linux it is still possible to choose the compiler by setting `CC` and `CXX`.

Everything can be built and run with the typical cargo commands:

```
cargo test --lib
cargo run --example smoke
cargo bench
```

Relevant Cargo features:

- `full-rust-bindings`
  - enables the higher-level Rust API such as `min_pk`, `min_sig`, aggregation helpers, and pairing helpers
- `avx`
  - enables the native AVX-512 build on x86-64 Unix targets
  - on Windows this maps only to the developer portable-AVX implementation;
    it is accepted in Debug and rejected in optimized builds
  - portable AVX is for differential testing and development, not deployment
  - on Linux native AVX requires Clang 21 or newer unless the native CMake build is explicitly configured to allow an unsupported compiler
- `portable`
  - enables the portable scalar build
- `force-adx`
  - forces ADX on supported x86_64 targets

Apple Silicon macOS is supported without an experimental override. Other
ARM64 and Morello/CHERI builds are retained experimental but are not fully
supported or release-qualified. They fail by default even when a prebuilt
native archive is present. Set `BLST_ALLOW_UNTESTED_ARM64=1` to attempt another
ARM64 Cargo build. A CHERI pure-capability build additionally requires
`BLST_ALLOW_UNTESTED_CHERI=1`. These environment variables are explicit risk
acknowledgements, not support or ABI guarantees.

Examples:

```
cargo test --lib
cargo test --lib --features full-rust-bindings
cargo run --example smoke
cargo run --example smoke --features full-rust-bindings
```

## Threading contract

The Rust crate does not create threads, discover a CPU count, own a thread
pool, or queue background work. Each operation runs to completion on the
calling thread. This keeps scheduling, resource limits, and oversubscription
under application control and gives the Rust API the same runtime-neutral
contract as the C and C++ APIs.

Applications can parallelize independent work with their own runtime. Pairing
work can be split across separately owned `Pairing` contexts and combined with
`Pairing::merge`. Multi-Miller-loop work can be split into paired slices and
the partial Fp12 products multiplied before final exponentiation. MSM work can
be split by matching point/scalar ranges and the resulting group elements
reduced by addition.

The crate's minimum supported Rust version is 1.65. Cargo enforces this through
the package's `rust-version` field.

Cross-target native builds require an explicit `CMAKE_TOOLCHAIN_FILE`. This is
deliberately fail-closed: target metadata alone is not enough to make CMake
produce an archive for that target. Native Apple Silicon builds remain fully
supported and need no override or toolchain file when Cargo itself runs on the
Apple Silicon host.

Prebuilt archives are never discovered implicitly. Troubleshooting builds may
set both `BLST_PREBUILT_LIB` and `BLST_PREBUILT_CONFIG_INCLUDE`; the caller is
responsible for ensuring that the archive was built for the exact Cargo target
and requested feature policy. The include directory must be the generated
directory that directly contains `config.h`.

## Layout
The Rust crate is split by surface area:

- `src/lib.rs`
  - crate entry point and top-level re-exports
- `src/primitives.rs`
  - thin helpers around the scalar C ABI
- `src/middle_product.rs`
  - `precompute(b)`, `middle_product(a, precomputed_b, out)`,
    `update(a0, precomputed_b, out)`, and `combine(mid, cnst, out)`
- `src/pairing.rs`
  - owned pairing context wrapper
- `src/pippenger.rs`
  - multi-point helpers used by the full high-level API

Additional narrative docs live in `docs/`:

- `docs/overview.md`
- `docs/middle_product.md`
- `docs/full_bindings.md`

## Usage
There are three primary ways to use the crate.

Low-level/raw use:
```rust
use blst::{blst_fp12, blst_p1, blst_p2};
```

Middle-product use:
```rust
use blst::{combine, middle_product, precompute_middle_product_b, update};
```

High-level signature use:

For minimal-pubkey-size operations:
```rust
use blst::min_pk::*;
```

For minimal-signature-size operations:
```rust
use blst::min_sig::*;
```

There are five structs with inherent implementations that provide the BLS12-381 signature functionality.
```
SecretKey
PublicKey
AggregatePublicKey
Signature
AggregateSignature
```

A simple example for generating a key, signing a message, and verifying the message:
```rust
use blst::min_pk::SecretKey;

let mut rng = rand::thread_rng();
let mut ikm = [0u8; 32];
rng.fill_bytes(&mut ikm);

let sk = SecretKey::key_gen(&ikm, &[]).unwrap();
let pk = sk.sk_to_pk();

let dst = b"BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_NUL_";
let msg = b"blst is such a blast";
let sig = sk.sign(msg, dst, &[]);

let err = sig.verify(true, msg, dst, &[], &pk, true);
assert_eq!(err, blst::BLST_ERROR::BLST_SUCCESS);
```

See:

- `examples/smoke.rs` for a runnable smoke example
- `src/lib.rs` and the submodules under `src/` for unit-test examples
- `benches/blst_benches.rs` for benchmark-style usage
