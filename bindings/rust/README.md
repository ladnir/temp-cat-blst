# category-blst Rust binding

The `blst` crate provides category-blst's Rust interface. The package and crate
names remain `blst` for compatibility with existing source-tree consumers.
The package version is `0.3.15+cat.1`: the `0.3.15` compatibility component
allows one `links = "blst"` package to satisfy established Cargo dependents,
and the build metadata identifies Category's hardened API variant.

New Category consumers should alias the package as `catblst` in Cargo and use
that namespace in source:

```toml
catblst = { package = "blst", version = "=0.3.15" }
```

The package name remains transitional so existing `links = "blst"` dependents
resolve one native library. The intended long-term Rust namespace is
`catblst`.

## Build
[bindgen](https://github.com/rust-lang/rust-bindgen) generates FFI bindings from
the C headers under `bindings/`. The local [build.rs](build.rs) then invokes
CMake inside Cargo's target area and links the resulting static library. On
Linux, set `CC` and `CXX` to select the C and C++ compilers.

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
  - on Windows this maps to the developer portable-AVX implementation
  - portable AVX is for differential testing and development, not deployment
    or native performance evidence; optimized use requires an explicit feature
    selection and is not part of the ordinary release matrix
  - on Linux native AVX requires Clang 21 or newer unless the native CMake build is explicitly configured to allow an unsupported compiler
- `portable`
  - enables the portable scalar build
- `force-adx`
  - forces ADX on supported x86_64 targets
- `no-threads`
  - deprecated no-op retained for older manifests; the production Rust API
    does not create worker threads regardless of this feature

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

The crate's minimum supported Rust version is 1.78. Cargo enforces this through
the package's `rust-version` field, and CI checks the locked supported feature
surface with that exact toolchain.

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
- `src/scalar.rs`
  - lean, allocation-free, `blstrs`-inspired types over the scalar backend
- `src/avx.rs`
  - native AVX representations; these are intentionally distinct from scalar
    curve types
- `src/middle_product.rs`
  - `precompute(b)`, `middle_product(a, precomputed_b, out)`,
    `update(a0, precomputed_b, out)`, and `combine(mid, cnst, out)`
- `src/pairing.rs`
  - owned pairing context wrapper
- `src/pippenger.rs`
  - multi-point helpers used by the full high-level API

Additional narrative docs live in `docs/`:

- `docs/overview.md`
- `docs/scalar.md`
- `docs/middle_product.md`
- `docs/full_bindings.md`

## Usage
There are four primary ways to use the crate.

Lean scalar curve use:
```rust
use catblst::Scalar;
use catblst::scalar::{G1Projective, G2Projective, Gt};

let exponent = Scalar::from(7u64);
let point = G1Projective::generator() * exponent;
let target = Gt::generator() * exponent;
```

This surface is inspired by `blstrs`, but is not a compatibility clone. It does
not implement the `ff`, `group`, or `pairing` ecosystems, create worker threads,
or hide heap allocation. Its wrappers contain category-blst's native scalar
representations and expose explicit raw borrows for kernel interoperation.

Low-level/raw use:
```rust
use catblst::{blst_fp12, blst_p1, blst_p2};
```

Middle-product use:
```rust
use catblst::{combine, middle_product, precompute_middle_product_b, update};
```

High-level signature use:

For minimal-pubkey-size operations:
```rust
use catblst::min_pk::*;
```

For minimal-signature-size operations:
```rust
use catblst::min_sig::*;
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
use catblst::min_pk::SecretKey;

let mut rng = rand::thread_rng();
let mut ikm = [0u8; 32];
rng.fill_bytes(&mut ikm);

let sk = SecretKey::key_gen(&ikm, &[]).unwrap();
let pk = sk.sk_to_pk();

let dst = b"BLS_SIG_BLS12381G2_XMD:SHA-256_SSWU_RO_NUL_";
let msg = b"blst is such a blast";
let sig = sk.sign(msg, dst, &[]);

let err = sig.verify(true, msg, dst, &[], &pk, true);
assert_eq!(err, catblst::BLST_ERROR::BLST_SUCCESS);
```

See:

- `examples/smoke.rs` for a runnable smoke example
- `src/lib.rs` and the submodules under `src/` for unit-test examples
- `benches/blst_benches.rs` for benchmark-style usage
