[![CI (current remote)](https://github.com/ladnir/blst/actions/workflows/ci.yml/badge.svg)](https://github.com/ladnir/blst/actions/workflows/ci.yml) [![CodeQL (current remote)](https://github.com/ladnir/blst/actions/workflows/codeql-analysis.yml/badge.svg)](https://github.com/ladnir/blst/actions/workflows/codeql-analysis.yml)

# category-blst

`category-blst` is Category Labs' performance- and assurance-focused BLS12-381
cryptographic library. It is maintained at
[category-labs/category-blst](https://github.com/category-labs/category-blst)
and is built from the BLST code lineage.

The short project name is **CatBlst** in prose and **CATBLST** in uppercase
identifiers. Existing `blst_*` C symbols, `blst` build artifacts, and the
[`blst`](https://crates.io/crates/blst) Rust crate name are currently retained
for upstream compatibility. Renaming those surfaces is a separate, deliberate
compatibility phase; new documentation should call the library
`category-blst`, not a BLST fork.

## Table of Contents

  * [Status](#status)
  * [Naming and compatibility](#naming-and-compatibility)
  * [General notes on implementation](#general-notes-on-implementation)
  * [Platform and Language Compatibility](#platform-and-language-compatibility)
  * [API](#api)
  * [Introductory Tutorial](#introductory-tutorial)
    + [Public Keys and Signatures](#public-keys-and-signatures)
    + [Signature Verification](#signature-verification)
    + [Signature Aggregation](#signature-aggregation)
    + [Serialization Format](#serialization-format)
  * [Build](#build)
    + [C and C++ library](#c-and-c-library)
  * [Language-specific notes](#language-specific-notes)
    + [Rust](#rust)
  * [Repository Structure](#repository-structure)
  * [Performance](#performance)
  * [License](#license)

## Status

**category-blst is under active development and is not yet release-qualified.**

The in-depth production-readiness audit, hardening sequence, and current
evidence are tracked under [`docs/audit`](docs/audit), starting with the
[audit status](docs/audit/STATUS.md) and
[hardening plan](docs/audit/HARDENING_PLAN.md).

NCC Group audited the upstream BLST codebase in January 2021; the report is
available [here](https://research.nccgroup.com/wp-content/uploads/2021/01/NCC_Group_EthereumFoundation_ETHF002_Report_2021-01-20_v1.0.pdf).
That report predates and does not cover Category Labs' changes.

Galois' formal-verification work on the upstream BLST foundation is available
[here](https://github.com/GaloisInc/BLST-Verification). Category Labs maintains
additional Lean proofs for the selected AVX construction under
[`proof/avx`](proof/avx); neither effort should be interpreted as a proof of
the entire current library.

The basic signature ciphersuites use the domains and hash-to-curve behavior
specified by
[BLS Signatures, draft-irtf-cfrg-bls-signature-07](https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-bls-signature-07)
and
[RFC 9380, Hashing to Elliptic Curves](https://www.rfc-editor.org/rfc/rfc9380.html).
Focused C++ and Rust vectors cover both G1-signature and G2-signature
orientations. The version-pinned BLS document is an active Internet-Draft and
work in progress, not an IETF standard.

The serialization formatting is implemented according to [the ZCash definition](#serialization-format).

## Naming and compatibility

- **category-blst** is the project, library, and repository name.
- **CatBlst** is the short human-readable/internal project name.
- **CATBLST** is the spelling reserved for new uppercase configuration and
  preprocessor identifiers.
- Existing **blst**, `blst_*`, and `BLST_*` names remain where compatibility or
  staged migration requires them. Their presence does not change the project
  name.

During migration, hosted CI may still run from the legacy
[`ladnir/blst`](https://github.com/ladnir/blst) repository; the badges above
intentionally point there. The Category Labs repository is canonical for the
project and will take over CI when the migration is complete.

## General notes on implementation

category-blst is intended as a low-level component for applications and
libraries that require high-performance BLS12-381 operations with explicit,
reviewable contracts. Scalar and native AVX-512 implementations are the primary
production paths. The project combines line-by-line review, differential and
conformance tests, automated analysis, and targeted formal verification; each
evidence source has documented boundaries.

The library deliberately abstains from dealing with memory management and
multi-threading. The C, native C++, and Rust APIs never create threads or own a
scheduler; each operation runs on its calling thread. Applications retain
control of parallel decomposition, resource limits, and oversubscription.
Another responsibility left to the application is random number generation.
All this is in the name of run-time neutrality, which makes integration into
more stringent environments like Intel SGX or ARM TrustZone straightforward.

## Platform and Language Compatibility

The supported hardware targets are x86_64 and Apple Silicon macOS. CatBlst
requires a little-endian target and rejects other byte orders at compile time.
Other ARM64 and Morello/CHERI sources are retained because they may work, but
they have not completed the CatBlst release-validation gate and are not fully
supported. They are disabled by default rather than silently shipping as
supported targets.

Apple Silicon macOS needs no experimental override. To attempt another ARM64
CMake build, set
`BLST_ALLOW_UNTESTED_ARM64=ON`. A CHERI pure-capability build additionally
requires `BLST_ALLOW_UNTESTED_CHERI=ON`. These options acknowledge an
unsupported build; they do not confer production-readiness or ABI guarantees.
Rust builds use the same option names as environment variables with a true
value such as `1` or `ON`.

Native AVX-512 is a production x86_64 backend. `BLST_AVX512_PORTABLE` is
emulation for differential tests and analysis, not a deployment fallback or
native performance evidence. It is omitted from the ordinary Release presets
and production CI matrix, but remains available through explicit configuration
when developers need optimized emulation. Production systems without the
required AVX-512/IFMA target use the scalar backend.

The supported language surfaces are the [C API](bindings/blst.h), the native
[C++ implementation and API](blst), and the [Rust binding](bindings/rust).

## API

The main scalar C API is defined in [bindings/blst.h](bindings/blst.h) together with the split headers under [bindings/](bindings). The API can be categorized as follows, with some example operations:
- Field Operations (add, sub, mul, neg, inv, to/from Montgomery)
- Curve Operations (add, double, mul, to/from affine, group check)
- Intermediate (hash to curve, pairing, serdes)
- BLS12-381 signature (sign, verify, aggregate)

Public preconditions, postconditions, mathematical domains, validation, and
timing guarantees use the common terminology in
[API Contract Vocabulary](docs/API_CONTRACTS.md). CatBlst permits fast expert
operations with narrow contracts; callers must establish documented
preconditions that the operation does not explicitly validate.

The split headers are:

- [bindings/blst_types.h](bindings/blst_types.h)
  - shared POD C types and layout definitions
- [bindings/blst_scalar.h](bindings/blst_scalar.h)
  - scalar/public C API
- [bindings/blst_avx.h](bindings/blst_avx.h)
  - expert AVX C API
- [bindings/blst_middle_product.h](bindings/blst_middle_product.h)
  - middle-product C API
- [bindings/blst_fft.h](bindings/blst_fft.h)
  - Fr FFT, prepared-domain, and isolated interpolation C API
- [bindings/blst_aux.h](bindings/blst_aux.h)
  - compatibility umbrella that includes the split surfaces

## Introductory Tutorial

Programming is understanding, and understanding implies mastering the lingo. So we have a pair of additive groups being mapped to multiplicative one... What does it mean? Well, this tutorial is not about explaining that, but rather about making the connection between what you're supposed to know about [pairing-based cryptography](https://en.wikipedia.org/wiki/Pairing-based_cryptography) and the interface provided by the library.

### Public Keys and Signatures

We have two elliptic curves, E1 and E2, points on which are contained in `blst_p1` and `blst_p2`, or `blst_p1_affine` and `blst_p2_affine` structures. Elements in the multiplicative group are held in a `blst_fp12` structure. One of the curves, or more specifically, a subset of points that form a cyclic group, is chosen for public keys, and another, for signatures. The choice is denoted by the subroutines' suffixes, `_pk_in_g1` or `_pk_in_g2`. The most common choice appears to be the former, that is, `blst_p1` for public keys, and `blst_p2` for signatures. But it all starts with a secret key...

The secret key is held in a 256-bit `blst_scalar` structure which can be instantiated with either [`blst_keygen`](https://tools.ietf.org/html/draft-irtf-cfrg-bls-signature#section-2.3), or deserialized with `blst_scalar_from_bendian` or `blst_scalar_from_lendian` from a previously serialized byte sequence. It shouldn't come as surprise that there are two uses for a secret key:

- generating the associated public key, either with `blst_sk_to_pk_in_g1` or `blst_sk_to_pk_in_g2`;
- performing a sign operation, either with `blst_sign_pk_in_g1` or `blst_sign_pk_in_g2`;

As for signing, unlike what your intuition might suggest, `blst_sign_*` doesn't sign a message, but rather a point on the corresponding elliptic curve. You can obtain this point from a message by calling `blst_hash_to_g2` or `blst_encode_to_g2` (see the [IETF hash-to-curve](https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve#section-3) draft for distinction). Another counter-intuitive aspect is the apparent g1 vs. g2 naming mismatch, in the sense that `blst_sign_pk_in_g1` accepts output from `blst_hash_to_g2`, and `blst_sign_pk_in_g2` accepts output from `blst_hash_to_g1`. This is because, as you should recall, public keys and signatures come from complementary groups.

Now that you have a public key and signature, as points on corresponding elliptic curves, you can serialize them with `blst_p1_serialize`/`blst_p1_compress` and `blst_p2_serialize`/`blst_p2_compress` and send the resulting byte sequences over the network for deserialization/uncompression and verification.

### Signature Verification

Even though there are "single-shot" `blst_core_verify_pk_in_g1` and `blst_core_verify_pk_in_g2`, you should really familiarize yourself with the more generalized pairing interface. `blst_pairing` is an opaque structure, and the only thing you know about it is `blst_pairing_sizeof`, which is how much memory you're supposed to allocate for it. In order to verify an aggregated signature for a set of public keys and messages, or just one[!], you would:
```
blst_pairing_init(ctx, hash_or_encode, domain_separation_tag);
blst_pairing_aggregate_pk_in_g1(ctx, PK[0], aggregated_signature, message[0]);
blst_pairing_aggregate_pk_in_g1(ctx, PK[1], NULL, message[1]);
...
blst_pairing_commit(ctx);
result = blst_pairing_finalverify(ctx, NULL);
```
**The essential point to note** is that it's the caller's responsibility to ensure that public keys are group-checked with `blst_p1_affine_in_g1`. This is because it's a relatively expensive operation and it's naturally assumed that the application would cache the check's outcome. Signatures are group-checked internally. Not shown in the pseudo-code snippet above, but `aggregate` and `commit` calls return `BLST_ERROR` denoting success or failure in performing the operation. Call to `finalverify`, on the other hand, returns boolean.

Another, potentially more useful usage pattern is:
```
blst_p2_affine_in_g2(signature);
blst_aggregated_in_g2(gtsig, signature);
blst_pairing_init(ctx, hash_or_encode, domain_separation_tag);
blst_pairing_aggregate_pk_in_g1(ctx, PK[0], NULL, message[0]);
blst_pairing_aggregate_pk_in_g1(ctx, PK[1], NULL, message[1]);
...
blst_pairing_commit(ctx);
result = blst_pairing_finalverify(ctx, gtsig);
```
What is useful about it is that `aggregated_signature` can be handled in a separate thread. And while we are at it, aggregate calls can also be executed in different threads. This naturally implies that each thread will operate on its own `blst_pairing` context, which will have to be combined with `blst_pairing_merge` as threads join.

### Signature Aggregation

Aggregation is a trivial operation of performing point additions, with `blst_p2_add_or_double_affine` or `blst_p1_add_or_double_affine`. Note that the accumulator is a non-affine point.

---

That's about what you need to know to get started with nitty-gritty of actual function declarations.

### Serialization Format

From the ZCash BLS12-381 specification

* Fq elements are encoded in big-endian form. They occupy 48 bytes in this form.
* Fq2 elements are encoded in big-endian form, meaning that the Fq2 element c0 + c1 * u is represented by the Fq element c1 followed by the Fq element c0. This means Fq2 elements occupy 96 bytes in this form.
* The group G1 uses Fq elements for coordinates. The group G2 uses Fq2 elements for coordinates.
* G1 and G2 elements can be encoded in uncompressed form (the x-coordinate followed by the y-coordinate) or in compressed form (just the x-coordinate). G1 elements occupy 96 bytes in uncompressed form, and 48 bytes in compressed form. G2 elements occupy 192 bytes in uncompressed form, and 96 bytes in compressed form.

The most-significant three bits of a G1 or G2 encoding should be masked away before the coordinate(s) are interpreted. These bits are used to unambiguously represent the underlying element:

* The most significant bit, when set, indicates that the point is in compressed form. Otherwise, the point is in uncompressed form.
* The second-most significant bit indicates that the point is at infinity. If this bit is set, the remaining bits of the group element's encoding should be set to zero.
* The third-most significant bit is set if (and only if) this point is in compressed form _and_ it is not the point at infinity _and_ its y-coordinate is the lexicographically largest of the two associated with the encoded x-coordinate.

## Build

The supported native build uses CMake and a C++20 compiler. Configure a preset,
build the required target, and run only the appropriate test scope. For example:

```sh
# Windows scalar Debug
cmake --preset x64-Debug-Scalar
cmake --build out/build/x64-Debug-Scalar --target blst_static
ctest --test-dir out/build/x64-Debug-Scalar

# Linux/WSL Release
cmake --preset linux
cmake --build out/build/linux --target blst_static
ctest --test-dir out/build/linux
```

Build products currently retain the compatibility name `blst`/`libblst`.
CMake options and generated configuration determine the selected scalar,
assembly, or AVX backend; do not assume that a binary built for one x86_64 CPU
is safe to deploy on another. Native AVX-512 requires the documented AVX-512
and IFMA target. `BLST_AVX512_PORTABLE` is developer emulation, not a production
fallback.

### C and C++ library

The `blst_static` target builds the native static library. Public C headers are
under [`bindings/`](bindings), while the supported native C++ implementation
and headers are under [`blst/`](blst). The C ABI remains `blst_*` during the
compatibility transition described above.

## Language-specific notes

### [Rust](bindings/rust)

[`blst`](https://crates.io/crates/blst) is the current compatibility name of
the Rust binding crate. Its supported build policy and feature contracts are
documented in the Rust binding [README](bindings/rust/README.md).

Low-level/raw use:
```
use blst::{blst_fp12, blst_p1, blst_p2};
```

To use the middle-product API:
```
use blst::{combine, middle_product, precompute_middle_product_b, update};
```

To use the high-level min-pk version:
```
use blst::min_pk::*;
```

To use the high-level min-sig version:
```
use blst::min_sig::*;
```

Useful Rust commands:
```
cargo test --lib
cargo test --lib --features full-rust-bindings
cargo run --example smoke
```

## Repository Structure

**Root** - Contains various configuration files, documentation, licensing, and a build script
* **Bindings** - Contains the files that define the compatibility C and Rust interfaces
    * blst.h - provides C API to blst library
    * blst_types.h - provides shared POD C types used across the split headers
    * blst_scalar.h - provides the scalar/public C API
    * blst_avx.h - provides the expert AVX C API
    * blst_middle_product.h - provides the middle-product C API
    * blst_aux.h - compatibility umbrella that includes the split C surfaces
    * **Rust** - folder containing the Rust crate, tests, benchmarks, and examples
    * **Vectors**
        * **Hash_to_curve**: folder containing test for hash_to_curve from IETF specification
* **blst** - contains the native C++ implementation and supported C++ headers for fields, curves, pairing, hashing, and the C ABI bridge
    * **Asm** - folder containing Perl scripts that are used to generate assembly code for different hardware platforms including x86 with ADX instructions, x86 without ADX instructions, and ARMv8, and [ABI](https://en.wikipedia.org/wiki/Application_binary_interface)[1]
* **Build** - this folder containing a set of pre-generated assembly files for a variety of operating systems and maintenance scripts.
    * **Cheri** - assembly code for use on [CHERI](https://www.cl.cam.ac.uk/research/security/ctsrd/cheri/) platforms
    * **Coff** - assembly code for use on Windows systems with GNU or LLVM toolchain
    * **Elf** - assembly code for use on Unix systems
    * **Mach-o** - assembly code for use on Apple operating systems
    * **Win64** - assembly code for use on Windows systems with Microsoft toolchain

[1]: See [refresh.sh](build/refresh.sh) for usage. This method allows for simple reuse of optimized assembly across various platforms with minimal effort.

## Performance
The [Rust binding](bindings/rust) provides benchmarks for a variety of signature-related operations.

## License

category-blst is licensed under the
[Apache License Version 2.0](LICENSE) software license. Retained upstream code
and third-party material remain subject to their recorded notices and
attributions.
