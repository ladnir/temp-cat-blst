# Rust Bindings Overview

The Rust crate is split by surface area rather than keeping everything in
`src/lib.rs`.

Current structure:

- `src/lib.rs`
  - crate entry point
  - raw generated bindings import
  - top-level re-exports
- `src/fft.rs`
  - explicit module for Fr FFT, prepared domains, general interpolation, and
    caller-scratch FFT-domain suffix interpolation
- `src/middle_product.rs`
  - safe `precompute(b)`, `middle_product(a, precomputed_b, out)`,
    `update(a0, precomputed_b, out)`, and `combine(mid, cnst, out)` wrappers
- `src/primitives.rs`
  - low-level scalar-ABI helpers for points, Fp12, and scalar hashing
- `src/pairing.rs`
  - owned pairing context and uniqueness helpers
- `src/pippenger.rs`
  - multi-point helpers for the full high-level binding surface

High-level policy:

- The scalar ABI is the stable public Rust boundary.
- AVX remains an expert/native C surface and is not reinterpreted as a different
  meaning for `blst_p1`, `blst_p2`, or `blst_fp12`.
- Middle product is intentionally isolated so it can evolve without dragging the
  rest of the crate with it.
- FFT and interpolation are intentionally module-scoped rather than re-exported
  from the crate root.

Contract policy:

- Safe Rust methods document Rust-specific panics/errors and retain the same
  mathematical preconditions as the C operation they wrap. Memory-safe to call
  does not mean curve, subgroup, GT, or protocol validity was checked.
- Raw functions from `src/bindings.rs` are generated unsafe FFI declarations.
  Their semantic contract is the matching declaration in
  `bindings/blst_scalar.h`, together with Rust's pointer, alignment, and
  initialization requirements.
- The generated file is not edited to add independent promises. Its generation
  and inherited documentation will be reviewed during the Rust binding phase.
- Shared terminology such as canonical scalar, G1/G2 element, GT element,
  validation, and variable-time is defined in `docs/API_CONTRACTS.md` at the
  repository root.

Documentation is currently incremental with the hand audit. Field arithmetic
and the reviewed G1/G2/GT multiplication and membership families are covered;
unaudited signature, pairing-state, decoding, FFT, and middle-product families
retain their existing documentation until their files are reviewed.

Recommended entry points:

- Low-level/raw use:
  - import concrete `blst_*` POD types from the crate root
  - call raw FFI functions or the small helper methods on `blst_fp12`
- Middle product:
  - `precompute_middle_product_b`
  - `middle_product`
  - `update`
  - `combine`
- FFT and interpolation:
  - use the explicit `blst::fft` module
  - prefer prepared domains for repeated transforms
  - provide caller-owned scratch for FFT-domain suffix interpolation
- High-level signature API:
  - enable `full-rust-bindings`
  - use `min_pk` or `min_sig`
