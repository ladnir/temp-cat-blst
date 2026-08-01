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
