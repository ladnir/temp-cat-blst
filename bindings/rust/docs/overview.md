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
- `src/scalar.rs`
  - lean, `blstrs`-inspired scalar-field, G1, G2, Miller-loop, and GT types
  - transparent access to the native scalar ABI without representation bridges
- `src/avx.rs`
  - distinct native AVX representations, operations, and explicit workspaces
- `src/pairing.rs`
  - owned pairing context and uniqueness helpers
- `src/pippenger.rs`
  - multi-point helpers for the full high-level binding surface

High-level policy:

- The scalar ABI is the stable public Rust boundary.
- The `scalar` and `avx` modules name distinct backend representations. They do
  not implicitly convert points or GT elements between backends.
- `Scalar` is shared at the crate root because both backends consume the same
  scalar field. Backend kernels may borrow its native `blst_fr` directly.
- AVX exposes both its expert/native C surface and a lean safe typed layer. The
  typed layer wraps the distinct `blst_avx_*` representations; it does not
  reinterpret `blst_p1`, `blst_p2`, or `blst_fp12`.
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

Documentation remains incremental with the hand audit. Field arithmetic,
reviewed G1/G2/GT operations, the high-level signature and aggregate surface,
pairing state, decoding, FFT, and middle product now state their reviewed
preconditions and postconditions. Unreviewed raw declarations inherit the C
contract and remain pending until their owning files are reviewed.

The exact lean scalar type surface and its deliberate omissions are documented
in [`scalar.md`](scalar.md).
The typed native AVX surface, timing rules, and reusable-workspace contracts are
documented in [`avx.md`](avx.md).

Recommended entry points:

- Lean scalar curve API:
  - import `Scalar` from the crate root
  - import `G1Affine`, `G1Projective`, `G2Affine`, `G2Projective`, `Gt`, and
    `MillerLoopResult` from `catblst::scalar`
  - the API is inspired by the useful typed surface of `blstrs`, but deliberately
    omits its `ff`, `group`, and `pairing` trait integration
  - operations run synchronously on the calling thread and allocate no heap
    storage
- Lean native AVX curve API:
  - import backend-specific `G1Projective`, `G2Projective`, and `Gt` from
    `catblst::avx`
  - allocate prepared-line, batched-pairing, MSM, and middle-product workspaces
    before entering a hot loop
  - retain native values until canonical serialization; do not convert through
    scalar point or GT storage
- Low-level/raw use:
  - import concrete `blst_*` POD types from the crate root
  - call raw FFI functions or the small helper methods on `blst_fp12`
- Middle product:
  - `precompute_middle_product_b`
  - `middle_product`
  - `update`
  - `combine`
  - use the stateless scalar or native-AVX range methods when the application
    owns parallel scheduling; [`middle_product.md`](middle_product.md) defines
    the required phase order and barriers
- FFT and interpolation:
  - use the explicit `catblst::fft` module
  - prefer prepared domains for repeated transforms
  - provide caller-owned scratch for FFT-domain suffix interpolation
- High-level signature API:
  - enable `full-rust-bindings`
  - use `min_pk` or `min_sig`
