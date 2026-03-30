# Full Rust Bindings

The `full-rust-bindings` feature enables the broader ergonomic Rust API on top
of the scalar C ABI.

Main pieces:

- `Pairing`
- `uniq`
- `min_pk`
- `min_sig`
- multi-point helpers in `pippenger.rs`

Use this surface when you want:

- `SecretKey`, `PublicKey`, and `Signature`
- aggregate verification helpers
- high-level serialization and validation routines

Use the default surface when you want:

- only raw `blst_*` POD types
- the C ABI directly
- the middle-product API without the rest of the signature wrapper layer
