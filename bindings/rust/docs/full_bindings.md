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

## Validation model

The wrapper types do not record whether a point has been validated. A
`PublicKey`, `Signature`, or aggregate may therefore contain raw unvalidated
storage, a decoded on-curve point, or a fully validated subgroup element. Safe
Rust prevents memory unsafety; it does not silently add curve or subgroup
checks.

- `deserialize`, `uncompress`, `from_bytes`, and Serde decoding establish a
  canonical on-curve representation. They do not establish subgroup
  membership.
- `PublicKey::key_validate` and `Signature::sig_validate` decode and establish
  the subgroup precondition. Public-key validation also rejects infinity;
  signature validation does so when requested.
- A `*_groupcheck` or `*_validate` argument set to `false` means the caller has
  already established the documented precondition. It does not mean that an
  unvalidated input is safe for the operation.
- Raw-to-wrapper `From` conversions perform no mathematical validation.

The shared definitions of on-curve, subgroup, validation, canonical scalar,
and variable-time are in `docs/API_CONTRACTS.md`.

## Aggregate verification flow

For inputs received from an untrusted boundary:

1. Decode and validate each public key and signature, or use an explicitly
   approved batch-validation procedure that establishes the same subgroup
   postcondition.
2. Aggregate only validated subgroup elements. Passing `false` to a validation
   flag is appropriate only after this step has already occurred.
3. Call the verification method whose protocol assumptions match the scheme.
   `fast_aggregate_verify` additionally requires a valid proof of possession
   for every public key.
4. For `verify_multiple_aggregate_signatures`, supply independent,
   unpredictable, nonzero random weights with sufficient soundness for the
   protocol. Sparse, repeated, attacker-controlled, or zero weights are not a
   substitute for proper batch-verification randomness.

All high-level pairing and aggregate operations are synchronous on the calling
thread. Methods that own a pairing context or MSM scratch document their
allocation behavior and publish no partial result on validation failure.

## Secret-key invariant

`SecretKey` always contains a canonical nonzero scalar. It deliberately has no
`Default` implementation because the all-zero scalar is invalid. Constructors
either derive that invariant or reject invalid input before publishing a value.
Cloning explicitly duplicates secret material; every owned copy is zeroized on
drop.
