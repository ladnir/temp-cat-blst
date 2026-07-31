# Phase 7C Rust Aggregate and Dependency Evidence

## Record

- Date: 2026-07-31
- Implementation commit: `fed5142`
- State: validated and committed
- Scope: upstream `b83544b` and `dece82e`
- Benchmarks: not run or audited; the benchmark target was compile-checked only
- CUDA: experimental and not exercised
- BTx/BTE: frozen and untouched

## Decision

Port `AggregatePublicKey::sub_aggregate` through the local
`sig_variant_impl!` architecture so both supported high-level variants receive
the same API. The operation copies the supplied projective aggregate, negates
the copy, and adds it to the receiver. It performs no allocation, does not
modify the supplied aggregate, and has the same validation expectations as
`add_aggregate`.

Remove the unnecessary direct `rand` development dependency. Tests and the
compile-only benchmark target import `RngCore` and `SeedableRng` from the
already-resolved `rand_core 0.6.4`; `rand_chacha` continues to provide the
deterministic generator. This introduces no new resolved dependency version.

## Validation

Following `VALIDATION_POLICY.md`, validation stayed within the Rust consumer
surface:

| Surface | Focused evidence | Result |
| --- | --- | --- |
| `min_pk` aggregate subtraction | Sum four keys, subtract the final three, compare with the first key, preserve the subtrahend, and subtract self to infinity | Passed |
| `min_sig` aggregate subtraction | Same macro-generated regression over G2 public keys | Passed |
| Native FFI link | Focused Rust tests rebuilt and linked the native library | Passed |
| Dependency graph | `cargo tree` dev graph and inverse `rand` lookup | Direct and transitive `rand` package absent; `rand_core 0.6.4` retained |
| Benchmark source compatibility | `cargo check --benches --features full-rust-bindings` | Passed; no benchmark executed |
| Formatting/integrity | Cargo formatting and Git whitespace | Passed |

No C/C++, AVX-512, sanitizer, CUDA, full Rust suite, or benchmark execution was
warranted. The new behavior is isolated to the feature-gated Rust wrapper and
is exercised for both macro instantiations.

## Residual validation

- An older Rust compiler was not invoked. The replacement dependency and
  version were already present in the resolved graph, so this packet does not
  raise the dependency-version floor.
- `Cargo.lock` is intentionally not versioned for this library crate; evidence
  uses the resolved Cargo graph rather than a lockfile diff.
- `sub_aggregate` does not group-check its input, matching `add_aggregate` and
  the existing aggregate API contract.
