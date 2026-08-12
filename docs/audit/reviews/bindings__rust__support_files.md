# File Review: Rust Tests, Example, and Publication Helper

## Record

| Field | Value |
| --- | --- |
| Paths | `bindings/rust/src/pippenger-test_mod.rs`; `bindings/rust/examples/smoke.rs`; `bindings/rust/tests/compressed_infinity.rs`; `bindings/rust/tests/conformance.rs`; `bindings/rust/publish.sh` |
| Roadmap phase / queue | Phase 18 Rust support-file closure |
| Reviewed commits | `16e5b4926e56821989a73e67ffc2e33832ffc262`; `5f6c42353192597d676c11e8ae2d5a211342fb09`; `51b23a5387fbbf4bfd0873da77ddf024eda63417` |
| Reviewer | Codex audit team |
| Date | 2026-08-12 |
| Disposition | in-scope support; no benchmark files |
| Profiles | scalar Release Rust; package shell policy |
| Related findings | `AUD-338`, `AUD-339` |

## Purpose and trust boundaries

- `pippenger-test_mod.rs` is included only below `#[cfg(test)]`. It provides
  deterministic G1/G2 oracle, boundary, validation, and alias-contract tests;
  it is not production code.
- `smoke.rs` demonstrates raw subgroup pairing, pairing middle product, and the
  feature-gated high-level signature flow. It is documentation by example and
  must not teach callers to skip validation for adversarial points.
- `compressed_infinity.rs` pins the distinction between successful canonical
  decoding and protocol validation that rejects infinity.
- `conformance.rs` checks the pairing-friendly-curves draft vector, canonical
  generators, and independent noble-curves signature vectors. Its raw field
  construction is test-only and does not expose a decoder.
- `publish.sh` is a maintainer tool with an external crates.io side effect when
  deliberately invoked without `--dry-run`. It must fail before Cargo when the
  source checkout is not the intended clean repository.

## Contracts and invariants

- Pippenger test scalars are contiguous little-endian byte strings. The tests
  cover zero bits, non-byte-aligned masking, 256 bits, short buffers, and widths
  above 256. Every tested point is a generator multiple or infinity.
- Exact C output/input aliasing is exercised with raw pointers, without
  constructing overlapping Rust references. Large deterministic buffers are
  allocated directly as `Vec<u8>` rather than through an array temporary.
- The smoke Miller loop receives subgroup points produced from canonical
  generator multiples and explicitly checks membership before affine use.
- Compressed infinity is accepted by the on-curve decoder and then rejected by
  `validate`; the test does not claim that infinity is a valid protocol key or
  signature.
- The conformance test's draft coefficients are imported in the documented BLST
  tower order. The production pairing is the draft pseudocode value cubed, as
  required by the cited draft's BLS12-381 definition.
- Publication requires a completely clean Git tree, a real symbolic `blst` link
  mapped to the exact repository root, and an isolated target directory unless
  the caller supplies one. The temporary link and target are removed on normal
  exit and HUP/INT/TERM. `--allow-dirty` permits only the helper-created source
  link after the explicit clean-tree check. MSYS must create a native symbolic
  link; copy-style emulation is rejected before Cargo.

## Checklist evidence

### Identity, lineage, and licensing

All five files were read completely and traced to their production owners. The
Pippenger macro retains upstream BLST lineage; the example and integration tests
are local support artifacts; `publish.sh` retains its upstream helper lineage
with local fail-closed hardening. All Rust sources now carry Apache-2.0/SPDX
headers. External conformance vectors identify the pairing draft and exact
noble-curves release/git head in source comments.

### Correctness and memory safety

The Pippenger oracle agrees with repeated single-point multiplication and
complete addition for both groups. Pointer arithmetic remains within the exact
allocated scalar vectors, and every raw output points to initialized owned
storage. The former aliased `&mut`/shared-reference FFI expression is gone.
Integration vectors pass unchanged. No support source accepts an arbitrary
pointer or length from an external caller.

### Side channels and secret handling

MSM tests use deterministic public pseudorandom scalars. The smoke signature
secret is a fixed test-only value. The files neither claim timing evidence nor
handle production secrets. No benchmark or performance conclusion is drawn.

### Failure behavior, resources, and concurrency

Safe Pippenger shape errors panic before FFI, as documented. Test allocations
are bounded constants. The publication helper uses one Cargo process, creates no
threads, fails closed on setup/source-state errors, and confines recursive
cleanup to the exact nonempty path returned by `mktemp -d`.

### Portability and build behavior

Rust sources format under the pinned Rust 1.65 formatter and document under
warning denial. The POSIX shell helper passes syntax checking. Windows scalar
Release tests link an existing reviewed native archive; no CMake/native rebuild
is part of this packet. Native Apple and broad supported-host execution remain
milestone evidence.

## Findings and remediation

- `AUD-338` is validated by the owned carrier, borrow, conversion, mutation,
  subgroup, validation, empty, and allocation documentation. No runtime check or
  production instruction was added.
- `AUD-339` is validated by fail-closed checkout discovery, complete clean-tree
  enforcement before link creation, exact root mapping, real-link enforcement,
  signal cleanup, native POSIX injected-no-op Cargo success, and MSYS safe
  refusal when native link permission is unavailable. No publication command
  contacted a registry.
- Missing support-file SPDX headers, test-only manifest classification, direct
  heap construction, and scalar-width regressions are recorded as review
  remediation rather than separate findings.

## Validation

- Pinned bindgen 0.65.1 and Rust 1.65 rustfmt regeneration: byte-exact checked
  artifact; scalar C/raw Rust declaration parity remains 253/253.
- `cargo test --release --features full-rust-bindings
  'multi_point::test_' --jobs 1`: 16/16 passed.
- `compressed_infinity`: 1/1 passed; scalar conformance: 2/2 passed.
- Full-binding smoke example passed.
- Full-binding rustdoc passed with warnings denied; Rust safe-API contract gate
  passed.
- `bash -n bindings/rust/publish.sh` passed. Dirty source was rejected before
  link creation; a clean native POSIX fixture with injected `/usr/bin/true`
  Cargo left no link or isolated target behind. A clean MSYS fixture without
  native-link permission refused and left no emulated directory copy.

No benchmark, native rebuild, full suite, sanitizer, CUDA, BTx/BTE audit,
portable-AVX Release, native AVX, Apple CI, push, or live Cargo publication ran.

## Residual risk and decision

The five files are complete at their recorded implementation commits. A real
registry dry run/package installation, Apple execution, sanitizer/fuzzing, and
the broad release matrix remain milestone evidence. Benchmark support remains
explicitly outside this review.
