# Phase 6A Legacy C++ Wrapper Removal Evidence

## Record

- Date: 2026-07-31
- Implementation commit: `c8f69195957eb8d3da90482c2c1a26eeabec6b7c`
- State: validated and committed
- Scope: `AUD-006`, upstream rows `e61acc2` and `d4ca0f2`
- Benchmarks: excluded and not run
- CUDA: disabled and not exercised

## Decision

The deleted `bindings/blst.hpp` was a class-oriented wrapper over the C ABI,
not the supported native C++ implementation. It duplicated the C surface while
adding ownership, exception, allocation, and lifetime contracts that were not
safe enough to retain. The project has no compatibility requirement for this
pre-production surface, so removing it is lower risk than repairing and
maintaining a second C++ API.

The upstream `e61acc2` and `d4ca0f2` changes apply only to that wrapper. They
remain `not_applicable`: importing polish or correcting its language-standard
claim would not address the wrapper's underlying design and safety issues.

Native C++ support remains in the canonical field, curve, and protocol types
under `blst/`. The generic-modulus `blst/blst_t.hpp` shim was initially retained
here, then retired by `AUD-078` at `a720184` after the field lifetime audit
confirmed that it duplicated low-level arithmetic, retained separate lifetime
violations, and had no production consumer.

## Removal review

The review found that the wrapper contained:

- constructors that built and discarded temporaries instead of initializing
  the destination object;
- empty-vector indexing through element zero;
- unchecked size shifts and hidden heap allocation in convenience paths;
- manual pairing-context and DST storage with unsafe implicit copy/move
  behavior;
- copyable secret-key storage without an explicit erasure contract; and
- raw enum exceptions and legacy SWIG/C++ compatibility machinery.

The implementation commit therefore:

- deletes `bindings/blst.hpp` and its entirely commented-out test;
- removes the dead test from CMake and the wrapper from Rust package contents;
- removes the wrapper-specific C++ branch from `bindings/blst.swg`, leaving
  its pre-existing C-API route until SWIG and unsupported bindings are removed;
- removes wrapper-specific manifest classification rules; and
- corrects the repository README to describe native C++ as the supported
  implementation and API.

A repository-wide non-documentation search found no remaining `blst.hpp` or
`BindingCpp` reference. Historical audit references remain so the decision and
upstream provenance are not erased.

## Validation

| Surface | Configuration | Result |
| --- | --- | --- |
| C and native C++ | Windows MSVC Debug, direct ADX scalar | 282 registered: 273 passed, nine expected AVX-only skips, zero failed |
| Native C++ packaging/build | Windows MSVC Debug, portable AVX-512 default | `blst_tests` configured and built successfully |
| Rust | `cargo test --features full-rust-bindings` | 71 unit tests and one compressed-infinity integration test passed |
| Rust package contents | `cargo package --list --allow-dirty` | No deleted wrapper path is packaged |
| Formatting/references | Cargo formatting, Git whitespace, repository search | Passed; no active wrapper reference remains |

## Residual validation

- SWIG is not installed on this host, so the retained C-API SWIG path was not
  parsed or regenerated. SWIG, Python, Java, and other unsupported binding
  surfaces remain planned removals rather than supported release paths.
- No ARM64 or WebAssembly build was run for this removal-only batch.
- This phase removes a file rather than marking a production implementation
  file complete. Native C++ implementation files remain subject to the normal
  file-by-file audit.
