# File Review: generated scalar Rust FFI

## Record

| Field | Value |
| --- | --- |
| Path | `bindings/rust/src/bindings.rs`; `build/generate-rust-bindings.sh`; `build/bindings_trim.pl` |
| Roadmap phase / queue | Phase 18 generated-Rust audit |
| Reviewed commit | `daa4eb3005fa4c2fef68fd23a84c34dd2518df7d` |
| Reviewer | Codex audit team |
| Date | 2026-08-12 |
| Disposition | generated production plus owning generators |
| Subsystem | scalar C ABI declarations for Rust |
| Profiles | Rust over the scalar ABI; selected scalar/native AVX implementation is opaque to the declarations |
| Generation/provenance | rust-bindgen 0.65.1 with Clang/libclang 18.1.3, Rust 1.65 rustfmt 1.5.1, the repository `rustfmt.toml`, then fail-closed Perl normalization |
| Upstream references | Upstream common base `f48500c1fdbefa7c0bf9800bccd65d28236799c1`; split CatBLST headers require local generation scope |
| Related files | `bindings/blst_scalar.h`, `src/fft.rs`, `src/pairing.rs`, `scripts/audit/check_rust_bindings.ps1` |

## Purpose, contracts, and invariants

The generated file is an unsafe declaration/layout artifact, not an independent
safe API. Semantic contracts come from `blst_scalar.h`; Rust additionally
requires valid pointers, alignment, lifetimes, and initialized outputs. The
generator consumes only that scalar owner header. AVX, FFT, and middle-product
split declarations remain in their handwritten modules.

`BLST_ERROR` and `BLST_UNIQ_RESULT` are transparent integer newtypes so unknown
future C values remain representable. Distinct pairing and uniqueness C tags
remain opaque Rust context types behind their public aliases. Layout tests
cover POD size, alignment, field offsets, result width, and unknown result
values. The generator allowlists public scalar names and excludes C-only enum
width assertions and system headers. The postprocessor requires exact rewrite
counts and aborts on missed or duplicate patterns; output is written through a
sibling temporary directory and moved only after success.

## Checklist evidence

- Identity/lineage/API: exact scalar-header ownership and 253/253 unique symbol
  parity; `blst_pairing_alignof` is generated and the three FFT-domain symbols
  are handwritten once.
- Correctness/layout/portability: 12 ABI-layout tests and the generated
  normality test pass; C99 is the recorded parser mode; result carriers are
  forward-safe `u32` newtypes.
- Memory/UB/side channels: declarations themselves execute nothing. Unsafe
  callers retain the C contract; safe wrappers are reviewed in owning files.
- Concurrency/resources/performance: no runtime behavior. Generator temporary
  storage is bounded by header output and cleaned only inside its validated
  sibling path.
- Tests: two consecutive generations to different output directories compare
  byte-for-byte; the exact inventory gate, Bash syntax, Perl syntax, and
  adversarial transform fixture pass; focused pairing/uniqueness Rust
  consumers pass against a cached native archive.

## Findings

`AUD-189`, `AUD-324`, `AUD-325`.

## Validation

- `check_rust_bindings.ps1`: 253 scalar C declarations and 253 generated Rust
  declarations, no missing/extra symbol.
- `bash -n build/generate-rust-bindings.sh`; `perl -c
  build/bindings_trim.pl`; `test_bindings_trim.sh`: passed.
- Generated ABI/layout tests: 13/13 including the normality check; focused
  pairing/uniqueness wrappers: 13/13.

## Automated evidence and residual risk

The path-scoped CI gate installs the exact bindgen and Rust formatter versions,
selects libclang 18, regenerates to a temporary path, and compares bytes without
compiling the native product. The equivalent local run passed with bindgen
0.65.1 and Clang/libclang 18.1.3. Cross-host libclang minor-version stability
and the first hosted execution remain release evidence, but no known
declaration, layout, or reproducibility gap remains.

At `16e5b49`, the audited complete-add exact-alias wording was added to the C
owner and regenerated through the same pinned pipeline. The artifact remains
byte-exact and declaration parity remains 253/253; no declaration, type, ABI,
or runtime code changed.

At `daa4eb30`, the four explicit pairing-check declarations inherit the
clarified complete-membership wording from the C owner. Declaration parity
remains 253/253; no signature, type, layout, ABI, or runtime code changed.

## 2026-09-03 scalar prepared-G2 range addendum

The generated artifact mirrors the five new scalar prepared-G2 declarations
and preserves `blst_scalar_prepared_g2_batch` as a distinct opaque C tag. The
integrity gate now requires exact 258/258 scalar C/raw Rust declarations and
checks that tag explicitly. Split FFT and middle-product declarations remain
outside this generated file. Semantic and safety contracts remain owned by
`bindings/blst_scalar.h` and the safe `scalar::PreparedG2Batch` wrapper.

## 2026-09-04 alpha-pairing range addendum

The generated artifact mirrors
`blst_scalar_prepared_g2_batch_pairing_range`. The integrity gate now requires
exact 259/259 scalar C/raw-Rust declarations. The raw function inherits its
pointer, range, native-representation, output-publication, and concurrency
contract from `bindings/blst_scalar.h`; the safe wrapper owns slice validation
and initialization state.
