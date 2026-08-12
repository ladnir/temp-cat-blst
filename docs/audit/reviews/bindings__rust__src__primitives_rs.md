# File Review: low-level Rust primitive conveniences

## Record

| Field | Value |
| --- | --- |
| Path | `bindings/rust/src/primitives.rs` |
| Roadmap phase / queue | Phase 18 low-level safe Rust primitives |
| Reviewed commit | `dfd0740` |
| Reviewer | Codex audit team |
| Date | 2026-08-12 |
| Disposition | in-scope production |
| Subsystem | scalar-ABI point multiplication, Fp12/GT, Miller loops, and scalar hashing |
| Profiles | safe Rust over the scalar C ABI; selected scalar or AVX native backend beneath it |
| Generation/provenance | Handwritten Rust conveniences over generated raw declarations with upstream BLST lineage |
| Related files | `bindings/blst_scalar.h`, generated `bindings.rs`, scalar/selected point and GT implementations, hash-to-field, and pairing wrappers |

## Purpose and trust boundary

This file adds allocation-free safe Rust calls over fixed-size raw C carriers.
Callers can construct the public carriers directly, so Rust memory safety does
not establish curve, subgroup, reduced-field, Miller-result, or GT validity.
The safe methods enforce slice shape and bit-count bounds that protect native
memory access, while mathematical preconditions remain explicit and unchecked
in Release.

Messages and DSTs entering scalar hashing are caller-controlled. The Rust
wrapper bounds the complete SHA-256/XMD transcript before FFI, reduces the
48-byte result modulo r, rejects zero through `Option`, and erases the expanded
intermediate.

## Contracts and invariants

- Point and GT scalar inputs are unsigned little-endian bit strings with at
  most 256 meaningful low bits. Checked multiplication is subgroup-optimized
  and modulo r; canonical secret-capable multiplication has the documented
  constant-time guarantee, while variable-time operations require public
  scalars.
- The two inherited `unchecked_mult` methods are the only wider multiplication
  exceptions. They require on-curve points and apply the exact low `nbits`
  scalar without subgroup endomorphisms or validation.
- Point equality requires on-curve inputs, does not require subgroup
  membership, and is not a validator.
- `blst_fp12::default()` is the multiplicative identity, not field zero.
  `Mul` and `MulAssign` are general Fp12 operations over reduced
  Montgomery-form coefficients; they neither require nor establish GT.
- GT scalar multiplication requires a previously established GT input. The raw
  Fp12 storage does not encode that property.
- Miller-loop inputs are previously validated G1/G2 subgroup points. Outputs
  are general Miller-result Fp12 values and become GT only after final
  exponentiation. `finalverify` accepts trusted Miller contributions.
- Every inherent method returning a cryptographic value, encoding, membership
  result, or verification decision is `#[must_use]`. Trait signatures remain
  unchanged.
- XMD transcript lengths remain below 2^64 bits, including fixed framing and
  the oversize-DST prehash. Invalid safe-Rust lengths panic before FFI rather
  than becoming indistinguishable from a legitimate zero scalar.
- All native result carriers are published only after the corresponding C
  operation writes them. No wrapper allocates heap storage, dispatches work,
  or permits an exception to cross the C ABI.

## Unsafe-boundary review

All 23 production `unsafe` blocks were traced to their declarations and native
implementations. Fixed-size result functions write every output before return;
the scalar slice helper establishes `nbits <= 256` and sufficient readable
bytes before multiplication. Static identity pointers are nonnull and
immutable. Fp12 exact aliasing used by `MulAssign` is part of the C contract.

`miller_loop_n` uses one contiguous head followed by a null continuation
sentinel. Scalar and selected AVX C traversal remain on that sentinel while
advancing the resolved point pointer, so the two-entry descriptor arrays are
not read out of bounds. Empty slices return identity without FFI.

No invalid raw mathematical value was found to induce Rust memory unsafety;
such values remain contract violations that can produce unspecified
mathematical output or a Debug assertion.

## Findings

`AUD-042`, `AUD-081`, `AUD-088`, `AUD-089`, `AUD-094`, `AUD-099`, `AUD-133`,
`AUD-160`, `AUD-177`, and `AUD-334` through `AUD-337`.

## Validation

- Scalar Debug native XMD packet: 7/7 focused tests pass.
- Scalar Release safe Rust primitives: 9/9 pass, including scalar/GT CT-VT
  agreement, Miller product agreement, XMD oracle, exact length boundaries,
  and invalid-domain panic.
- Full-binding rustdoc passes with warnings denied; `cargo fmt --check` passes.
- Raw binding regeneration with bindgen 0.65.1 and Rust 1.65 rustfmt is
  byte-identical and retains exact 253/253 scalar declaration parity.
- The group-contract structural gate covers all 243 public C group symbols.
- Optimized MSVC inlines the public XMD length checks and emits no call to the
  helper from the C entry point. No benchmark was warranted because the guard
  is a fixed set of public-length comparisons before hashing begins.

No full suite, benchmark, sanitizer sweep, CUDA, BTx/BTE audit, portable-AVX
Release build, native AVX run, Apple CI, or push was performed.
The scalar Debug target rebuilt serially, and frozen BTx/BTE compiled unchanged
only as an existing target dependency.

## Residual risk and decision

The file is complete at the reviewed commit. Mathematical validation remains a
documented caller responsibility by deliberate compatibility and performance
policy. Higher-level validation-state types and batch subgroup checks remain
separate API-design work. Broader sanitizer, fuzzing, Apple, and supported
matrix evidence remain release gates rather than requirements for this focused
change packet.
