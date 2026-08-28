# File Review: `blst/avx/curve/point.h`

## Record

| Field | Value |
| --- | --- |
| Path | `blst/avx/curve/point.h` |
| Roadmap phase / queue | Phase 16 AVX curves and protocols |
| Reviewed commit | `713a012583657102c29a4287783585f822baecc7` |
| Reviewer | Codex audit team |
| Date | 2026-08-28 |
| Disposition | in-scope production |
| Subsystem | AVX G1/G2 point representation, arithmetic, encoding, validation, and signing |
| Profiles | portable AVX Debug, scalar-default C/C++, AVX-default C/C++, native AVX-512/IFMA Release |
| Generation/provenance | VROOM AVX backend adapted to CatBLST homogeneous points and selected API |
| Upstream references | VROOM import lineage; no direct upstream BLST AVX counterpart |
| Related files | `blst/avx/curve/ec.h`, `blst/avx/curve/scalar_mult.h`, `blst/c_api.cpp`, `blst/scalar/curve/point.h`, `bindings/blst_scalar.h`, `test/avx/Avx_c_api_contract_tests.cpp` |

## Purpose and trust boundary

This header owns the AVX homogeneous G1/G2 projective and affine types and the
selected native operations used by native C++ and AVX-default C builds. Point
coordinates are public for ordinary validation, encoding, and aggregation;
signing and public-key derivation additionally consume a secret scalar. The
file must preserve subgroup, scalar, timing, alias, and failure-publication
contracts without converting whole points to the scalar backend.

## Contracts and invariants

- Basic point primitives accept on-curve points; optimized multiplication,
  signing, and advanced algorithms require prime-order subgroup points.
- `sign`, `sign2`, `to_pk`, and `to_pk2` require a canonical nonzero secret key
  and perform no Release validation. Their scalar-dependent schedules are
  fixed for valid keys.
- The two inherited `unchecked_mult` methods apply exact scalars to on-curve
  points without subgroup endomorphisms.
- Projective storage is homogeneous: affine coordinates are `X/Z` and `Y/Z`;
  any zero Z is infinity. Affine infinity is the all-zero representation.
- Exact whole-object aliasing is supported where documented. Conditional
  operations and `times_minus_z` require distinct objects.
- Raw and standard decoding publish output only on success. Decoding establishes
  on-curve canonical storage but not subgroup membership unless stated.
- `sign2` and `to_pk2` optionally return full-width pre-decorated `X || Y`
  buffers, not directly conforming compressed encodings.

The public C declarations are documented in `bindings/blst_scalar.h`. The
selected native C++ declarations are in this file. Safe Rust hashes to the
signature subgroup and owns a validated secret-key type; generated raw Rust
inherits the C contract.

## Dependency and call-path map

`blst/curve/point.h` selects these types when `BLST_DEFAULT_AVX` is enabled.
The C boundary converts exact public storage to and from these types. Curve
formulas come from `ec.h`; GLV/GLS and fixed-window scalar multiplication come
from `scalar_mult.h`. Scalar point code is an independent differential oracle,
not a production whole-point fallback. Pairing, mapping, aggregation, and MSM
consume these point types through separately reviewed entry points.

## Upstream/specification comparison

The AVX implementation descends from the VROOM import rather than upstream
BLST. Public C spellings and the two exact-scalar exceptions retain BLST
compatibility. CatBLST adds selected native C++ parity, explicit domain and
timing contracts, direct AVX observation, and failure-atomic decoding.

## Checklist evidence

### 1. Identity and purpose

The complete declaration/definition file and its supported selected callers
were inspected. Generic helper visibility is intentional under `AUD-167`.

### 2. Contracts and invariants

Point domains, homogeneous representation, scalar/key ranges, timing intent,
aliasing, output sizes, optional outputs, and decode publication are explicit.

### 3. Upstream and specification lineage

VROOM lineage is recorded in the file header and manifest. There is no direct
upstream BLST AVX implementation to merge mechanically.

### 4. Correctness

Exact-scalar multiplication, decode publication, equality, signing, key
derivation, optional output, and helper aliasing were compared with scalar
oracles. Projective equality normalizes RNS cross-product differences before
zero testing and separately requires matching infinity status.

### 5. Memory and lifetime safety

Types own fixed-size coordinate storage. No allocation occurs. Conversion
boundaries use typed values and explicit copies. Decoders publish local results
only on success.

### 6. Undefined and implementation-defined behavior

Control-bit and distinct-object constraints use Debug assertions. Supported
builds require 64-bit limbs and little-endian storage as recorded by the
platform audit.

### 7. Side channels and cryptographic hygiene

Secret signing uses fixed-schedule GLV/GLS multiplication and a fixed
addition-chain final inversion. Temporary scalar bytes and scalar-field
inversion values are wiped. Direct equality and zero predicates avoid
inversions and scalar point conversion. Native optimized assembly inspection
is still required before the signing fix is validated.

### 8. Randomness and secret lifecycle

N/A for entropy: this file consumes caller-owned keys. Secret-derived temporary
arrays and inversion values are cleared; key generation is owned elsewhere.

### 9. API and ABI

The selected AVX types now provide parity for `sign`, `sign2`, `to_pk`,
`to_pk2`, and both inherited `unchecked_mult` exceptions. No public C symbol or
storage layout changed.

### 10. Concurrency

Point operations use caller-owned state and immutable constant-initialized
curve owners. No mutable global pool, guarded first-use initialization, or
runtime dispatch is introduced.

### 11. Resource exhaustion and denial of service

All reviewed operations have fixed-size stack storage and bounded scalar work.
No heap allocation or input-sized recursion occurs.

### 12. Portability and build behavior

Scalar-default and AVX-default portable-AVX Debug builds compile and run. Native
AVX-512/IFMA runtime and code generation remain evidence gates.

### 13. Performance-sensitive structure

Whole-point scalar signing fallbacks and equality inversions are removed.
Multiplication and coordinate arithmetic remain AVX-native; only one final Z
inversion crosses the field boundary. `sign2` shares that normalized result and
`times_minus_z` avoids an unconditional protective point copy.

### 14. Tests and evidence

Both focused portable-AVX Debug configurations pass 12/12 contract tests.
Coverage includes low-order exact multiplication, decode failure atomicity,
point domains, equality/infinity representations, signing/key derivation,
identity, exact aliasing, optional outputs, C routing, encoding finalization,
and `times_minus_z`.

### Specialized additions

Secret-key assertions and scalar differential tests cover G1 and G2. Native
code-generation inspection is retained because portable emulation cannot prove
AVX-512 register, spill, or instruction shape.

## Findings

`AUD-088`, `AUD-089`, `AUD-093`, `AUD-094`, `AUD-160`, `AUD-164` through
`AUD-173`, `AUD-225`, `AUD-228`, and `AUD-229`.

## Validation

```text
cmake --build out/build/audit-phase1-avxportable --target blst_avx_contract_tests --parallel 1
out/build/audit-phase1-avxportable/test/blst_avx_contract_tests.exe
cmake --build out/build/audit-phase3a-avxportable --target blst_avx_contract_tests --parallel 1
out/build/audit-phase3a-avxportable/test/blst_avx_contract_tests.exe
cargo fmt -- --check
scripts/audit/check_public_api_documentation.ps1
scripts/audit/check_group_api_contracts.ps1
scripts/audit/check_audit_manifest.ps1
git diff --check
```

Result: both focused test configurations passed 12/12; documentation,
formatting, audit-manifest, and whitespace checks passed. Structural inspection finds
no whole-point scalar conversion in AVX signing or key derivation.

The `AUD-229` follow-up replaces guarded construction of the generator,
endomorphism, A-prime, and Frobenius owners with constant-initialized MatrixNoK
values. The separately checked curve-literal header regenerates from canonical
BLS12-381 coordinates and the existing RNS parameters. The focused portable
target passes 35/35, including direct scalar comparisons for both generators,
both negative generators, both A-prime values, and Fp12 identity. MSVC
production objects contain no dynamic-initialization symbols. Peach Clang
21.1.8 native Release passes the consolidated 79-test AVX suite and all five
MatrixNoK audit gates. The optimized probe contains no selected-vector wrapper
or conversion symbol and retains the accepted hot-owner stack bounds.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Prior supported baseline | No packet-specific diagnostic | Template mathematics, constant-time behavior, and code-generation shape require direct review |
| Clang Static Analyzer / clang-tidy | Prior baseline | No packet-specific diagnostic | Does not prove group formulas or constant time |
| Clippy | Safe Rust declaration review; `cargo fmt -- --check` | Documentation updated | No Rust implementation change |
| ASan/UBSan | Not run in this packet | pending milestone evidence | User-gated broad run |
| Fuzzing | Not run in this packet | pending milestone evidence | Planned campaign; does not prove secret timing |

## Residual risk and decision

Complete for this file at the reviewed source commit. `AUD-169` is validated by
native AVX-512/IFMA optimized code generation and runtime evidence. Sanitizer,
fuzz, statistical-leakage, and broad release matrices remain subsystem gates;
they do not reopen the completed line-by-line review. Apple is a supported
scalar platform and does not compile this x86 AVX owner. No benchmark was run.

## 2026-08-28 native validation, codec, and batch-normalization addendum

Commits `4e5a47a5` and `40ac7710` replace whole-point scalar fallbacks with
native equality, curve equations, subgroup relations, and canonical codecs.
For homogeneous storage, the native projective checks evaluate
`Y^2 Z = X^3 + b Z^3`; every zero-Z representation remains infinity. The G1
membership check evaluates `sigma^2(P) = -[z^2]P`, and the G2 check evaluates
`psi(P) = -[z]P`. The accepted `CurveMembership` chain proves that each
relation implies prime-order membership under the concrete BLS12-381
hypotheses recorded by `AUD-231`. Affine entry points first check the curve
equation and then apply the same projective relation.

Canonical serialization is deliberately a representation boundary. Each
native field coordinate is converted to canonical scalar-field limbs only for
sign selection or big-endian output. Decoding rejects noncanonical field
elements before conversion to native storage, reconstructs compressed Y, and
publishes output only after all checks succeed. G1 and G2 encodings retain the
scalar CatBLST byte order and infinity flags, so scalar-only nodes remain wire
compatible.

Commit `713a0125` adds caller-buffered batch normalization for public G1 and G2
points. The forward pass stores prefix products and skips zero-Z inputs; the
reverse pass uses one variable-time inversion and emits the all-zero affine
infinity for each skipped input. Output and input lengths must match, and the
workspace must contain at least one field element per input. No allocation or
global scratch is introduced.

On Peach, Clang 21.1.8 passes 43 focused native tests covering point contracts,
both batch-normalization paths, single-point multiplication, MSM, and the C
boundary. The Rust typed G1/G2 round-trip matches scalar wire encodings. Direct
compilation of `CurveMembershipCardinality` imports and checks the complete
operation, group bridge, generator, dispatch, and membership proof chain.
`AVX-P4-002` and `AVX-P4-003` are rebound to the exact current `point.h` owner.
