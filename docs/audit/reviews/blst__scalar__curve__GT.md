# File Review: `blst/scalar/curve/GT.h` and `GT.cpp`

## Record

| Field | Value |
| --- | --- |
| Paths | `blst/scalar/curve/GT.h`; `blst/scalar/curve/GT.cpp` |
| Roadmap phase / queue | Phase 15 scalar curves and protocols, sixth packet |
| Reviewed commit | `09722691117bfbae36baa27da0c6f6481fcd58ad` |
| Reviewer | Codex audit team |
| Date | 2026-08-05 |
| Disposition | in-scope production |
| Subsystem | scalar target group, final exponentiation, and Miller loops |
| Profiles | scalar host; selected scalar/AVX facade; CUDA excluded |
| Generation/provenance | C++ translation of upstream BLST pairing/target-group code with CatBLST native APIs and batching |
| Related files | AVX GT, C boundary, pairing aggregation, FFT/middle product, C/Rust contracts |

## Purpose and trust boundary

These files implement scalar GT arithmetic, target-group scalar
multiplication, final exponentiation, direct and prepared Miller loops, and the
canonical pairing generator. Pairing inputs are trusted on-curve affine points;
these hot operations do not validate curves or subgroups. A Miller result is a
general Fp12 pairing accumulator and is not a GT element until final
exponentiation.

## Contracts and invariants

- `GT` means the prime-order target subgroup. Its optimized arithmetic assumes
  that membership was established by construction or validation.
- `MillerResult` is the distinct pre-final domain. It exposes only identity,
  multiplication, and conjugation before `GT::final_exp`.
- `final_exp` requires a valid nonzero Miller/general-Fp12 value and produces
  GT. Zero remains outside the legacy void API's supported domain.
- Direct Miller inputs must be valid on-curve affine points. Infinity contributes
  the multiplicative identity. No curve or subgroup validation is performed.
- The pointer-array overload retains BLST's compact continuation convention;
  malformed nonempty heads fail closed at the public boundary.
- Prepared scalar line evaluation stores only `(-2x, 2y)` in
  `PreparedG1LineScale`; it does not manufacture an invalid affine point.
- GT conditional negation is constant-time for a public zero-or-one condition
  bit and selects only the Fp12 half changed by conjugation.
- Checked GT decoding rejects noncanonical coefficients and non-GT values
  without modifying output. `deserialize_unchecked` still requires canonical
  coefficients and defers only membership; its candidate is not an ordinary GT
  operand until a later aggregate membership check succeeds.
- The scalar and AVX generator is exactly final-exp(pairing(G1 generator,
  G2 generator)) and uses static compile-time storage.

## Dependency and call-path map

Pairing aggregation, the C API, FFT/middle product, BTE compatibility
consumers, tests, and optional benchmark sources call the selected scalar/AVX
surface. The C and generated Rust ABI continue to transport Miller and GT
values in `blst_fp12` storage, with the semantic role documented at each
operation. Ordinary C boundaries use direct field conversion. The middle-
product selected-backend bridge retains a serialize/decode conversion and uses
the canonical unchecked form because the source was already established as GT.

## Checklist evidence

1. Every declaration, definition, helper, generated constant, selected-backend
   caller, and retained host test was traced.
2. GT, Miller, prepared-line, pointer-array, zero/infinity, alias, and timing
   contracts are explicit.
3. Pairing schedules and field formulas were compared with the retained
   upstream lineage; CatBLST batching and selected-backend adaptations were
   reviewed separately.
4. Direct, multi, prepared-line, identity, bilinearity, final verification, and
   generator equivalence paths pass focused tests.
5. Fixed scalar batches remain bounded at 16 entries; no allocation, runtime
   dispatch, thread pool, or new dynamic initialization was introduced.
6. Faux affine storage and whole-wrapper byte copies were removed from the
   prepared-line and identity initialization paths.
7. Constant-time GT conditional negation remains branchless. Miller scheduling
   depends only on public pair counts and public infinity structure.
8. This subsystem consumes no randomness and owns no long-lived secret state.
9. Native C++ distinguishes Miller from GT; C/Rust compatibility storage is
   documented and converted explicitly at boundaries.
10. All state is automatic, caller-owned, or immutable static storage.
11. Pair batches are fixed-size; malformed C pointer heads fail closed.
12. Scalar, portable-AVX Debug, and native AVX-512 production compilation pass.
13. Hot pairing schedules, line formulas, and batch widths are unchanged;
   boundary serialization copies and unnecessary cneg work were removed.
14. Focused evidence is recorded below; broad milestone gates remain deferred.
15. Canonical decoding checks all twelve coefficients before publication;
    checked decoding performs one GT membership test, while trusted backend
    conversion structurally calls only `deserialize_unchecked`.

## Findings

`AUD-081`, `AUD-085`, `AUD-099`, `AUD-103`, and `AUD-131` through `AUD-134`.

## Validation

- MSVC scalar Debug and Release pass the dedicated GT decoder target 5/5.
- Selected portable-AVX Debug passes the same decoder target 5/5.
- Focused scalar pairing/GT/C API/precomputed-line/middle-product packet passes
  24/24.
- MSVC portable-AVX Debug production library and five directly affected test
  translation units compile serially.
- Native AVX-512/IFMA production library compiles serially.
- Generator equivalence and absence of a native dynamic initializer were
  established in commit `f7100402573f429a1f2132d85dd03ba67829cce7`.
- No benchmark, full suite, sanitizer sweep, CUDA build, Apple CI, or
  portable-AVX Release build ran.

## Residual risk and decision

The scalar files are completely reviewed at the recorded commit. C/Rust must
retain representation-only `blst_fp12` storage for compatibility, and explicit
native raw-field construction remains an intentional expert escape hatch.
`AUD-085` is validated by canonical-boundary, membership, failure-atomic, and
deferred-aggregate tests. Experimental CUDA was deliberately not migrated to
the new native C++ decoder/Miller contracts and remains outside the scalar/AVX
production-readiness gate.
