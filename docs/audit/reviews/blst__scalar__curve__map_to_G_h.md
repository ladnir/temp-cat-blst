# File Review: map-to-group headers

## Record

| Field | Value |
| --- | --- |
| Paths | `blst/scalar/curve/map_to_G.h`; `blst/avx/curve/map_to_G.h`; `blst/curve/map_to_G.h` |
| Roadmap phase / queue | Phase 15 scalar curves and protocols, seventh packet |
| Reviewed commit | `427980fb80687582022cc72bb402e1a8d40d779f` |
| Reviewer | Codex audit team |
| Date | 2026-08-05 |
| Disposition | in-scope production |
| Subsystem | RFC 9380 hash/encode/map to BLS12-381 G1 and G2 |
| Profiles | scalar; selected portable/native AVX facade; C/C++/generated Rust contracts |
| Generation/provenance | C++ translation of upstream BLST hash-to-curve formulas with CatBLST selected-backend wrappers |
| Upstream/specification | pinned BLST lineage; RFC 9380 BLS12-381 XMD:SHA-256 SSWU suites |
| Related files | `blst/scalar/fields/hash_to_field.h`, `blst/c_api.cpp`, `bindings/blst_scalar.h`, `bindings/rust/src/bindings.rs`, `test/Map_to_group_audit_tests.cpp` |

## Purpose and trust boundary

The scalar header implements simplified SWU, the G1/G2 isogenies, cofactor
clearing, and message-to-subgroup composition. The selected AVX header is an
explicit compatibility bridge because no native AVX SSWU implementation exists.
Messages, augmentations, and DSTs are public byte strings. Direct map inputs are
trusted reduced Montgomery Fp/Fp2 values; these void APIs do not validate field
representations in Release.

## Contracts and invariants

- `Encode_to_curve` implements the RFC 9380 nonuniform suite and maps one field
  element; `Hash_to_curve` implements the random-oracle suite and maps two.
- Successful operations produce on-curve prime-order G1/G2 points after
  cofactor clearing. No Release postcondition check is performed.
- Byte pointers cover their stated lengths and may be null only for zero length.
  Operations are allocation-free and nonthrowing for valid inputs.
- Direct `map_to_curve` inputs are reduced Montgomery values. A null second
  value selects the one-element form.
- `clear_cofactor(out, p)` requires distinct objects because it consumes `p`
  after beginning to write `out`; Debug asserts this invariant.
- `isogeny_map_to(out, p)` supports exact aliasing. The supported template
  groups are compile-time restricted to G1 and G2.
- Hash-to-field writes name actual Fp/Fp2 subobjects and statically verify that
  the destination array contains exactly `count * extension_degree` outputs.
- Generated raw Rust declarations inherit the authoritative C contract. There
  is no separate safe Rust map-to-group wrapper in the reviewed surface.

## Dependency and call-path map

The public C encode/hash functions call the scalar implementation directly
because public ABI storage is scalar. Native selected C++ callers use the scalar
or AVX facade according to `BLST_DEFAULT_AVX`. The AVX facade converts
homogeneous and Jacobian projective coordinates with multiplication only; it
does not normalize through affine form. Signature/key aggregation consumes the
resulting subgroup points. Experimental CUDA retains automatic local constants
in device compilation and is outside this release gate.

## Upstream/specification comparison

The fixed SSWU, isogeny, and cofactor schedules were preserved. Exact RFC 9380
empty-message outputs for all four BLS12-381 G1/G2 NU/RO suites are now release
tests. CatBLST intentionally restores immutable host static coefficient storage,
uses typed C++ subobject access, hides raw mapping helpers from the selected
namespace, and keeps the non-native AVX bridge explicit.

## Checklist evidence

1. All declarations, specializations, constants, selected exports, C callers,
   and retained tests in the three headers were traced.
2. Field representation, optional second input, subgroup output, aliasing,
   pointer, allocation, exception, and validation contracts are explicit.
3. The schedules were checked against retained BLST lineage and RFC 9380 vectors.
4. Four exact KATs, zero-field exceptional SSWU, exact isogeny aliasing, and
   scalar/AVX projective round trips pass.
5. Flattening `Fp2[]` through `Fp*` was removed; typed array bounds are checked
   at compile time and no heap storage is used.
6. Every hash-to-field write names its live field subobject; no cross-array
   pointer arithmetic remains.
7. Map scheduling and conditional selection retain fixed structure. Inputs are
   public protocol data; this packet introduced no secret-indexed access.
8. No randomness or long-lived secret state is owned here.
9. C, selected native C++, and generated Rust contracts are recorded in
   `PUBLIC_API_DOCUMENTATION.csv`.
10. Host constants are immutable static data with no dynamic initialization.
11. Work is fixed by the suite and uses no input-sized allocation.
12. Scalar Debug/Release, portable-AVX Debug execution, and native AVX compile
   checks pass. CUDA was deliberately excluded.
13. Host stack frames fell substantially and all AVX bridge inversions and the
   redundant C scalar/AVX round trip were removed.
14. Focused evidence is listed below; broad milestone gates remain deferred.

## Findings

`AUD-135` through `AUD-140`.

## Validation

- MSVC scalar Debug and Release map-to-group targets pass 6/6 each.
- MSVC portable-AVX Debug passes 7/7, including inversion-free G1/G2 projective
  bridge and identity cases; its focused AVX C API contract target passes 5/5.
- Native AVX-512/IFMA Debug compiles serially. The local MSVC configuration is
  compiler-unsupported by policy, so that binary was not executed.
- Optimized MSVC stack frames changed from 4,240 to 1,696 bytes for G1 isogeny,
  2,272 to 1,040 for G2 isogeny, 640 to 448 for G1 SSWU, and 1,024 to 736 for
  G2 SSWU. G1 no longer emits `__chkstk`; coefficient initialization is absent.
- Structural scans find no Fp flattening cast, affine AVX bridge, internal
  selected-facade export, or C encode/hash selected-backend round trip.
- No benchmark, full suite, sanitizer sweep, CUDA build, BTx/BTE work, Apple CI,
  or portable-AVX Release build ran.

## Automated evidence

| Tool/profile | Version and configuration | Result or finding IDs | Blind spots / N/A reason |
| --- | --- | --- | --- |
| CodeQL | Existing local C++ baseline | No additional alert in this surface | Mathematical conformance and projective equivalence require direct tests |
| Clang Static Analyzer / clang-tidy | Existing repository baselines | No packet-specific blocker | MSVC code generation and template domain checks were inspected directly |
| Clippy | N/A | Generated raw FFI inherits C declaration | No safe Rust implementation exists here |
| ASan/UBSan | Deferred milestone evidence | N/A | Focused typed-lifetime repair is covered by compile-time structure and KATs |
| Fuzzing | Deferred subsystem campaign | N/A | Exact RFC vectors do not replace later arbitrary-input fuzzing |

## Residual risk and decision

The three mapping headers are complete for the supported scalar/AVX host gate at
the recorded commit. `hash_to_field.h`, `c_api.cpp`, and the public C header were
changed only in this focused family and remain pending their complete file
reviews. Native AVX hardware execution, sanitizers, fuzzing, and the full
supported matrix remain milestone evidence. Experimental CUDA intentionally
uses automatic device constants and is not production-qualified.
