# `blst/avx/fields/fp12.h` Review

## Scope

This record covers the complete retained selected AVX Fp12 arithmetic owner
through commit `793b542`. It excludes BTx/BTE and CUDA.

## Contracts

- `CommonElement<N>` is twelve flattened base-field coordinates representing
  six Fp2 coefficients in the scalar tower's established order. `N` is a
  compile-time producer bound, not runtime validation.
- General multiplication, square, Frobenius, conjugation, sparse Miller
  multiplication, and nonzero inversion accept valid Fp12 values.
- Cyclotomic square has the narrower cyclotomic/GT precondition. It performs no
  membership check in Release.
- Inversion requires nonzero input and inherits the supplied base inverter's
  timing. The selected provider is deliberately variable-time.
- Operations are allocation-free and fixed-shape with respect to field data.
  Returned-value operations support exact caller-side storage reuse. The
  sparse output-reference operation reads all inputs before assigning output.
- The two cyclotomic halves write into one final owner. Its input stays distinct
  until both halves complete, so exact caller-side input/output reuse cannot
  expose a partially overwritten value.

## Findings

- `AUD-268`: unselected versioned general and cyclotomic schedules enlarged the
  proof and maintenance surface.
- `AUD-269`: stateless constants were modeled as inherited runtime policy and
  selected aliases were duplicated.
- `AUD-270`: the concrete inversion decomposition lacked direct proof and
  differential evidence.
- `AUD-271`: returned half-arrays and generic concatenation survived optimized
  MSVC lowering as large stack temporaries and a copy.
- `AUD-272`: split reduction used loose bound 200 and dead tuple/index plumbing.
- `AUD-273`: provenance, representation, timing, alias, and domain contracts
  were incomplete.

## Validation

- The structural validator rejects every removed schedule, alias, constructor,
  dead helper, loose bound, and `concat` restoration.
- A serial two-worker portable-AVX Debug rebuild succeeds. Five focused Fp12,
  GT, final-exponentiation, Miller, and pairing differentials pass.
- `Fp12Inverse.lean` proves the source-shaped inverse pipeline from exactly the
  supplied base-inverter contract. Direct dense inverse/product/round-trip and
  individual-versus-batch final-exponentiation tests pass.
- Accepted correspondence rows `AVX-P3-002` through `AVX-P3-005` are rebound,
  and new `AVX-P3-008` owns inversion. All 21 rows and 61 mapped hashes validate.
- MSVC `/O2 /arch:AVX512` compile-only inspection shows the selected square
  wrapper has a 40-byte frame, two split-helper calls, and no coordinate-copy
  instructions or `memcpy`. Before repair it had a 3,168-byte frame and copied
  a twelve-coordinate temporary through `memcpy`.

## Residual risk

The selected base inverter remains variable-time and zero remains outside its
contract. The full narrow-to-standard probe still has a final 1,536-byte copy
at the generic returned-value assignment boundary; it is not in the selected
cyclotomic owner and should be considered with that generic boundary's own
file review. Native compiler/IFMA semantics remain P0 evidence. No benchmark,
remote execution, broad suite, sanitizer, CUDA, BTx/BTE audit, or portable-AVX
Release build ran for this packet.
