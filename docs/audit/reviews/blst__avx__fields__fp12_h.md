# `blst/avx/fields/fp12.h` Review

## Scope

This record covers the complete retained selected AVX Fp12 arithmetic owner
through commit `8201c087`. It excludes BTx/BTE and CUDA.

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
- Dense multiplication writes coordinates 0--5 and 6--11 through two fixed
  non-inlined helpers. Both helpers write into one fresh result and do not
  inspect field values.
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
- `AUD-401`: monolithic dense multiplication retained excessive live AVX-512
  state and generated material spill traffic under Clang and GCC.

## Validation

- The structural validator rejects every removed schedule, alias, constructor,
  dead helper, loose bound, and `concat` restoration.
- A serial two-worker portable-AVX Debug rebuild succeeds. Five focused Fp12,
  GT, final-exponentiation, Miller, and pairing differentials pass.
- `Fp12Inverse.lean` proves the source-shaped inverse pipeline from exactly the
  supplied base-inverter contract. Direct dense inverse/product/round-trip and
  individual-versus-batch final-exponentiation tests pass.
- `Fp12Raw.lean` models the two dense-multiplication writers and their ordered
  join. The source-shaped multiplication, sparse-Miller, and cyclotomic modules
  type-check after the split.
- Correspondence rows `AVX-P3-002` through `AVX-P3-005` and `AVX-P3-008` are
  rebound to commit `8201c087` and the exact `fp12.h` hash. `AUD-402` tracks
  eight unrelated stale mappings that still block the repository-wide gate.
- MSVC `/O2 /arch:AVX512` compile-only inspection shows the selected square
  wrapper has a 40-byte frame, two split-helper calls, and no coordinate-copy
  instructions or `memcpy`. Before repair it had a 3,168-byte frame and copied
  a twelve-coordinate temporary through `memcpy`.
- Serial Peach measurements show an 8.7%--9.3% Clang reduction and a
  10.2%--10.7% GCC reduction for dense multiplication across all four retained
  input-bound combinations. Representative Clang consumers improve by
  2.3%--3.6%.

## Residual risk

The selected base inverter remains variable-time and zero remains outside its
contract. The full narrow-to-standard probe still has a final 1,536-byte copy
at the generic returned-value assignment boundary; it is not in the selected
cyclotomic owner and should be considered with that generic boundary's own
file review. Native compiler/IFMA semantics remain P0 evidence. Sanitizers,
CUDA, BTx/BTE audit, and portable-AVX Release remain outside this packet.
