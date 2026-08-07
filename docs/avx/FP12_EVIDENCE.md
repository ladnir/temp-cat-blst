# MatrixNoK Fp12 and GT Candidate Evidence

## Packet boundary

- Candidate implementation: `1011916baf69a4ecb4b0e012c30a74051c8ec850`
- Sparse/final-exponent follow-up: `3e1f99a`
- Upstream algorithm commit: `d6ecd534c54a41702ea80ae75b52250b6d51b372`
- Compared VROOM head: `26bd9aee2168f167583cf096add355547a7b56d4`
- Generic field owner: `blst/avx/fields/fp12.h`
- MatrixNoK policy: `blst/avx/fields/matrixnok_fp12.h`
- GT exponent owner: `blst/avx/gt/glv4.h`

This packet adds an isolated production-shaped Fp12 and GT candidate. It does
not change a public symbol, type, or production dispatch. The existing Fp12
template remains source-compatible through a default constants policy, while
the candidate selects a stateless MatrixNoK policy whose fixed values are
checked-in, aligned, and independently regenerated.

## Correctness design

The focused native suite compares decoded candidate results with scalar
CatBLST. It covers general Fp12 multiplication and squaring, multiplication by
one, Frobenius powers one through three, conjugation, cyclotomic squaring on a
valid GT element, and constant- and variable-time 256-bit GT exponentiation.
The follow-up also covers the Miller sparse operand layout, exact in-place
aliasing, and the complete hard final-exponentiation chain on three independent
inputs. The inverse is computed by scalar CatBLST and imported separately, so
the arithmetic chain is tested independently of a future MatrixNoK inversion
boundary. Together with the previously retained Fp/Fp2/G1 cases, all 14 focused
`AvxMatrixNoKPrimitives.*` tests pass under Clang 21.1.8 on native
AVX-512IFMA.

The constant validator independently reconstructs the MatrixNoK Fp12 offset,
one-third, negative-two-thirds, and Frobenius representations from the fixed
moduli. It verifies every generated value and its 64-byte alignment.

The first generic Fp12 translation exposed `AUD-216`: the old template chose
wide multiplication offsets from `ELEMENT_BITS` alone. The incumbent and
MatrixNoK rings both use 50-bit limbs but use different RNS moduli, so that
dispatch compiled successfully and silently performed arithmetic with the
wrong basis. Constants are now selected by an explicit ring-parameter policy;
the two 50-bit rings cannot accidentally share fixed values.

## Native code generation

Clang 21.1.8 compiled the probes with `-O3 -DNDEBUG -march=native`. The selected
Fp12 multiplication body is approximately 14.4 KiB (`0x381e`) versus VROOM's
`0x383d`. Both contain exactly 492 low IFMA operations, 456 high IFMA
operations, and 12 `mulx` instructions. The candidate contains no dynamic
initialization guard, vector-conversion wrapper, or parameter-context load.

Code-generation review also exposed `AUD-217`: forced unrolling of the outer
eight-row R2 change-base loop expanded the candidate body to approximately
21.9 KiB and increased register pressure. Retaining the fixed inner arithmetic
but allowing one fixed-count outer loop restored VROOM's arithmetic counts and
code size. The scalar batch buffer is explicitly 64-byte aligned.

## Native performance gate

One Google Benchmark process ran at a time on Peach's AMD Ryzen 9 7950X. The
same Clang 21.1.8 Release binary contained the incumbent, this candidate, and
actual VROOM commit `26bd9aee`. Google Benchmark reported CPU scaling enabled,
so these are selection results rather than final release claims.

| Operation | Incumbent | CatBLST MatrixNoK | VROOM | vs incumbent | vs VROOM |
| --- | ---: | ---: | ---: | ---: | ---: |
| Fp12 multiply | 362 ns | 303 ns | 304 ns | -16.3% | -0.3% |
| Cyclotomic square | 237 ns | 207 ns | 198 ns | -12.7% | +4.5% |

Representative inverse GT FFT results were:

| Elements | Incumbent | CatBLST MatrixNoK | VROOM | vs incumbent | vs VROOM |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 256 | 35.2 ms | 30.3 ms | 30.0 ms | -13.9% | +1.0% |
| 1,024 | 191 ms | 165 ms | 160 ms | -13.6% | +3.1% |

The candidate materially beats the incumbent and matches VROOM for the core
multiply. The cyclotomic primitive and 1,024-element FFT remain outside the
preferred 2% VROOM tolerance and are recorded as residual optimization work;
they do not invalidate the correctness or the decision to continue the port.

## Rejected alternatives

- Splitting correction operations into a separate batch phase worsened Fp12
  multiplication from approximately 303 ns to 323 ns.
- Caching cyclotomic constants in local references worsened the specialized
  square without helping multiplication.
- Storing constants behind static reference-returning accessors reduced one
  apparent frame but regressed the measured cyclotomic path.

All rejected arithmetic variants were removed. The retained multiply reads
the wide offset once per operation; the cyclotomic path keeps its direct
compile-time constant uses.

## Exact correspondence inputs

- `fp12.h` SHA-256:
  `e54d1ef61eb8ed7ee36bbcac4d4acd00327de8368462c847414248a6d17194d7`
- `matrixnok_fp12.h` SHA-256:
  `4cf4bd71b4387f7e52f6f1cf105d62d4943f7469800702ca0f85720528e4bcd9`
- `matrixnok_fixed.h` SHA-256:
  `73ba31d8a347208ac84173790a61929a34f443a0d47783cffa21559661ac026f`
- `matrixnok_fp.h` SHA-256:
  `11c172dccd76dba2d923ea205201b32a1e6bdf679f63f3d1560016d4dcc017d7`
- `generated_matrixnok_constants.hpp` SHA-256:
  `d1bd0045e8d7f223cb332842bd4c9023daefc715079c0649fbce953cabdfea43`
- `glv4.h` SHA-256:
  `8765c031be885fed24b8743c7ee677ca02640c3f60423bd9d8ae53fa87cc1d31`

## Residual gates

- Production dispatch remains unchanged until the retained G2, MatrixNoK
  inversion boundary, Miller-loop, and pairing slices pass their gates. Sparse
  Fp12 and the hard final-exponentiation arithmetic are now validated.
- The cyclotomic and 1,024-element FFT deltas remain explicit optimization
  targets; no benchmark claim is made under stable release conditions yet.
- Fp12/RNS arithmetic and every bound transition remain formal proof
  obligations under `PROOF_PLAN.md`.
- Sanitizers, fuzzing, broad suites, Apple CI, CUDA, BTx/BTE, and portable-AVX
  Release were outside this focused packet.
