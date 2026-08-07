# MatrixNoK Miller and Pairing Preview Evidence

## Packet boundary

- Candidate implementation: `eaa7822431954980694341dc77cdcd8b12cce7eb`
- Upstream algorithm commit: `d6ecd534c54a41702ea80ae75b52250b6d51b372`
- Compared VROOM head: `26bd9aee2168f167583cf096add355547a7b56d4`
- Candidate boundary: `blst/avx/pairing/matrixnok_pairing.h`
- Generic owners: `blst/avx/pairing/{miller,pairing}.h`
- Conversion/inversion boundary: `blst/avx/fields/conversion_inversion.h`

This packet makes the existing Miller-loop and pairing templates accept the
fixed MatrixNoK Fp12 policy. Existing one-parameter instantiations retain their
default types and production behavior. The new aliases are deliberately
isolated: no public symbol, selected type, or production dispatch changes.

## Correctness design

The focused test constructs the G1 and G2 generator, double, and triple with
scalar CatBLST, then converts only the affine inputs into MatrixNoK storage. It
checks four independent pairing combinations: `(G1,G2)`, `(2G1,G2)`,
`(G1,2G2)`, and `(3G1,3G2)`.

For each combination it performs two comparisons:

- candidate Miller output is converted to scalar Fp12, passed through the
  independent scalar final exponentiation, and compared with scalar CatBLST's
  Miller-loop/final-exponentiation result; and
- the complete candidate pairing, including the RNS-to-scalar inversion
  boundary and MatrixNoK final exponentiation, is compared with the same
  scalar GT result.

Local portable-AVX Debug compiled the focused translation unit. On Peach,
Clang 21.1.8 compiled and ran all 17 `AvxMatrixNoKPrimitives.*` tests under
native AVX-512IFMA. The complete retained primitive, field, G1, G2, Fp12/GT,
final-exponentiation, Miller, and pairing set passed.

## Native code generation

Clang 21.1.8 compiled the focused Release objects with `-O3 -DNDEBUG
-march=native`. The candidate Miller body is 89,530 bytes with an 8,128-byte
frame. The incumbent is 85,031 bytes with a 5,888-byte frame, and pinned VROOM
is 85,732 bytes with a 4,096-byte frame. The candidate is therefore about 4.4%
larger than VROOM and uses about 4 KiB more local stack, but remains bounded.

The candidate Miller section has no allocation, exception-runtime,
gather/scatter, selected-vector conversion, or out-of-line `_vv` reference.
Its only calls are fixed-size `memset` and `memcpy` lowering. The full candidate
pairing benchmark wrapper has a 7,552-byte frame versus 6,912 bytes for the
incumbent; the nested Miller frame remains the peak observed local stack.

## Native performance gate

Exactly one Google Benchmark process ran on Peach's AMD Ryzen 9 7950X. One
Clang 21.1.8 Release executable contained the incumbent, CatBLST MatrixNoK
candidate, and actual pinned VROOM implementation. Seven repetitions used a
minimum case time of 0.2 seconds and reported these means:

| Operation | Incumbent | CatBLST MatrixNoK | VROOM | vs incumbent | vs VROOM |
| --- | ---: | ---: | ---: | ---: | ---: |
| Miller loop | 69.454 us | 55.066 us | 58.131 us | -20.7% | -5.3% |
| Full pairing | 166.404 us | 141.775 us | 139.671 us | -14.8% | +1.5% |

Coefficients of variation were 0.42%, 0.14%, 0.38%, 0.16%, 0.02%, and 0.01%
for incumbent/candidate/VROOM Miller followed by incumbent/candidate/VROOM
pairing. CPU scaling was enabled, so these are selection measurements rather
than final release claims.

The candidate passes the agreed gate: it materially beats the incumbent and
the full pairing remains inside the 2% non-regression tolerance relative to
VROOM. Subtracting the Miller means leaves approximately 86.7 us for the
candidate inversion/final-exponent portion versus 81.5 us for VROOM. This is a
diagnostic decomposition, not an independently timed measurement, but it
identifies that boundary as the next optimization target.

Pinned VROOM's complete pairing depends on support symbols from its vendored
BLST object. The ordinary CatBLST benchmark target therefore enables only the
VROOM Miller comparison. For this one-off evidence build,
`BLST_VROOM_PAIRING_PREVIEW=1` enabled the VROOM pairing case, and a copied
VROOM `server_small.o` had only its colliding `blst_div_by_z` and
`blst_div_by_zz` symbols renamed with `objcopy`. This object was used solely to
link the benchmark comparison; no VROOM support object enters CatBLST or the
ordinary benchmark target. The guarded committed source was separately
compiled and linked after the measurement, and all six intended benchmark
names were discovered. No second benchmark process ran.

## Exact correspondence inputs

- `miller.h` SHA-256:
  `31c02b7af7e10448195bbae90e3d65181160eb67007a2e9708cf89572d993a97`
- `pairing.h` SHA-256:
  `8491cbca8e7c927a661fe338e201ae2a6492359e210ec17679e07b87914be7d7`
- `matrixnok_pairing.h` SHA-256:
  `f2d9cb9eb40e6cc1ca2192de6d6d3a54c990a94a72890a2f79ed42981e220936`
- `conversion_inversion.h` SHA-256:
  `a94d17ee0ec4829f283d9daeb2e3fb84b6055cb3d7e682b05caa45caf6bac075`
- `matrixnok_fp.h` SHA-256:
  `2a83c539980f7cf3219245920f5182501a46cfdd5c15fad8d7e3f31d915c5155`
- focused test SHA-256:
  `540acb591acdd2f9f4c6b1ba85ab09758f464ad851d01db382d5c8e6b631c9f0`
- CatBLST pairing benchmark SHA-256:
  `9a02d873f67507634e0d4fd842a87fc8970ea5050577fd296b641fc2df74b74c`
- VROOM preview benchmark SHA-256:
  `465ad10a9db4fc8dc8a051dbcd0ac35f676c5554e2fea26754eaf24e1a36cca3`

## Residual gates

- MatrixNoK remains isolated from production dispatch.
- The inversion/final-exponent portion should be reviewed against VROOM before
  freezing the implementation for proof.
- The candidate Miller's larger code and stack should be reconsidered only if
  the next integration profile exposes an instruction-cache or stack issue;
  current end-to-end timing is favorable.
- Formal arithmetic/bound proofs and exact C++/proof correspondence remain
  pending.
- Broad suites, sanitizers, fuzzing, Apple CI, CUDA, BTx/BTE, and portable-AVX
  Release were outside this focused packet.
