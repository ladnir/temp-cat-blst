# MatrixNoK Fp12 and GT Candidate Evidence

## Packet boundary

- Candidate implementation: `1011916baf69a4ecb4b0e012c30a74051c8ec850`
- Sparse/final-exponent follow-up: `3e1f99a`
- Cyclotomic split follow-up: `a1c8f9f652b5ae7d8e77966c454f312e3079bd95`
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

## Cyclotomic split follow-up

Review of inversion and hard final exponentiation found that CatBLST and pinned
VROOM use the same exponentiation chain. A focused diagnostic measured the
candidate inversion at 2.794 us and the final exponent after inversion at
77.926 us, making the cyclotomic schedule and its generated code the material
part of the remaining pairing delta.

The original MatrixNoK cyclotomic body used 9,498 bytes of text, a 4,352-byte
frame, 279 `vmovdqa64` instructions, and 235 stack-pointer references. The
selected schedule preserves the exact formulas and 258 high plus 292 low IFMA
operations, but emits two fixed noinline halves:

| Half | Text | Frame | `vmovdqa64` | stack references |
| --- | ---: | ---: | ---: | ---: |
| first | 4,785 bytes | 1,408 bytes | 117 | 59 |
| second | 3,938 bytes | 512 bytes | 88 | 23 |

Neither helper contains a call relocation, gather/scatter, allocation, or
exception-runtime reference. The incumbent remains on its original monolithic
source path and its representative optimized symbol remains exactly `0x2844`
bytes. A source-sharing cleanup was rejected because it grew that incumbent
symbol by about 7%.

The fresh split was selected from complete-pairing timing, where the final
exponentiation is exercised repeatedly in its real call graph. The standalone
cyclotomic benchmark object in that process predated the header change, so its
207 ns result remains evidence for the old monolithic candidate only and is
not attributed to the split implementation.

## Rejected alternatives

- Splitting correction operations into a separate batch phase worsened Fp12
  multiplication from approximately 303 ns to 323 ns.
- The existing version-2 batch schedule compiled and passed the focused native
  tests, but grew the body to 10,630 bytes and increased `vmovdqa64` and stack
  references to 423 and 341 without reducing IFMA work. It was rejected before
  timing.
- Sharing the split helper bodies with the incumbent grew incumbent code by
  about 7%; deliberate formula-local duplication retains unchanged incumbent
  code generation.
- Caching cyclotomic constants in local references worsened the specialized
  square without helping multiplication.
- Storing constants behind static reference-returning accessors reduced one
  apparent frame but regressed the measured cyclotomic path.

All rejected arithmetic variants were removed. The retained multiply reads
the wide offset once per operation; the cyclotomic path keeps its direct
compile-time constant uses.

## Exact correspondence inputs

- `fp12.h` SHA-256:
  `f91f32711620f47988dddef464266bea595e620bc4fbcf5e8c8ca41432762ae2`
- `matrixnok_fp12.h` SHA-256:
  `ec3963f55e2e792fd4dbc149d6b4b6d94d6de892926a1cbed49350be6938d4f8`
- `matrixnok_fixed.h` SHA-256:
  `73ba31d8a347208ac84173790a61929a34f443a0d47783cffa21559661ac026f`
- `matrixnok_fp.h` SHA-256:
  `11c172dccd76dba2d923ea205201b32a1e6bdf679f63f3d1560016d4dcc017d7`
- `generated_matrixnok_constants.hpp` SHA-256 at the measured selection commit:
  `d1bd0045e8d7f223cb332842bd4c9023daefc715079c0649fbce953cabdfea43`

The current repaired constant hash is recorded in
`proof/avx/CORRESPONDENCE.csv`; the eight ConvertTo estimate changes do not
alter Fp12 arithmetic or this historical performance measurement.
- `glv4.h` SHA-256:
  `8765c031be885fed24b8743c7ee677ca02640c3f60423bd9d8ae53fa87cc1d31`

## Formal correspondence

Rows `AVX-P3-002` through `AVX-P3-005` bind the current source hashes rather
than the historical performance-selection hashes above. Proof commit `256b980`
completes the selected split cyclotomic square: all twelve literal output
polynomials equal the general Fp12 square under exactly the scalar
`is_cyclotomic` Frobenius relation. It also checks both generated fraction
constants, all 32 RNS lanes, the generated offset's semantic zero, and the
first/second split physical maxima of 28/888 and 101/2041 M1/RNS units against
the 800/5760 reduction boundaries. The theorem intentionally accepts zero and
does not replace the cyclotomic precondition with a narrower GT assumption.

Proof commit `21887f9` completes `AVX-P3-005`. It proves the concrete flat
conjugation is the general Fp12 `p^6` map and transcribes the full final-
exponentiation schedule over its precise nonzero-unit domain. Exact integer
identities verify the negative BLS parameter, the 2/3/9/32/15 addition chain,
and the repository's intentional cubed final exponent. The source ledger is
35 multiplications, 315 cyclotomic squares, nine conjugations, Frobenius counts
1/2/1, and one inversion. Existing
`AvxMatrixNoKPrimitives.Fp12SparseMultiplyAndFinalExpMatchScalar` and complete-
pairing differential evidence cover the concrete inversion provider; that
provider is explicitly outside the new algebra theorem.

Proof commit `82f26c1` completes `AVX-P3-006` for public GT multiplication. It
proves the exact canonical unsigned radix-`z` split and odd-Frobenius
conjugation basis, the literal sixteen-entry table, fixed 64-round CT executor,
unsigned VT skip executor, binary CT fallback, fixed 256-bit basic VT path, and
public scalar-mask/dispatch composition. The CT ledger is 76 total Fp12
multiplications including table construction and final normalization, 64
cyclotomic squares, 960 table-candidate comparisons, and 11,520 coefficient
selections. The proof domain is GT and composes the accepted general-Fp12,
Frobenius, conjugation, and cyclotomic results rather than broadening them.

## Residual gates

- Production dispatch remains unchanged. Retained G2, inversion, Miller-loop,
  and complete-pairing correctness and selection gates now pass, but proof and
  release evidence remain outstanding.
- The split closes the measured complete-pairing performance gap. A fresh
  standalone split-cyclotomic measurement and stable release-condition timing
  remain future evidence, not prerequisites for beginning the proof.
- General Fp12 multiplication/squaring, Frobenius maps, cyclotomic square,
  final exponentiation, public GT exponentiation, and the balanced compiled
  FFT exponentiation path are accepted by `AVX-P3-002` through `AVX-P3-007`.
  Proof commits `51462b9` and `effdb3d` establish the exact signed split,
  bytecode compiler/executor, 63-byte valid-input bound inside the 70-byte
  owner, 62-square/62-executor-product maxima, and canonical-`Fr` FFT
  corollary. The C++ `Fr::from_mont` representation bridge and P0
  compiler/ISA lowering remain explicit residual boundaries.
- Sanitizers, fuzzing, broad suites, Apple CI, CUDA, BTx/BTE, and portable-AVX
  Release were outside this focused packet.
