# MatrixNoK Fp2 Candidate Evidence

## Packet boundary

- Candidate implementation through: `19f37a1`
- Upstream algorithm commit: `d6ecd534c54a41702ea80ae75b52250b6d51b372`
- Compared VROOM head: `26bd9aee2168f167583cf096add355547a7b56d4`
- C++ owner: `blst/avx/fields/matrixnok_fp2.h`
- Base-field owner: `blst/avx/rns/matrixnok_fp.h`

This packet adds an isolated BLS12-381 Fp2 candidate. It does not change
production dispatch, the incumbent AVX representation, or a public API.
Elements represent `c0 + u*c1`, where `u^2 = -1`.

This document is historical performance and differential evidence for the
direct comparison-only owner `matrixnok_fp2.h`. The formal P3 proof instead
targets the selected production-shaped `FP2Ring<FpRingCandidate>` in
`fp2.h`; the two implementations have different delayed-reduction schedules
and must not share a correspondence claim merely because their field formulas
agree.

Multiplication retains VROOM's simple four-product formula but combines those
products into two wide output accumulators before reducing or changing base:

- real: `a*c - b*d`;
- imaginary: `a*d + b*c`.

Squaring issues three products and doubles the `a*b` accumulator:

- real: `a^2 - b^2`;
- imaginary: `2*a*b`.

Both operations are allocation-free, have fixed schedules, and support exact
output aliasing through wrappers that compute the complete result before the
assignment. The fixed BLS12-381 base and extension-ring owners are stateless,
empty C++ types. Each has the language-mandated minimum object size of one
byte and owns no parameter data or pointer.

## Correctness evidence

`AvxMatrixNoKPrimitives.Fp2MultiplyAndSquareMatchScalar` uses scalar
`blst::Fp2::mul` and `blst::Fp2::sqr` as the semantic oracle. It covers:

- zero, one, `p-1`, and dense fixed base-field components;
- all 64 pairings of eight Fp2 values;
- all eight one-shot squares;
- exact left, right, and square output aliasing;
- a 64-operation dependent multiplication chain;
- a 16-operation dependent squaring chain; and
- the `< 2m_i` postcondition in both RNS bases after every retained result.

The focused test passes under local MSVC 19.50 portable-AVX Debug and Peach
Clang 21 native AVX-512IFMA Release. The native build enables both native and
portable support, while the candidate selects native operations.

## Optimized code generation

Clang 21.1.8 compiled the isolated probes with C++23, `-O3`, `-DNDEBUG`, and
`-march=native` on Peach. Exact operation bodies report:

| Operation | Instructions | Frame | Calls | Conditional branches | IFMA low/high | Immediate rotations |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| multiply | 352 | 0 bytes | 0 | 0 | 54 / 50 | 14 |
| square | 348 | 0 bytes | 0 | 0 | 52 / 48 | 14 |

No `_vv` wrapper or candidate operation remains out of line. The formula,
reductions, and fixed MatrixNoK permutations are visible in one optimized
body. The stateless refactor removes the former 32-byte frame and 21
instructions from each Fp2 probe. Later Fp12 and pairing probes must still
recheck aggregate stack pressure.

## Clang 21 performance comparison

CatBLST and the actual pinned VROOM headers were compiled into separate
translation units of the same executable. One process ran the six Fp2 cases
with randomized interleaving, nine repetitions, and 0.2 seconds minimum per
repetition. Ring construction and operand conversion occurred outside the
timed loop.

Mean CPU time was:

| Operation | CatBLST prior 50-bit | CatBLST MatrixNoK | VROOM MatrixNoK | Candidate vs prior | Candidate vs VROOM |
| --- | ---: | ---: | ---: | ---: | ---: |
| multiply | 41.4 ns | 33.6 ns | 34.7 ns | -18.8% | -3.2% |
| square | 41.1 ns | 33.6 ns | 34.1 ns | -18.2% | -1.5% |

The candidate coefficients of variation were 0.05% for multiplication and
0.07% for squaring. This clears the performance gate for retaining the simple
four-product design as the field-tower baseline. A three-product Karatsuba
variant is not selected without separate correctness and timing evidence: its
wider sums and modular corrections may consume the saved IFMA work.

## Stateless-context resolution

Commit `19f37a1` replaces the 5120-byte runtime-composed base context with
empty fixed-parameter operation types. Reducer constants, shifted correction
vectors, and matrices are compile-time-derived, immutable, cache-line-aligned
data. There is no heap allocation, mutable global state, function-local
initialization guard, or runtime constructor.

Focused compile-time checks require both ring types to be empty and trivially
constructible/copyable. Runtime tests compare the fixed helpers directly with
the former generic matrix and reducer helpers, including single and batch
conversion and a bound-24 reduction. Native symbol inspection finds no
`__cxx_global_var_init`, `_GLOBAL__sub_I`, or `__cxa_guard` symbol. `AUD-214`
is therefore validated.

## Exact correspondence inputs

- `matrixnok_fp.h` SHA-256:
  `3070239f8369e2e397026cdb5731b7af33a0b4d01a2dfccfd719cea0ca32e4bc`
- `matrixnok_fixed.h` SHA-256:
  `5938d841014795e97715107ab8a1345e036be2323755890c7377e03b63309a2d`
- `matrixnok_fp2.h` SHA-256:
  `ec900342d50040d0268660bcfa4efea3251a25a725d012990f1373aa787f0e31`
- `matrixnok_change_base.h` SHA-256:
  `74eec0fc27747cb7c27ee4d062ce66c6f601d2283d267659e4a586162de0a032`
- `generated_matrixnok_constants.hpp` SHA-256 at the measured selection commit:
  `d1bd0045e8d7f223cb332842bd4c9023daefc715079c0649fbce953cabdfea43`

The current repaired constant hash is recorded in
`proof/avx/CORRESPONDENCE.csv`; the eight ConvertTo estimate changes do not
alter Fp2 arithmetic or this historical performance measurement.

Any change to these sources reopens the code-generation, performance, and
proof correspondence.

## Residual gates

- The direct comparison-only microkernel retains its historical bound-8 and
  bound-24 proof obligations if it is ever reconsidered for selection.
- The selected generic `FP2Ring<FpRingCandidate>` formulas, tight bound-8
  accumulators, exact bound-24 combined values, and conservative
  `800 + 32 = 832` reducer branch are now modeled in
  `proof/avx/CatBLST/AVX/MatrixNoK/Fp2.lean` and bound to current source hashes
  by the in-review P3 correspondence row.
- The variadic strided layout and independent per-item batch R1 recurrence are
  now checked. A valid split-`M1` Fp2 input formally demonstrates that the
  correction-free R1 representative can be `-1`; the existing R2 schedule is
  also proved correct for that case. Complete the physical delayed-M1 signed
  composition and prove the R2 fixed-point margin across every reachable
  negative representative before accepting P3.
- Extend the isolated tower to Fp12 and measure complete pairing workloads
  before changing production dispatch.

No broad build, sanitizer run, CUDA work, BTx/BTE work, portable-AVX Release
build, or production dispatch change occurred in this packet.
