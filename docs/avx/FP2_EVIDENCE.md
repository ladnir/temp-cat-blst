# MatrixNoK Fp2 Candidate Evidence

## Packet boundary

- Candidate implementation: `1ca9af2`
- Upstream algorithm commit: `d6ecd534c54a41702ea80ae75b52250b6d51b372`
- Compared VROOM head: `26bd9aee2168f167583cf096add355547a7b56d4`
- C++ owner: `blst/avx/fields/matrixnok_fp2.h`
- Base-field owner: `blst/avx/rns/matrixnok_fp.h`

This packet adds an isolated BLS12-381 Fp2 candidate. It does not change
production dispatch, the incumbent AVX representation, or a public API.
Elements represent `c0 + u*c1`, where `u^2 = -1`.

Multiplication retains VROOM's simple four-product formula but combines those
products into two wide output accumulators before reducing or changing base:

- real: `a*c - b*d`;
- imaginary: `a*d + b*c`.

Squaring issues three products and doubles the `a*b` accumulator:

- real: `a^2 - b^2`;
- imaginary: `2*a*b`.

Both operations are allocation-free, have fixed schedules, and support exact
output aliasing through wrappers that compute the complete result before the
assignment. The extension-ring object contains only an eight-byte pointer to
its owning base ring.

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
| multiply | 373 | 32 bytes | 0 | 0 | 54 / 50 | 14 |
| square | 369 | 32 bytes | 0 | 0 | 52 / 48 | 14 |

No `_vv` wrapper or candidate operation remains out of line. The formula,
reductions, and fixed MatrixNoK permutations are visible in one optimized
body. The small frame is accepted for this isolated candidate; later Fp12 and
pairing probes must recheck aggregate stack pressure.

## Clang 21 performance comparison

CatBLST and the actual pinned VROOM headers were compiled into separate
translation units of the same executable. One process ran the six Fp2 cases
with randomized interleaving, nine repetitions, and 0.2 seconds minimum per
repetition. Ring construction and operand conversion occurred outside the
timed loop.

Mean CPU time was:

| Operation | CatBLST prior 50-bit | CatBLST MatrixNoK | VROOM MatrixNoK | Candidate vs prior | Candidate vs VROOM |
| --- | ---: | ---: | ---: | ---: | ---: |
| multiply | 41.6 ns | 34.2 ns | 34.7 ns | -17.8% | -1.4% |
| square | 41.2 ns | 33.5 ns | 34.2 ns | -18.7% | -2.0% |

The candidate coefficients of variation were 0.11% for multiplication and
0.08% for squaring. This clears the performance gate for retaining the simple
four-product design as the field-tower baseline. A three-product Karatsuba
variant is not selected without separate correctness and timing evidence: its
wider sums and modular corrections may consume the saved IFMA work.

## Context-size finding

The focused test records `sizeof(FpRingCandidate) == 5120`, while
`sizeof(Fp2RingCandidate) == 8`. The base object currently owns copies of
fixed RNS matrices and two reducer contexts assembled by its constructor.
Construction was correctly outside the arithmetic benchmark, but a naive
production integration could introduce excessive object/stack size and
repeated fixed-parameter setup.

`AUD-214` therefore blocks production selection. The retained arithmetic can
continue to inform field-tower design, but the base owner must become
stateless or refer to one explicitly owned immutable fixed-parameter context.
The repair must not add hidden allocation, mutable global state, function-local
initialization guards, or per-operation construction.

## Exact correspondence inputs

- `matrixnok_fp.h` SHA-256:
  `556e108b6b726c407d7e2b68b0517e5b1f422c114c123bf995c761915f456423`
- `matrixnok_fp2.h` SHA-256:
  `cc0c058afea7230de2a3e9546c6a3c37db595b02f87c1c15cb4758361e7c544f`
- `matrixnok_change_base.h` SHA-256:
  `74eec0fc27747cb7c27ee4d062ce66c6f601d2283d267659e4a586162de0a032`
- `generated_matrixnok_constants.hpp` SHA-256:
  `f8b2047be64327830ecd596c87d5fe69e6af2591df3fad311cddbd9ff92f67be`

Any change to these sources reopens the code-generation, performance, and
proof correspondence.

## Residual gates

- Repair and remeasure the 5120-byte fixed-parameter owner before production
  integration.
- Prove both Fp2 formulas over the proved base-field representation.
- Prove the bound-8 wide accumulators, bound-24 combined second-base values,
  and absence of unintended 64-bit IFMA accumulator wrap.
- Bind the eventual Lean declarations to the exact retained C++ and generated
  constant hashes.
- Extend the isolated tower to Fp12 and measure complete pairing workloads
  before changing production dispatch.

No broad build, sanitizer run, CUDA work, BTx/BTE work, portable-AVX Release
build, or production dispatch change occurred in this packet.
