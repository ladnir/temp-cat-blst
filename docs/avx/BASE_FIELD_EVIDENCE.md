# MatrixNoK Base-Field Candidate Evidence

## Packet boundary

- Candidate implementation through: `19f37a1`
- Upstream algorithm commit: `d6ecd534c54a41702ea80ae75b52250b6d51b372`
- Compared VROOM head: `26bd9aee2168f167583cf096add355547a7b56d4`
- C++ owner: `blst/avx/rns/matrixnok_fp.h`
- Fixed parameters: `blst/avx/generated_matrixnok_constants.hpp`

This packet adds an isolated BLS12-381 base-field multiplication candidate.
It has no production caller and does not change CatBLST dispatch, public types,
or the incumbent AVX field representation.

The candidate uses two eight-lane, 50-bit RNS bases. Multiplication reduces the
first-base product, accumulates the fixed MatrixNoK r1 permutation into the
second-base product, reduces the combined value, and expands the result back
to both bases. The fixed source bounds are:

- standard residues: `[0, 2m_i)`;
- elementwise products: bound 4;
- no-correction r1 contribution: bound 16; and
- combined second-base accumulator: bound 20, below the configured bound-40
  reduction capacity.

These source bounds are checked by C++ types and assertions but remain formal
proof obligations until the corresponding Lean theorems are complete.

## Differential oracle

`AvxMatrixNoKPrimitives.BaseFieldMultiplyMatchesScalar` converts candidate
outputs back through the fixed convert-from matrix, then decodes them with
explicit portable scalar Montgomery primitives. It deliberately does not
compare raw RNS lanes: the no-correction representation may encode the same
field value with noncanonical residues.

Coverage includes:

- `0`, `1`, `2`, `p-1`, dense fixed canonical values, and exact output aliasing
  with either input;
- every pair from eight boundary/fixed values;
- 64 deterministic dense random canonical pairs;
- a 64-step dependent multiplication chain; and
- the declared `< 2m_i` postcondition for both RNS bases after every operation.

The scalar oracle uses `cc_from_mont_n`, `cc_redc_mont_n`, and
`cc_mul_mont_n` directly, so a combined native-plus-portable configuration does
not silently compare the candidate with an assembly-dispatched copy of itself.

## Local portable Debug

MSVC 19.50.35728 built only `blst_avx512_compat_tests`, with one job, in a
portable-only Debug configuration. Native AVX-512, assembly, CUDA, and
benchmarks were disabled. All four MatrixNoK tests and the independent
constant validator passed. No portable-AVX Release build was run.

## Peach native AVX-512IFMA

The AMD Ryzen 9 7950X Peach host built the focused native target serially with
Clang 21.1.8, Release, and `-march=native`. CUDA, assembly, broad tests, and
unrelated benchmarks were disabled. The combined native-plus-portable
compatibility target passes the MatrixNoK scalar differential tests at batch
sizes 1, 2, 4, 8, and 12. The final alignment packet additionally passes the
cache-line-alignment and complete base-field multiplication tests 2/2.

Earlier GCC 15.2.0 execution remains useful cross-compiler correctness
evidence, including its rejection of the previously invalid `int`/`size_t`
batch-template deduction. GCC measurements are not used for performance
selection.

## Optimized code generation

The selected Clang 21 code uses VROOM's fixed unrolled rotation schedule. Each
r1 permutation has four `valignq` and three `vshufi64x2` immediate rotations,
not variable-index `vpermq`; `_vv` wrappers are fully inlined and no conversion
objects or lane round trips survive. The batch-8 body has the same 144
`vpmadd52luq`, 128 `vpmadd52huq`, 32 `valignq`, and 24 `vshufi64x2`
instruction counts as the pinned VROOM reference.

The first comparison exposed that CatBLST's generated 8x8 permutation table
landed at address 48 modulo 64, making every 64-byte row load cross a cache
line. Commit `cd28256` gives the table explicit 64-byte alignment and makes
both the independent constant validator and native test enforce that layout.
The linked table now lands at zero modulo 64 and Clang emits aligned row loads.

Earlier probes also caught and repaired an out-of-line candidate
multiplication, batch-of-one stack adapters, and ambiguous aggregate handling.
The retained returning operation is fully inlined; the output-reference
wrapper delegates to it and preserves exact aliasing with either input.

## Clang 21 performance comparison

Commit `47124e0` adds an optional benchmark-only reference translation unit to
the existing `avx_group_ops_benchmarks` executable. Configuration rejects any
VROOM source tree except commit `26bd9aee2168f167583cf096add355547a7b56d4`.
CatBLST and the actual VROOM headers therefore run in the same Clang 21.1.8
binary, process, Google Benchmark harness, core, and randomly interleaved run.
Both use the same fixed field operands. One benchmark process ran at a time.

The final nine-repetition, 0.2-second run reports mean CPU time:

| Operation | CatBLST | VROOM | CatBLST delta |
| --- | ---: | ---: | ---: |
| multiply | 19.098 ns | 19.073 ns | +0.133% |
| batch 1 | 19.106 ns | 19.059 ns | +0.251% |
| batch 2 | 31.315 ns | 31.479 ns | -0.522% |
| batch 4 | 52.984 ns | 53.082 ns | -0.184% |
| batch 8 | 104.812 ns | 102.501 ns | +2.254% |
| batch 12 | 175.312 ns | 187.689 ns | -6.595% |

The six-case geometric mean favors CatBLST by 0.816%. A separate 15-repetition,
0.3-second batch-8 confirmation measured 104.393 ns versus 102.548 ns
(+1.799%). Batch 8 is therefore recorded as a borderline compiler-layout
microcase rather than hidden by averaging. It does not justify further source
complexity before field-tower measurements show whether that arity is a
production bottleneck.

Against CatBLST's prior 50-bit field path, the same Clang profile measured the
candidate at 19.1 ns versus 22.8 ns for one multiplication, with wins at every
measured batch size. The older 52-bit path measured 24.5 ns. GCC diagnostic
runs were discarded and are not part of this selection evidence.

The final stateless-context run measured 19.4 ns for CatBLST MatrixNoK, 19.1
ns for pinned VROOM, and 22.8 ns for the prior CatBLST 50-bit path. Removing
the runtime-owned context therefore costs about 1.6% against the earlier
context-backed candidate, within the 2% primary-case gate, while retaining a
14.9% improvement over the incumbent. The simpler ownership and initialization
model is accepted on that evidence.

## Exact correspondence inputs

- `matrixnok_fp.h` SHA-256:
  `3070239f8369e2e397026cdb5731b7af33a0b4d01a2dfccfd719cea0ca32e4bc`
- `matrixnok_fixed.h` SHA-256:
  `5938d841014795e97715107ab8a1345e036be2323755890c7377e03b63309a2d`
- `matrixnok_change_base.h` SHA-256:
  `74eec0fc27747cb7c27ee4d062ce66c6f601d2283d267659e4a586162de0a032`
- `changebase.hpp` SHA-256:
  `75fcd7be6cb157e6e1ca676d8cecf702f3ae07fd899c4ba23d83568f071748af`
- `avx512ifma.hpp` SHA-256:
  `9da2bd1b3fc11be3a33cdf4b336206233e9580225cc7f92425bd8737b87c64b1`
- `generated_matrixnok_constants.hpp` SHA-256:
  `3b4349e9e9afd763b1297d533e84dcb2ad9c04a92d04bfa557a48a67a3b34302`

The fixed-constant hash and regeneration record remain in
`CONSTANT_EVIDENCE.md`. Any change to these sources or constants reopens the
code-generation and proof correspondence.

## Residual proof and selection gates

- Prove lane extraction and IFMA semantics for the exact operations used.
- Prove the fixed conversion matrices and both Montgomery reductions.
- Prove that omission of the r1 quotient correction preserves the represented
  value under the bound-20 input established by multiplication.
- Prove the `< 2m_i` result bounds and absence of unintended 64-bit IFMA
  accumulator wrap.
- Bind the proved declarations to the source and constant hashes above.
- The isolated Fp2 candidate now passes its operation-level gate; extend the
  candidate through Fp12 and pairing-level workloads before production
  selection.
- Revisit the batch-8 microcase only if field-tower profiles show it on a
  primary hot path.

No broad suite, sanitizer run, CUDA work, or production dispatch change
occurred in this packet. Formal proof remains deferred until the field-tower
candidate is selected on correctness and performance evidence.
