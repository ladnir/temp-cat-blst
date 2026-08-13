# MatrixNoK Base-Field Candidate Evidence

## Packet boundary

- Candidate implementation through: `adf2ed5df84de1e4bf705561e52d7941f04833a4`
- Retained-owner closure: `01bcec83fa45d80ee4bb39bdad2f8da0b96c5dd7`
- Upstream algorithm commit: `d6ecd534c54a41702ea80ae75b52250b6d51b372`
- Compared VROOM head: `26bd9aee2168f167583cf096add355547a7b56d4`
- C++ owner: `blst/avx/rns/matrixnok_fp.h`
- Fixed parameters: `blst/avx/generated_matrixnok_constants.hpp`

This document began as the isolated BLS12-381 base-field candidate packet.
The construction is now the selected production AVX field representation.
Commit `01bcec8` closes the retained-owner dispatch gap: ordinary one-shot and
fixed-batch `modmul` operations call the same exact-bound `multiply` kernels
that were compared with VROOM. Higher formulas still use delayed reduction
when they intentionally aggregate several products.

The selected field uses two eight-lane, 50-bit RNS bases. Multiplication reduces the
first-base product, accumulates the fixed MatrixNoK r1 permutation into the
second-base product, reduces the combined value, and expands the result back
to both bases. The fixed source bounds are:

- standard residues: `[0, 2m_i)`;
- elementwise products: bound 4;
- no-correction r1 contribution: bound 16; and
- combined second-base accumulator: bound 20, below the configured bound-40
  reduction capacity.

These source bounds are checked by C++ types and assertions and are now also
proved for the selected multiplication pipeline in Lean.

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

The native-MSVC compatibility follow-up at `33f6a8b` removes a separate source
of out-of-line work in the quotient accumulator. Runtime 64-by-64 products now
use `_umul128`, while constant evaluation retains portable limb arithmetic.
The optimized production `modmul` body contains eight direct `mul`
instructions, zero wide-integer calls, and zero conditional jumps; the prior
body contained eight calls to the custom compatibility `operator*`. This is
code-generation evidence only; no new benchmark ran.

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
  `4a09c34fb4fe97ad23ced3e4e3736290d0912171b2205b16c6d7066121d6dd3f`
- `matrixnok_fixed.h` SHA-256:
  `0ca412ccef04a3aa45c87826fb25d02e61697404728e865de5395ecf24583c04`
- `matrixnok_change_base.h` SHA-256:
  `d706b60bfbbc3d77338bbf4b3015e3cdebdc30af48728d96c395922b51107d88`
- Historical generic `changebase.hpp` comparison owner SHA-256:
  `75fcd7be6cb157e6e1ca676d8cecf702f3ae07fd899c4ba23d83568f071748af`
  (removed at `b1036f7`; the selected fixed owner was unchanged and its test
  oracle is now an independent scalar IFMA-accumulator model)
- `avx512ifma.hpp` SHA-256:
  `d7dda9a8054fec44d3c16980bee316f57084729c14ec1927fadb50671ae53b28`
- `generated_matrixnok_constants.hpp` SHA-256:
  `1643256db5038b0a643ded8053dc833ab2613d20dbc6de8122296f0a4699361a`

The fixed-constant hash and regeneration record remain in
`CONSTANT_EVIDENCE.md`. Any change to these sources or constants reopens the
code-generation and proof correspondence.

The selected-vector owner was reaccepted at `6cc9e753`. Its template now
requires the exact eight lanes modeled by `CatBLST.AVX.Vec8`; no retained
production or test instantiation used another width. Focused portable checks
and optimized native selection, equality, and IFMA inspection passed. The
arithmetic bodies mapped by `AVX-P0-001` retain their established semantics.

The proof-preparation pass at `adf2ed5` names these bounds and dimensions,
documents accumulator and scratch mutation semantics at their owners, and
adds `MATRIXNOK_CONSTRUCTION.md` as the readable construction map. A proposed
narrowing of the delayed-path C++ type allowance from 32 to the exact R1 bound
16 changed optimized higher-tower scheduling, so it was not selected without a
new benchmark gate. Final Clang 21 probe assembly is byte-identical to the
performance-selected `e6b2995` baseline. The clean native-plus-portable focused
target passes all 17 MatrixNoK tests; local portable Debug passes the four
directly affected RNS/base-field tests.

## Formal result and residual selection gates

The pinned Lean workspace proves the selected Montgomery reducers, exact
cyclic R1 schedule, correction-free R1 relation, bound-20 fused M2 stage,
shifted R2 expansion, fixed quotient, accumulator no-wrap conditions, and the
`selected_base_field_pipeline` for nonnegative raw lifts. The theorem accepts the actual
standard lane domain `<2m_i`, returns both bases below `<2m_i`, and bounds the
raw QR output representative below `18p` within the C++ `40p` allowance.

Canonical import is mapped separately. Formal analysis found `AUD-221`: the
inherited ceiling estimate could select a negative lift outside the next
conversion's contract. Commit `5f0399e` rounds only ConvertTo estimates down.
Lean now proves that a canonical import selects a representative in `[0,2p)`,
preserves the field value, and satisfies the existing R2 input contract. The
exact former witness is proved to select its positive complement. Follow-up
commit `f91268d` maps the actual split IFMA streams and selected reducer calls
line by line. It proves canonical ConvertTo below `16m²` while retaining the
implementation's `800m²` reducer instantiation, derives both imported residue
bases, and closes the reachable import/multiply/export theorem without assumed
intermediate postconditions. The P2 C++ correspondence rows are accepted.
Lean still does not prove compiler lowering, the x86 ISA, or the scalar
byte-weaving/REDC transcription. The later `AUD-259`/`AUD-263` packet
independently proves the complete machine envelope and confirms branchless
MSVC carry lowering; scalar REDC semantics, IFMA semantics, and native Clang
lowering remain external.

The abstract negative-R2 counterexample remains useful as a boundary theorem:
R2 is not promised for arbitrary negative representatives. Repaired canonical
import no longer reaches that domain. `AUD-221` is validated: all 18 focused
tests pass locally and under Peach Clang 21 native AVX-512IFMA, and the import
probe retains the same 209-instruction shape with only eight changed immediate
constants.
The later Fp2, Fp12, curve, GT, and pairing integration and proof records
supersede the original candidate-only next steps. The retained batch-8
microcase remains a profiling-driven optimization question, not a correctness
gate. Commit `01bcec8` additionally makes all raw proof-domain and alias
preconditions executable in Debug and confirms the production `modmul`
spelling has no NDEBUG wrapper branch. No new benchmark was required because
the routed kernel is the already selected and measured implementation.
