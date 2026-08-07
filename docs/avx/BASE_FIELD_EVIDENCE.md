# MatrixNoK Base-Field Candidate Evidence

## Packet boundary

- Candidate implementation commit: `f18c547bc972442ac30e987e06b47c79bf28fe88`
- Upstream algorithm commit: `d6ecd534c54a41702ea80ae75b52250b6d51b372`
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

Exact commit `f18c547bc972442ac30e987e06b47c79bf28fe88` was transferred as a Git
bundle and built on the AMD Ryzen 9 7950X Peach host with GCC 15.2.0, Release,
`-march=native`, native and portable AVX enabled, and assembly, CUDA, and
benchmarks disabled. The target was built with one job.

All 14 native-plus-portable compatibility tests passed: four MatrixNoK tests
and ten direct `_cc` versus `_mm512` intrinsic tests. The existing `AUD-102`
GCC ignored-vector-attribute warnings remain visible.

GCC also rejected two batch calls whose `int` template argument had been left
to deduction from a `size_t` `std::array` extent. Explicit fixed batch counts
repair that cross-compiler defect with no generated-code cost.

## Optimized code generation

`audit_avx_matrixnok_fp_multiply` was compiled with GCC 15.2.0 using
`-O3 -DNDEBUG -march=native`. Its arithmetic body contains:

- 23 `vpmadd52huq` and 25 `vpmadd52luq` instructions;
- seven fixed `vpermq` instructions;
- no candidate or `_vv` calls;
- no scalar-dependent branch or address;
- no vector conversion or lane-memory round trip; and
- no stack data load or store.

The normal hardened compiler configuration reserves a 192-byte logical frame
and emits a canary check because inlined `AVXVector` temporaries are
`std::array`-backed. Inspection-only compilation with
`-fno-stack-protector` emits no call, branch, or stack memory access. CatBLST
does not disable stack protection in this packet. Whether the canary overhead
matters is deferred to the approved serial benchmark gate.

Earlier probes caught and repaired an out-of-line candidate multiplication,
batch-of-one stack adapters, and a by-value aggregate return. The accepted hot
interface writes directly to caller-owned storage, supports exact output
aliasing with either input, and keeps batch APIs for actual batching.

## Exact correspondence inputs

- `matrixnok_fp.h` SHA-256:
  `857546fbfd27eeb4f72850398e7b61812a7b93156ca5f4687c0e63ab08bce822`
- `changebase.hpp` SHA-256:
  `75fcd7be6cb157e6e1ca676d8cecf702f3ae07fd899c4ba23d83568f071748af`
- `avx512ifma.hpp` SHA-256:
  `a2a23a7d7e7f921c0e243af192dce540928c07ffc64c166f7c1a6bc048b6aa82`

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
- Integrate the candidate into Fp2/Fp12 only after the base-field theorem is
  accepted.
- Run one approved serial performance comparison before production selection.

No benchmark, broad suite, sanitizer run, CUDA work, or production dispatch
change occurred in this packet.
