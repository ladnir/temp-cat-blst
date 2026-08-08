# AVX-512 Modernization and Proof Project

## Objective

Replace CatBLST's older VROOM-derived AVX-512 arithmetic with the faster
BLS12-381 MatrixNoK design, preserve CatBLST's public contracts and hardening,
and produce both native execution evidence and a machine-checked mathematical
proof for the retained arithmetic kernels.

This is a separate Phase 16 project. The audited scalar backend is the
correctness oracle and is not refactored to accommodate the candidate.

## Source pins

- CatBLST starting commit: `df92cd496a33e956463c4afdb1bd3567bcb55d5a`
- Original VROOM import: `31e98bbf7f32e1345299d4276dd6a5bc44dc0b69`
- MatrixNoK algorithm commit: `d6ecd534c54a41702ea80ae75b52250b6d51b372`
- Reviewed VROOM head: `26bd9aee2168f167583cf096add355547a7b56d4`
- Upstream repository: <https://github.com/SimonLangowski/VROOM>

The 31-commit delta is classified in `VROOM_UPSTREAM_MANIFEST.md`, and the
selected dependency graph is recorded in `MATRIXNOK_PORT_MAP.md`. Only
`d6ecd534` changes executable CPU arithmetic. CatBLST will translate the
BLS12-381 subset; it will not merge VROOM wholesale.

## Scope

Primary scope:

- AVX-512IFMA lane primitives used by the selected backend;
- BLS12-381 MatrixNoK parameters and fixed permutation tables;
- RNS conversion, Montgomery multiplication, reduction, and bounds;
- Fp, Fp2, Fp12, curve, Miller-loop, final-exponentiation, and pairing code
  needed by CatBLST's retained AVX API; and
- compile-time selection, differential tests, native code generation, and
  formal specifications for those operations.

Explicitly excluded:

- BN254;
- GPU/CUDA code;
- VROOM's vendored BLST copy, examples, baselines, benchmarks, result logs,
  artifact scripts, and dependency installers;
- runtime GMP or runtime parameter generation in production; and
- changes to frozen BTx/BTE.

## Implementation boundary

The candidate is introduced side-by-side behind an internal build option until
it passes correctness and performance gates. Existing CatBLST APIs and semantic
types remain the integration boundary. Upstream fallback intrinsics are
translated through CatBLST's `_cc512_*` emulation and `_vv512_*` selected-vector
layer; native code must call `_mm512_*` without conversion objects or lane
round-trips.

Precomputed constants are checked in as fixed compile-time data with provenance
and reproducible validation. Production initialization performs no search,
allocation, GMP computation, or mutable global setup.

Focused evidence for the prerequisite vector and bounded-reduction packet is
recorded in `PRIMITIVE_EVIDENCE.md`.
The first fixed-parameter import and its independent regeneration evidence are
recorded in `CONSTANT_EVIDENCE.md`.
The isolated 8x8 MatrixNoK change-of-base kernel is recorded in
`CHANGE_BASE_EVIDENCE.md`.
The first complete side-by-side base-field multiplication candidate and its
portable/native/code-generation evidence are recorded in
`BASE_FIELD_EVIDENCE.md`.
The isolated Fp2 multiplication and squaring candidate, same-process VROOM
comparison, and stateless fixed-context resolution are recorded in
`FP2_EVIDENCE.md`.
The end-to-end G1 and GT FFT previews, including the corrected G1 input gate,
are recorded in `FFT_PREVIEW_EVIDENCE.md`.
The production-shaped MatrixNoK G1 candidate, independent scalar oracle,
native code generation, and same-process incumbent/VROOM comparison are
recorded in `G1_EVIDENCE.md`.
The isolated MatrixNoK G2 candidate, scalar differential coverage, shared
signed-expansion repair, and preliminary native code generation are recorded
in `G2_EVIDENCE.md`. Its correctness and performance gates pass; direct final
table construction closes `AUD-218` by reducing peak CT GLS stack from about
102 KiB to about 55 KiB without a performance regression.
The production-shaped Fp12/GT candidate, parameter-bound constant repair,
native code generation, and same-process incumbent/VROOM comparison are
recorded in `FP12_EVIDENCE.md`.
The isolated MatrixNoK Miller-loop and full-pairing preview, scalar
differentials, native code generation, and serial incumbent/VROOM comparison
are recorded in `PAIRING_PREVIEW_EVIDENCE.md`. Miller is 20.7% faster than the
incumbent and 5.3% faster than VROOM. After splitting the register-heavy
MatrixNoK cyclotomic square, full pairing is 21.2% faster than the incumbent
and 5.5% faster than VROOM. Performance selection is complete; the candidate
remains isolated while formal proof and release gates are pending.

## Work packets

1. **Manifest and mapping** — classify every upstream commit and map each
   selected upstream file to its CatBLST owner.
2. **Primitive and RNS candidate** — import MatrixNoK constants, permutation,
   conversion, multiplication, and reduction behind an internal candidate
   option.
3. **Field tower** — port Fp2/Fp12 operations and validate each operation
   against scalar CatBLST.
4. **Curve and pairing** — port curve formulas, Miller loop, final
   exponentiation, and pairing while preserving documented subgroup and
   validation contracts.
5. **Selection** — inspect native code generation and run one approved serial
   benchmark gate on Peach. Promote the candidate only if it is correct and
   materially faster without unacceptable regression.
6. **Formal proof** — prove the selected arithmetic and bounds described in
   `PROOF_PLAN.md`, using `MATRIXNOK_CONSTRUCTION.md` as the readable
   construction map; bind exact constants and source hashes into the evidence.
7. **Release evidence** — complete native, emulated, sanitizer, fuzz,
   constant-time, API, and audit records.

## Validation profiles

- Local portable-AVX Debug: fast deterministic differential and invariant
  testing; never a performance claim.
- Peach native AVX-512IFMA: Clang 21 is the authoritative code-generation and
  serial performance compiler; GCC remains cross-compiler correctness evidence.
- Scalar CatBLST: independent semantic oracle at the same public contracts.
- Formal model: exact integer/RNS theorems, independent of compiler behavior.

No full build, sanitizer sweep, benchmark, or multi-process remote job is a
per-commit default. Resource-intensive gates require explicit approval and run
serially.

## Completion criteria

The project completes only when the selected BLS12-381 backend has:

- no unclassified VROOM delta;
- operation-complete scalar differential coverage;
- native AVX-512IFMA execution and accepted optimized code generation;
- accepted performance evidence from one benchmark process at a time;
- proved arithmetic correctness and overflow/bound invariants for the retained
  primitive/RNS/field kernels;
- explicit residual assumptions for intrinsics, compiler lowering, and the
  hand-maintained C++/proof correspondence; and
- updated Phase 16 file reviews and release records at exact commits.
