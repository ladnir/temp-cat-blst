# MatrixNoK G2 Candidate Evidence

## Packet boundary

- Candidate implementation: `2fc0ae98f962f5ab2259ed2d368f6328713b7310`
- Upstream algorithm commit: `d6ecd534c54a41702ea80ae75b52250b6d51b372`
- Compared VROOM head: `26bd9aee2168f167583cf096add355547a7b56d4`
- C++ curve owner: `blst/avx/curve/matrixnok_g2.h`
- Shared complete-add owner: `blst/avx/curve/matrixnok_curve.h`
- RNS/field owners: `blst/avx/{rns/matrixnok_fp.h,fields/fp2.h}`

This packet adds an isolated production-shaped G2 candidate. It changes no
public type, symbol, or production dispatch. The candidate applies CatBLST's
complete homogeneous schedule and existing GLS executor to the fixed,
stateless MatrixNoK Fp/Fp2 owners.

## Correctness design

The candidate converts only at the test boundary. Scalar CatBLST is the
independent oracle. Focused native coverage includes:

- generator and infinity conversion;
- complete add, double, equal inputs, opposite inputs, both infinity operand
  positions, and a 16-step deterministic addition chain; and
- constant- and variable-time GLS multiplication on the generator and its
  double for zero, one, two, `r-1`, fixed dense values, and 16 deterministic
  dense canonical scalars.

All 16 `AvxMatrixNoKPrimitives.*` tests pass under Clang 21.1.8 on Peach's
native AVX-512IFMA host. This includes the retained base-field, Fp2, G1, Fp12,
GT, sparse multiplication, and final-exponentiation cases, so the shared
addition change did not regress an earlier candidate slice. The focused
portable-AVX Debug G2 test translation unit also compiles locally.

The first G2 complete-add run exposed a second manifestation of `AUD-215`.
The generic Fp2 adapter attempted to expand negative m2-only cross terms
component-wise without first preserving their signed integer representative.
G1 and G2 now share one complete-add schedule, and `FP2Ring::prep_expand`
delegates each component to the base ring's bound-aware lift. The initially
failing self-add and distinct-add cases pass after this repair.

## Native code generation

Clang 21.1.8 compiled focused probes with `-O3 -DNDEBUG -march=native`.
The complete G2 add is fully inlined: its `0x3771`-byte body contains 744 IFMA
instructions and no call. No selected-vector conversion wrapper or explicit
lane round trip appears in the candidate boundary.

The G2 double audit wrapper is 74 bytes and calls the existing specialized
`PointDouble` template body; its return adapter uses a 896-byte logical frame.
The constant-time GLS wrapper is 95 bytes and calls the existing specialized
GLS executor; keeping that large executor out of line matches the incumbent
structure. A detailed secret-dependent control-flow/table-access inspection
and performance comparison remain selection gates.

The shared G1 complete-add helper still compiles as one inlined body with no
call. Its current instruction count is recorded for remeasurement rather than
silently treating earlier code-generation evidence as immutable after the
later fixed-R2 and shared-helper changes.

## Exact correspondence inputs

- `matrixnok_curve.h` SHA-256:
  `9fec0e6abbdf727011a7b4f859f0ae551235ea9aee9d9f58c49ff61e12d13556`
- `matrixnok_g1.h` SHA-256:
  `cb0651d44a8390a624d846f78e51fc24f94afa98b349b6a6599f54bc18c5dfe6`
- `matrixnok_g2.h` SHA-256:
  `f9f42440acbb729f80cd5fc690e9dd20de54c624617aaf3dad211fe9af4b3b75`
- `scalar_mult.h` SHA-256:
  `b1a3a30196a5e1a24bcaddf9963791aba9115fdfe5ecc1abcef1b4c9d8004396`
- `fp2.h` SHA-256:
  `a474cbbbb3bb3231698e7f33c653e73e9966adecceea0bd94eae44d70c834128`
- `matrixnok_fp.h` SHA-256:
  `36143f479eee5460595bf632b5a9e803388a6549c505d9939000e7c853bfeadf`

## Residual gates

- One approved serial same-process benchmark must compare G2 add, double, CT
  and VT GLS, and a representative end-to-end G2 workload with both the
  incumbent and pinned VROOM implementations.
- Optimized CT GLS still needs an explicit secret-dependent branch,
  table-address, spill, and stack-frame review at the final selected source.
- The GLS owner still parses a small fixed constant set during construction;
  production selection should decide whether generated compile-time constants
  are worthwhile based on measured construction and integration behavior.
- MatrixNoK inversion, Miller-loop, pairing integration, and formal
  arithmetic/bound proofs remain open.
- Broad suites, sanitizers, fuzzing, Apple CI, CUDA, BTx/BTE, and portable-AVX
  Release were outside this focused packet.
