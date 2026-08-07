# MatrixNoK G2 Candidate Evidence

## Packet boundary

- Candidate implementation: `2fc0ae98f962f5ab2259ed2d368f6328713b7310`
- Focused benchmark source: `b83f967328834c727ed381bd8b868b270083f0bf`
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
structure.

The optimized CT table scan contains no gather instruction or scalar-indexed
table selection. Its conditional branches compare fixed public loop counters
against 4, 17, 64, or 129; Booth selection remains a full fixed scan using
conditional moves. The multiply owner has one fixed 16-entry table-building
loop and fixed calls to the constant-time split, endomorphism, executor, copy,
and wipe helpers. No scalar-dependent branch or address was found.

The candidate CT multiply owner has a `0x18100` (98,560-byte) frame and calls
an executor with a `0x1700` (5,888-byte) frame, giving peak local stack use
above 102 KiB. The incumbent has the same owner frame and pinned VROOM uses
`0x183c0` (99,264 bytes), so this is not a candidate regression. It is recorded
as `AUD-218` because such a large inherited frame remains a production risk.
The candidate executor is also large (`0x40b2` bytes versus VROOM's `0x071e`)
because the complete MatrixNoK add schedule is expanded into the fixed table
scan. Native timing is favorable, but mixed-workload instruction-cache effects
remain worth watching.

The shared G1 complete-add helper still compiles as one inlined body with no
call. Its current instruction count is recorded for remeasurement rather than
silently treating earlier code-generation evidence as immutable after the
later fixed-R2 and shared-helper changes.

## Native performance gate

One Google Benchmark process ran at a time on Peach's AMD Ryzen 9 7950X. One
Clang 21.1.8 Release executable contained the incumbent, CatBLST candidate,
and actual pinned VROOM implementation. Every case used a valid subgroup
generator, a dependent operation chain, and the same dense canonical 256-bit
scalar. Seven repetitions produced the following means:

| Operation | Incumbent | CatBLST MatrixNoK | VROOM | vs incumbent | vs VROOM |
| --- | ---: | ---: | ---: | ---: | ---: |
| G2 add | 409 ns | 291 ns | 309 ns | -28.9% | -5.8% |
| G2 double | 258 ns | 219 ns | 217 ns | -15.1% | +0.9% |
| G2 CT GLS | 52.226 us | 40.969 us | 42.372 us | -21.6% | -3.3% |
| G2 VT GLS | 41.195 us | 31.360 us | not available | -23.9% | n/a |

Reported coefficients of variation were 0.01% through 0.30%. CPU frequency
scaling was enabled, so the figures are selection evidence rather than final
release claims. The candidate wins every incumbent comparison and remains
within 1% of or faster than VROOM for all common cases. No benchmark-only
arithmetic alternative was retained.

The ordinary benchmark target attempted to rebuild unrelated `blst_static`
code and stopped at the pre-existing `AUD-210` AVX pairing-context assertion.
The focused G2 and pinned-VROOM translation units were therefore compiled
directly from the current headers and linked with the existing focused native
partial archive for unchanged scalar/assembly support. Candidate, incumbent,
and VROOM G2 kernels all reside in the newly compiled benchmark objects; no
stale G2 kernel object contributed to these measurements.

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

- `AUD-218` tracks the inherited greater-than-102-KiB peak CT GLS stack use;
  direct final-table construction should be evaluated separately and retained
  only with correctness and serial performance evidence.
- The GLS owner still parses a small fixed constant set during construction;
  production selection should decide whether generated compile-time constants
  are worthwhile based on measured construction and integration behavior.
- MatrixNoK inversion, Miller-loop, pairing integration, and formal
  arithmetic/bound proofs remain open.
- Broad suites, sanitizers, fuzzing, Apple CI, CUDA, BTx/BTE, and portable-AVX
  Release were outside this focused packet.
