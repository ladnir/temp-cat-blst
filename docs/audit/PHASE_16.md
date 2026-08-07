# Phase 16: AVX Curves and Protocols

## Scope

Phase 16 performs the file-by-file hand audit of the native and portable AVX
curve, pairing, map-to-group, and multi-scalar implementation. Scalar code is
used as an independent differential oracle where appropriate. Native AVX-512
code generation and runtime evidence is required before performance-sensitive
files close; portable-AVX Debug supplies fast local correctness feedback.

BTx/BTE remains frozen, CUDA remains experimental, and neither is audited by
this phase. Portable-AVX Release, broad suites, sanitizers, fuzzing, remote
machines, and benchmarks remain explicit milestone or permission gates.

## AVX point review

The owning `blst/avx/curve/point.h` review follows the completed selected
facade review. It traces the AVX homogeneous representation, scalar boundary
conversions, fixed-window exact multiplication, GLV/GLS multiplication,
signing and key derivation, basic point operations, validation, encoding, and
visible helper functions through native C++, the selected C API, and focused
tests.

Initial findings are:

- `AUD-169`: signing and public-key derivation select AVX types but execute
  scalar point multiplication and redundant normalizations.
- `AUD-170`: equality and zero observation use scalar conversions and
  inversions; the dormant cross-product helper lacks infinity parity.
- `AUD-171`: `to_pk2` declarations lost the inherited pre-decorated-output
  contract and in AVX C++ state the wrong size and format.
- `AUD-172`: `times_minus_z` requires distinct input and output objects but did
  not state or assert that invariant.
- `AUD-173`: AVX selection removes the inherited native C++ `sign2` method and
  forces the C boundary to duplicate normalization and output logic.

The approved remediation keeps multiplication and coordinate arithmetic
AVX-native, retains a fixed-schedule scalar base-field inversion only as the
normalization boundary until an independently justified native inversion
exists, removes observation conversions, restores exact output documentation,
and adds zero-Release-overhead Debug invariants.

## Validation policy

The edit loop uses only the existing selected portable-AVX Debug target and
structural documentation checks. Focused tests cover G1/G2 scalar agreement,
canonical and all-zero infinity representatives, projective scaling, signing,
public-key derivation, optional outputs, pre-decorated finalization, exact
aliasing where promised, and rejected aliasing where prohibited.

Native AVX-512/IFMA runtime and optimized code-generation inspection remain
required before `AUD-169` and the file review are validated. No benchmark runs
without explicit approval.

## Focused remediation evidence

Implementation commit `5e4ab8dec6acbbfc117713f62a0f789676c9fd4c`
keeps GLV/GLS multiplication and homogeneous point arithmetic in AVX storage.
Only the final secret-derived Z coordinate crosses to the audited scalar
base-field inversion, which uses the fixed addition-chain `flt_reciprocal`;
the scalar input and inverse are wiped before return. `sign2` reuses that one
normalized product for affine and byte outputs, and the selected C boundary now
delegates to it.

Both `audit-phase1-avxportable` (scalar default) and
`audit-phase3a-avxportable` (AVX default) built the focused contract target
serially and passed 12/12 tests. The cases cover scalar/AVX G1 and G2 signing,
identity, exact input/output aliasing, public-key derivation, optional outputs,
the selected C route, pre-decorated output finalization, equality across scaled
homogeneous representations and both infinity forms, and distinct-object
`times_minus_z`. Rust formatting, public documentation, group-contract,
manifest, and `git diff --check` checks pass.

`AUD-170` through `AUD-173` are validated by this focused evidence.
`AUD-169` is fixed but remains open for validation until native AVX-512/IFMA
optimized code generation and runtime are inspected. No benchmark was run.

## MatrixNoK base-field candidate

The Phase 16 modernization project now has a complete, isolated base-field
multiplication slice through implementation commit
`cd282561dd6bccf69cf52bf35faed456fc6b6374`. It is not selected by production
dispatch. Local portable Debug and Peach native AVX-512IFMA differential tests
pass against explicit portable scalar Montgomery primitives, including field
boundaries, dense deterministic inputs, dependent chains, and exact aliasing.

Code-generation review caught an out-of-line wrapper, batch-of-one stack
adapters, variable-index rotations, and a cache-line-crossing generated
permutation table. The resulting Clang 21 arithmetic is fully inlined, uses
VROOM's immediate rotation schedule, and loads 64-byte-aligned table rows.

The same Clang 21 executable compiles both CatBLST and actual VROOM commit
`26bd9aee`. Same-process randomized measurements favor CatBLST by 0.816% over
six base-field cases and by 6.595% at batch 12. Batch 8 is borderline: the
final full run is +2.254%, while the longer focused confirmation is +1.799%.
The prior CatBLST 50-bit and 52-bit single-multiply paths measure 22.8 ns and
24.5 ns versus 19.1 ns for MatrixNoK.

`AUD-211` through `AUD-213` record the repaired code-generation,
cross-compiler, and table-layout defects. Mathematical correctness and bound
safety remain open until the selected field-tower design reaches the first
Lean milestone in `docs/avx/PROOF_PLAN.md`. Exact evidence and correspondence
hashes are in `docs/avx/BASE_FIELD_EVIDENCE.md`. No broad build was run.

## MatrixNoK Fp2 candidate

Implementation through commit `19f37a1` adds an isolated Fp2 candidate without
changing production dispatch. The four-product multiply and three-product
square combine terms into two wide outputs before reduction and base change.
Local portable Debug and Peach native AVX-512IFMA pass scalar differential,
boundary, exact-alias, output-bound, and dependent-chain tests.

The base and Fp2 owners are now empty, one-byte C++ objects backed by
compile-time-derived immutable fixed operations. Direct tests establish
correspondence with the former generic matrices and reducers. Native symbol
inspection finds no runtime constructor or guard, and Clang 21 emits fully
inlined fixed-schedule Fp2 bodies with no calls, conditional branches, or
stack frame.

In the final same-process randomized run, multiplication and squaring are both
33.6 ns: 18.8% and 18.2% faster than CatBLST's prior 50-bit path. Multiplication
is 3.2% faster and squaring 1.5% faster than pinned VROOM. Stateless base-field
multiplication is 19.4 ns, about 1.6% behind the former context-backed
candidate but 14.9% faster than the incumbent. `AUD-214` is validated. Exact
evidence and residual proof obligations are recorded in
`docs/avx/FP2_EVIDENCE.md`.
