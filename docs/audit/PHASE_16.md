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

## MatrixNoK FFT preview

A benchmark-only same-process adapter now exercises the pinned VROOM Fp12
stack under CatBLST's FFT roots and GLV4 schedules. Its 256-element GT
forward/inverse round trip passes. On Peach with Clang 21, inverse GT FFT is
1.19x to 1.20x faster than the incumbent across 256 through 16,384 elements.
This is sufficient to prioritize the Fp12/GT port but does not select the
candidate for production. The adapter is fixed at `71bfbfd`.

The apparent G1 blocker was a benchmark-adapter defect: a mistyped decimal
generator coordinate created an off-curve input, while a same-implementation
reference check consumed the same bad value and masked the cause. Commit
`edaf5462` removes the copied coordinate, imports CatBLST's canonical
generator, and adds an independent homogeneous curve-equation gate. The
two-element butterfly and 256-element forward/inverse round trip now pass.

On Peach with Clang 21, forward G1 FFT is 1.10x to 1.12x faster than the
incumbent across 256 through 16,384 elements. The G1 packet should now adapt
MatrixNoK arithmetic to CatBLST's audited point owners, retain independent
scalar or mathematical oracles, and use this result as its performance target.
Exact evidence is in `docs/avx/FFT_PREVIEW_EVIDENCE.md`.

## MatrixNoK G1 candidate

Implementation commit `7e731d9` replaces the benchmark-only G1 preview with an
isolated production-shaped candidate. It reuses the audited complete
homogeneous formulas and GLV executors through a typed curve policy; public
symbols and production dispatch remain unchanged.

The initial delayed-expansion translation exposed `AUD-215`: a negative
small-m2 cross term must be lifted by a sufficient multiple of the field
modulus before base expansion or the RNS multiple is lost. The fixed helper
derives the minimum power-of-two offset from the static lower bound. Generated
offset constants are independently regenerated and validated.

On Peach, all 11 focused native AVX-512IFMA primitive tests pass. The G1 oracle
covers complete add/double and exceptional operands, deterministic chains,
both CT and VT GLV on two bases with boundary and dense canonical scalars, and
a forward/inverse FFT round trip against scalar CatBLST. MSVC also compiles the
native code-generation probe, and Clang's selected add contains no call or
vector conversion wrapper.

One serial five-repetition process measured the candidate G1 FFT 9.4% to 9.9%
faster than the incumbent. It ranges from 1.6% faster to 1.7% slower than the
pinned VROOM implementation, satisfying the 2% primary-case tolerance. Direct
add/double remain 12.6% and 11.3% slower than VROOM and are recorded as future
optimization targets, not concealed by the end-to-end result. Exact evidence,
source hashes, rejected alternatives, and residual proof gates are in
`docs/avx/G1_EVIDENCE.md`.

## MatrixNoK Fp12 and GT candidate

Implementation commit `1011916baf69a4ecb4b0e012c30a74051c8ec850`
adds an isolated production-shaped Fp12/GT candidate. The generic Fp12 owner
now receives a parameter-specific constants policy, while its default policy
preserves the incumbent source surface and production behavior.

The first translation exposed `AUD-216`: fixed wide-arithmetic constants were
selected by the shared 50-bit limb width even though the incumbent and
MatrixNoK rings have different RNS moduli. Explicit parameter binding and an
independent constant validator close that silent cross-basis correctness bug.
Code-generation review also closed `AUD-217` by retaining a fixed outer loop
instead of forcibly cloning all eight R2 rows, restoring VROOM-equivalent IFMA
counts and code size.

All 13 focused native AVX-512IFMA tests pass, including general Fp12,
Frobenius, cyclotomic, and CT/VT GT exponent cases against scalar CatBLST. One
serial Clang 21 process measured Fp12 multiplication at 303 ns, versus 362 ns
for the incumbent and 304 ns for pinned VROOM. The 256-element inverse GT FFT
is 13.9% faster than the incumbent and 1.0% slower than VROOM; the
1,024-element case is 13.6% faster than the incumbent and 3.1% slower than
VROOM. The latter and the 4.5% cyclotomic-square gap remain explicit residual
optimization work. Exact evidence is in `docs/avx/FP12_EVIDENCE.md`.

Follow-up commit `3e1f99a` keeps the existing sparse Fp12 kernel and hard
final-exponentiation schedule, but makes their compile-time field policy
explicit. Native differential coverage now includes sparse line
multiplication, exact aliasing, and three complete final-exponentiation chains;
all 14 focused cases pass. The incumbent AVX GT translation unit recompiles
unchanged. MatrixNoK inversion, G2, Miller, and pairing integration remain the
next correctness gates.

## MatrixNoK G2 candidate

Implementation commit `2fc0ae9` adds an isolated production-shaped G2 curve
policy and specializes the existing GLS executor through a compile-time curve
parameter. Public symbols and production dispatch are unchanged.

The initial G2 add exposed that `AUD-215` also applied through Fp2: negative
m2-only cross terms were expanded without preserving their signed integer
representative. G1 and G2 now share one complete-add schedule and Fp2 delegates
the bound-aware lift to its base ring component-wise.

All 16 focused native MatrixNoK tests pass on Peach under Clang 21, including
G2 complete exceptional cases and CT/VT GLS against scalar CatBLST on two
bases, boundary values, and deterministic dense scalars. The focused portable
Debug translation unit compiles locally. Optimized G2 add is fully inlined
with no calls; double and GLS retain the incumbent out-of-line template-body
structure.

One serial seven-repetition process measured candidate G2 add, double, CT GLS,
and VT GLS 28.9%, 15.1%, 21.6%, and 23.9% faster than the incumbent. Candidate
add and CT GLS are 5.8% and 3.3% faster than pinned VROOM; double is 0.9%
slower. Optimized inspection finds fixed control flow and full-scan table
selection with no gather or scalar-indexed address.

The inspection also records `AUD-218`: CT GLS peaks above 102 KiB of local
stack due to inherited source/final table duplication. The incumbent has the
same 98,560-byte owner frame and VROOM is slightly larger, so this is not a
candidate regression, but it should be redesigned and remeasured before final
release selection. Exact source hashes, timing, code-generation facts, and
residual gates are in `docs/avx/G2_EVIDENCE.md`.
