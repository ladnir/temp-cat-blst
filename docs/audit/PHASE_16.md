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

The inspection also recorded `AUD-218`: CT GLS initially peaked above 102 KiB
of local stack due to inherited source/final table duplication. Direct
construction in the final four-table owner reduces the multiply frame from
98,560 to 50,176 bytes and peak local stack to about 55 KiB. All 16 focused
native tests pass, the fixed executor is unchanged, and one serial run improves
CT GLS from 40.969 to 39.911 us. `AUD-218` is validated. Exact source hashes,
timing, code-generation facts, and residual gates are in
`docs/avx/G2_EVIDENCE.md`.

## MatrixNoK Miller and pairing preview

Implementation commit `eaa7822431954980694341dc77cdcd8b12cce7eb`
parameterizes the existing Miller, pairing, and inversion templates over the
Fp12 policy while preserving their incumbent defaults. The new MatrixNoK
aliases are isolated and do not alter public symbols or production dispatch.

Local portable-AVX Debug compilation and Peach Clang 21 native execution pass
all 17 focused MatrixNoK tests. The new pairing test compares four independent
G1/G2 combinations with scalar CatBLST at both the post-Miller/final-exponent
boundary and the complete pairing result.

One serial same-process benchmark measured MatrixNoK Miller at 55.066 us,
20.7% faster than the incumbent and 5.3% faster than pinned VROOM. Complete
pairing measured 141.775 us, 14.8% faster than the incumbent and 1.5% slower
than VROOM. The preview passes the agreed 2% selection tolerance. Optimized
inspection finds bounded peak stack of about 8 KiB and no allocation,
exception-runtime, gather/scatter, or selected-vector conversion in the
candidate Miller body.

Exact evidence is in `docs/avx/PAIRING_PREVIEW_EVIDENCE.md`. MatrixNoK remains
non-dispatched.

Follow-up commit `a1c8f9f652b5ae7d8e77966c454f312e3079bd95` compared the
inversion and hard final-exponent paths with pinned VROOM. The exponentiation
algorithm is already the same; optimized inspection instead found excessive
spill pressure in the monolithic MatrixNoK cyclotomic square. Two fixed
noinline six-output halves preserve all formulas and IFMA work while reducing
the candidate frame from 4,352 bytes to 1,408 and 512 bytes. The incumbent
continues to compile through its original monolithic path and retains its exact
representative symbol size.

All 17 focused native tests pass on the final source and the portable-AVX Debug
translation unit compiles. One serial seven-repetition process measured full
pairing at 132.084 us, 21.2% faster than the incumbent and 5.5% faster than
pinned VROOM. Performance selection is therefore complete. Formal proof,
production dispatch review, and release evidence remain pending.

## MatrixNoK proof preparation and semi-freeze

Commits `ddb26f7`, `5890075`, and `adf2ed5` prepare the performance-selected,
still non-dispatched candidate for formal correspondence. The new
`docs/avx/MATRIXNOK_CONSTRUCTION.md` records the selected call graph, two-base
representation, import/export and multiplication pipelines, field-tower
coordinate order, fixed-matrix accumulator semantics, exact bound ledger, and
the source duplication retained for latency/throughput code generation.

Source constants now name and assert the exact direct-multiply bounds 4, 16,
and 20 and generated base redundant multiple 40. The generic delayed-reduction
path deliberately retains its conservative C++ type allowance 32: narrowing it
to 16 changed higher-tower optimized assembly and was rejected without a new
benchmark gate. Lean should prove the stronger mathematical contribution 16;
the weaker C++ type claim is not the theorem target.

The pass found `AUD-219`: the focused compatibility target compiled its scalar
differential oracle but did not link the scalar definitions. It now compiles
the four required scalar owners directly, independent of the configured
production-default backend. MSVC portable Debug passes the four affected RNS
and base-field cases; a clean Peach Clang 21 native-plus-portable build passes
all 17 MatrixNoK tests. Optimized code-generation probes are byte-identical to
the performance-selected `e6b2995` baseline. No benchmark or broad suite ran.

The correspondence ledger is rebound to exact commit `adf2ed5`, current source
and generated-constant hashes, and a new fixed-Montgomery-reduction obligation.
Mapped arithmetic is now semi-frozen: executable changes reopen the associated
record and require focused correctness/code-generation evidence, plus a
benchmark only when generated code materially changes.

The first Lean P0 slice now builds with the pinned Lean 4.32.1/Mathlib 4.32.1
workspace. It defines the 64-bit/eight-lane model and abstract IFMA operations,
proves the 52-bit product split and generic no-wrap accumulation rule, and
separately transcribes the portable C++ 26-bit-half schedule. Machine-checked
theorems prove that both portable low/high results equal the abstract IFMA
digits for all 52-bit inputs. Correspondence remains `in_review` because the
C++ transcription, compiler/unsigned-word interpretation, and native x86 IFMA
semantics are explicit residual assumptions rather than Lean-proved facts.

P1 now has a deterministic full-literal binding from the generated C++ header
to Lean, guarded by an exact source hash and a no-executable CTest check. The
first constant theorem packet proves basis dimensions, odd/distinct/pairwise-
coprime moduli, exact products, dynamic range, required inverse and multiplier
relations, table shapes and ranges, and the cyclic R1 permutation map. The
second packet defines the CRT generator in Lean, validates its inverse
witnesses and identity residue patterns, regenerates all four fixed table
families, bounds each ceiling estimate, and proves the 128-bit quotient
accumulator envelope. P1 is still open for the complete kernel-specific bound
ledger and field-tower encoding relations; the checked-in data snapshot and
Python validator are not counted as substitutes for these Lean claims.

## MatrixNoK canonical-import repair

Formal P2 work found `AUD-221`: VROOM's ConvertTo ceiling estimates can choose
a negative CRT representative for a valid canonical BLS12-381 input, while the
immediately following R2 and ConvertFrom schedules require a nonnegative
representative. The exact witness fails both import/export and
multiply-by-one under the inherited constants.

Commit `5f0399e` changes only ConvertTo's eight fixed-point quotient constants,
each by minus one, and updates the independent generator to use floor for that
table. R1, R2, ConvertFrom, and all runtime schedules remain unchanged. Lean
independently regenerates the repaired constants, proves the estimate exact or
one low, proves the representative lies in `[0,2p)` and preserves its field
value, composes it with the R2 contract, and proves the exact former witness
now selects the positive complement. The proof target builds without trust
escapes.

The enabled exact portable-Debug regression passes. Follow-up commit `1e286e2`
fixes `AUD-222`: scalar almost-inverse output has a private `<2rR` contract,
while public `Fr::redc` correctly retains its narrower `<rR` contract. The
private inverse path now checks its exact bound in Debug and performs the same
backend REDC-plus-`to_mont` sequence as Release, adding no Release work. All 18
focused MatrixNoK tests now pass locally.

Peach validation at exact commit `7a5f094` uses Clang 21.1.8 with `-O3` and
`-march=native`, native AVX-512 plus the portable compatibility layer, and no
CUDA, benchmark, or sanitizer target. All 18 focused MatrixNoK tests pass in
one process. A minimal canonical-import probe compiled at pre-repair commit
`73cf354` and at `7a5f094` has the same 209-instruction disassembly after
normalizing immediates (normalized SHA-256
`fc6f14e0917b310a0188aeccda1f03b2c0f65ae9dec79d5f81ca1cc6beba03c9`).
The raw diff contains exactly eight changed `movabs` immediates, each reduced
by one. No benchmark ran because neither repair changes the runtime instruction
schedule. `AUD-221` is validated.

Follow-up proof commit `e01ee59` closes the remaining reachable base-field P2
boundary. `CanonicalExport` models the exact eight-coefficient wide
radix-`2^52` result produced by non-shifted ConvertFrom and recombined by
`from_vec`; proves the quotient numerator stays below `2^128`; proves both
high and low IFMA accumulators stay below `2^64`; proves canonical
import/export; and composes repaired import, selected MatrixNoK multiplication,
and export to the canonical field product. The complete 3,024-job Lean target
passes without trust escapes. This proof-only commit changes no production
code and therefore requires no new benchmark or native test run.
