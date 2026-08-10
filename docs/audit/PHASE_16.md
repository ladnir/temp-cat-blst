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

Correspondence follow-up `f91268d` removes the last assumed-residue boundary.
It proves exact reconstruction of the real split IFMA accumulators for R1,
shifted R2, and ConvertTo; validates ConvertTo's selected `800/200` reducer
with a tighter canonical-input bound below `16m²`; derives both actual import
bases; and specializes the full multiplication/export theorem to those exact
producers. The source review found no additional C++ defect. The 3,024-job
Lean target and both deterministic constant checks pass, with no admitted
claims or project axioms. No C++ build, benchmark, CUDA, or broad suite ran.

## Selected Fp2 P3 proof boundary

The initial selected-Fp2 packet proves the source-shaped `FP2Ring` formulas,
prepared negation, variadic strided layout, pointwise batched R1 recurrence,
and the exact selected delayed-reducer bounds. A follow-up signed-domain review
constructed an abstract uncompleted expression `-M1`, whose correction-free
R1 representative is `-1`, and proved that R2 preserves that boundary value.
This does not describe the selected physical path: its negative static range
causes `FpRingCandidate::complete()` to add `max_neg_offset` before reduction.

Proof commit `8cbdc4c` models the semantic generator input for that completion,
proves the completed representative lies in `[0, 22p)`, and proves exact R2
correction throughout that range. Commits `fa2af08` and `aa206e8` prove the
generated wide offset lanes in M1 and M2, the raw real/imaginary multiply and
square producers and bounds, and the complete physical M1/R1-change-base/M2/R2
path. Selected-Fp2 correspondence row `AVX-P3-001` is accepted at the exact
source hashes. Fp12 and later P3 layers remain pending. The focused 3,013-job
Fp2 target, 3,026-job root target, and trust-word scan pass.

## General Fp12 P3 proof boundary

`AVX-P3-002` starts with the general field operations rather than mixing in
Frobenius constants or GT-only cyclotomic identities. The first Lean packet
fixes the selected C++ layout as six consecutive Fp2 coefficients split into
two Fp6 halves, proves the `v^3 = 1 + u` and `w^2 = v` tower interpretation,
and proves both the flattened multiplication and optimized square schedules.
The largest scalar coordinate has 22 signed products. Its simple semantic
ledger is 88 units. The literal pass corrects the physical accounting: direct
unsigned subtraction is at most 40 units, a single-offset coordinate is below
`169m^2`, and square coordinates that double the offset-bearing expression are
below `242m^2`. The actual general-field call uses MatrixNoK `MAX_ADD=800`;
the local `ready<200>` helper is cyclotomic-only. A compile-only forced
`ready<200>` probe correctly rejected `Bounds<128,210>`, confirming the type
system and clearing the suspected C++ defect. Lean checks all sixteen offset
lane encodings against `9 * 1440 * p^2` and proves their physical bounds.

The physical reduction chain is now composed through selected M1 reduction,
R1 change-base, M2 reduction, a completed representative in `[0,18p)`, and
exact R2 correction. `Fp12Raw.lean` now transcribes the twelve literal scalar
polynomials, proves both multiply and square reconstruct the six Fp2 outputs,
and checks their operation counts and exact offset multiplicities. The literal
ledger proves every unsigned subtraction is dominated by its offset and bounds
the multiply and square RNS representatives by `4608p^2` and `5760p^2`.
The generalized completed-representative theorem composes those values through
the full selected reduction chain, so `AVX-P3-002` is accepted at the recorded
source hashes. The remaining reviewed boundaries are the already accepted
Fp2 prepared-negation semantics and common P0 compiler/IFMA semantics.
This packet changes no production code and runs no
benchmark, CUDA, product build, or broad test suite.

## Fp12 Frobenius P3 proof boundary

`AVX-P3-003` proves the three selected Frobenius schedules for general Fp12.
The proof gives the exact concrete tower a commutative-ring structure, proves
the BLS12-381 modulus prime from a Lean-checked Lucas factor tree, and derives
the literal map-1 schedule as the `p`-power map. It factors the literal map-2
and map-3 code into coefficient-parametric linear templates, checks their
twelve basis images against repeated map 1, and concludes universal `p^2` and
`p^3` exponentiation theorems. This proof split avoids oversized kernel terms;
it does not alter or simplify the production schedule.

The physical packet binds all ten canonical coefficients and all 160 generated
M1/M2 lanes to the exact `to_mont_exact` construction, their MatrixNoK
representatives, standard P2 encodings, and canonical field values. Frobenius
forms at most two products per reduced coordinate: standard inputs below `2m`
and constant lanes below `m` give a tight `<4m^2` source bound versus the
literal `ready<MAX_ADD=800>` boundary. The proof also models
`standard_negate`'s `2m-value`, generated `40p` standard-zero offset, and
`reduce_auto<2>` conditional subtraction, then composes the resulting lanes
through the already accepted completed-reduction theorem. No new C++ defect or
contract weakening was found. Source owners and hashes are recorded in
`CORRESPONDENCE.csv`; the common compiler/IFMA boundary remains P0. Focused
Lean targets and deterministic constant regeneration are the validation gate;
no product build, benchmark, CUDA, BTx/BTE, sanitizer, Apple CI, or broad suite
is part of this packet.

The same focused investigation found `AUD-223`: internal `FP2Ring::mul` and
`square` exposed the one-element batch returned by `batch_reduce_expand`
instead of returning its element. Commit `a01a39b` selects `[0]`, adding no
runtime work, and the focused portable-Debug completion-boundary regression
matches scalar Fp2 for multiplication and squaring. No benchmark, CUDA, broad
suite, or portable-AVX Release build ran.

## Fp12 cyclotomic-square P3 proof boundary

Proof commit `256b980` accepts `AVX-P3-004` for the policy-selected split
cyclotomic square. The Lean model retains the exact first-six/second-six source
polynomials, factors three and six, concatenation order, scalar Frobenius
predicate, and general-square schedule. It proves a source-independent residual
identity for all twelve lanes and specializes it only with the exact scalar
`is_cyclotomic` relation. Thus zero is intentionally accepted, general Fp12 is
not, and no GT membership assumption is smuggled into the theorem.

Both generated constants are reconstructed from their canonical values and
`to_mont_exact`; all 32 M1/M2 lanes satisfy the standard P2 contracts. The
source-count physical ledger proves the first split is bounded by 28 M1 units
and 888 RNS units. The offset-bearing second split peaks at 101 M1 units and
2041 RNS units, with its sole subtraction dominated by the generated 80/1440
lower offset; both are below the selected 800/5760 reduction boundaries. The
post-reduction factors three and six fit the declared common bound six.

The focused `CatBLST` Lean target, generated-constant reproducibility check,
trust-word scan, and `git diff --check` pass. Existing portable and Peach native
scalar-differential evidence remains bound by the recorded unchanged source
hashes. No product build, benchmark, CUDA, BTx/BTE, sanitizer, Apple CI, broad
suite, or portable-AVX Release build ran. GT exponentiation and the relevant
final-exponentiation formulas remained P3 at that checkpoint.

## Fp12 final-exponentiation P3 proof boundary

Proof commit `21887f9` accepts `AVX-P3-005`. The new concrete theorem proves
the selected flat conjugation is Frobenius six for every Fp12 value, preserving
the distinction between its general-field meaning and inversion on a nonzero
cyclotomic value. The group-level source transcription retains the exact
`mul_n_sqr` loop, constants 2/3/9/32/15, negative BLS parameter convention,
easy part, and all named locals and Frobenius placements in
`final_exp_after_inverse`. The checked parameter identities show that the
complete schedule implements CatBLST/BLST's intentional cubed production
pairing exponent `3 * ((p^12 - 1) / r)`.

The physical ledger counts 35 Fp12 multiplications, 315 cyclotomic squares,
nine conjugations, one Frobenius-1, two Frobenius-2, one Frobenius-3, and one
inversion. These are serial calls into accepted P3 operations; the final-
exponentiation header introduces no fused accumulator or wider intermediate
bound. Its `mul_n_sqr` helper explicitly recasts to `CommonElement<6>`, while
the other multiplication calls return accepted normalized outputs. The
precise residual boundary is the documented input
contract: the Miller result must be nonzero and the supplied inversion provider
must return its mathematical inverse. The concrete inverter remains supported
by focused scalar-differential execution rather than a new formal inversion
proof.

The 3,067-job focused root Lean target, constant self-test, generated-Lean
reproducibility check, trust-word scan, and `git diff --check` pass. No product
build, benchmark, CUDA, BTx/BTE, sanitizer, Apple CI, broad suite, or portable-
AVX Release build ran. GT binary and GLV4 exponentiation were the remaining P3
packet at that checkpoint.

## GT exponentiation P3 proof boundary

Proof commit `82f26c1` accepts `AVX-P3-006` for public GT multiplication. Lean
models the production CT/public-VT scalar preparation as reduction modulo `r`,
the exact divide-by-`z^2` then divide-by-`z` split, all digit bounds, and exact
recomposition. It proves `p = -z (mod r)`, so the source's conjugated odd
Frobenius basis is precisely the unsigned radix basis.

The functional packet constructs the literal sixteen table entries and the
descending 64-round square/select/multiply loop. Its invariant proves the CT
result, while a separate literal `started` state machine proves leading-zero
and zero-window skipping leaves the public-scalar VT result unchanged. The
binary CT select loop, fixed 256-bit basic VT loop, preceding `nbits` mask, and
public 129-through-256 dispatch are included. All claims are on GT; general
Fp12 and the broader cyclotomic relation are not admitted as inputs.

The source ledger records 16 table entries, 11 table-build multiplications,
three Frobenius maps, two conjugations, 64 fixed loop iterations, 64 squares,
64 loop multiplications, one normalization multiplication, 960 candidate
comparisons, 11,520 coefficient selections, and `CommonElement<6>`. The VT
path has at most 63 loop squares and 63 loop multiplications after its first
nonzero joint bit. Focused root Lean, constant regeneration, trust-word scan,
and `git diff --check` pass. No product build, benchmark, CUDA, BTx/BTE,
sanitizer, Apple CI, broad suite, or portable-AVX Release build ran.

The remaining P3 item is the separate balanced signed bytecode compiler and
executor used by FFT: prove four signed digits plus sign mask, the 70-byte
schedule bound, chunked square runs, executor semantics, and its FFT call-site
preconditions. Existing differential coverage is evidence, not a formal
substitute for that row.

Proof commits `51462b9` and `effdb3d` close `AVX-P3-007`. The arithmetic packet
proves the four centered divisions, `z^4` carry fold, tail at most two,
Frobenius-basis reconstruction, and strict signed 64-bit digits. The compiled
packet then proves sign-mask bits, absolute magnitudes, high-to-low discovery,
literal opcode nibbles, fifteen-square chunks, trailing runs, and executor
semantics. For every 256-bit scalar the compiled executor returns the original
GT exponent. The physical proof bounds the schedule by 63 of 70 bytes, squares
and post-initialization products by 62 each, and total products by 74 including
table construction and normalization. A canonical-`Fr` corollary covers FFT
prepared roots and inverse scales; the C++ `from_mont` representation bridge
and P0 lowering remain explicit residual boundaries.

All P3 rows are accepted at exact source and generated-constant hashes. The
next formal packet is P4 curve/Miller-loop algebra. No product build,
benchmark, sanitizer, CUDA, BTx/BTE, Apple CI, broad suite, or portable-AVX
Release build ran for this proof-only closure.

## Miller doubling P4 proof boundary

Proof commit `8e0dc2d` begins `AVX-P4-001` against exact `miller.h` and `GT.h`
hashes. A new Fp2 algebra bridge packages the accepted P3-001 formulas without
changing their representation. The literal shared doubling schedule then
proves its homogeneous output remains on `y^2 = x^3 + b`, its line contains
the input and negative doubled output, and prepared G1 line scaling equals the
direct Miller construction. No C++ discrepancy was found.

Follow-up proof commit `adc9ffc` proves literal `flat_mulxy00z0` sparse
multiplication, its 72-product ledger and 2880-unit RNS bound, complete mixed
addition geometry, the exact 68-operation BLS parameter schedule, curve
preservation across that schedule, and final Frobenius-six conjugation. The
focused 3,079-target root Lean build passed. No C++ discrepancy was found.

Proof commits `e992ce1` and `4dbb3cd` close the physical and end-to-end gaps.
Every accumulated Fp2 coordinate fits the proved four-product `16m^2`
envelope; generated completion, M1/R1/M2/R2 reduction, batch permutation, and
both exact loop ledgers compose at that bound. Literal flat square and sparse
multiplication—including the omitted initial identity square—compose with the
binary BLS-x loop and its curve invariant. The 3,082-target focused Lean build,
clean trust-word scan, and three portable-Debug scalar differentials pass.
`AVX-P4-001` is accepted at the recorded hashes. P0 compiler/ISA lowering is
explicitly external, as for accepted P2/P3 rows. No broad product or other
permission-gated run was performed.

## MatrixNoK production integration

Production AVX point, GT, pairing, G1/GT FFT, and middle-product dispatch now
selects the exact MatrixNoK arithmetic owners through zero-overhead semantic
aliases.
The incumbent `BoundedRing` engine remains only as an in-tree historical
comparison implementation and is not reachable from those production owners.
The exact graph, integration repairs, completed checks, and remaining native
release gates are recorded in `docs/avx/MATRIXNOK_INTEGRATION.md`.

Integration found `AUD-225`: raw RNS equality is not semantic equality under
MatrixNoK's redundant representatives. The Peach native run then found
`AUD-228`: the same mistake affected public infinity observation and secret
normalization, including signing the identity. Canonical export is now limited
to those semantic boundaries; the secret path reuses the conversion already
required for inversion, and every optimized arithmetic-probe frame is
unchanged from the preceding commit.

`AUD-226` closes the proof-to-production binding gap. A structural validator
pins the selected aliases, the generated Lean snapshot plus
`GtGeneratorConstant.lean` binds the new GT literal's 192 lanes, and accepted
`AVX-P4-002` records final exact hashes. Peach Clang 21.1.8 builds the native
Release production library at `6ce0d44`, passes all 34 focused tests and all
three MatrixNoK audit CTests, and emits no out-of-line `_vv` or conversion
symbol in the optimized probe. No benchmark or broad run was performed.

Native object inspection opened `AUD-229`: several AVX owners and public
literals still use function-local dynamic initialization guards. Removing them
without perturbing hot code or theorem bindings is the next integration cleanup
packet.

The `AUD-229` implementation at `fe09e26` replaces those derived first-use
objects with constant-initialized owners and separately generated curve
literals. The curve artifact is independently regenerated from canonical
BLS12-381 values and the existing MatrixNoK parameters without changing the
already proved constants header. All 35 focused portable-Debug tests pass, and
MSVC production C API objects contain no thread-safe-static or dynamic-init
symbols. Because `point.h`, `GT.h`, `gt_context.h`, and the Fp2 value
constructor are exact-hash proof inputs, the five affected correspondence rows
were deliberately held `in_review` pending one bounded native Clang reacceptance;
the evidence below reaccepts them.
No benchmark or broad build ran.

The first native reacceptance run found `AUD-230` in the audit infrastructure,
not the arithmetic: exact correspondence hashes had been captured from a
Windows CRLF worktree and failed against Peach's LF checkout. Correspondence
hashing now canonicalizes ordinary CRLF to Git's LF text and rejects bare CR;
all CRLF-derived ledger hashes have been replaced by their committed LF blob
hashes. Local validation covers all 19 rows and 51 mapped hashes. Native
correctness had already passed 35/35 when this gate stopped the run. The
corrected audit passes on that clean Peach checkout at `63e2feb`.

The bounded Peach reacceptance validates `AUD-229` and `AUD-230`. Clang 21.1.8
Release passes 35/35 focused native tests and all five audit gates. The old AVX
C API object had 231 guard/atexit relocations and 11 relevant guarded owners;
the new object has zero of each. Its text shrinks from 701,649 to 682,541 bytes
and its BSS from 5,393 bytes to zero. Public generator and Fp12-one accessors
reduce to direct eight-byte loads. The selected arithmetic code-generation
probe is byte-identical before and after, its stack-usage rows are unchanged,
and GLV/GLS owner frames remain `0x380`/`0x680`; Miller-loop frames are
unchanged or smaller. Therefore no benchmark ran. The five source-bound proof
rows are reaccepted at their canonical hashes.

### AUD-227 public curve proof packet (accepted at `f1cfde9`)

`AVX-P4-003` closes the selected public G1/G2 arithmetic packet without
broadening the accepted claim. `CurveComplete.lean` transcribes the shared
MatrixNoK complete-add schedule and public Algorithm-9 doubling schedule. An offline
polynomial certificate, rechecked by Lean, proves complete addition preserves
the short-Weierstrass curve over a commutative ring; identity, inverse,
commutativity, add-self/double agreement, and selected Fp/Fp2 instances are
also proved. `CurveScalarSplit.lean` proves the exact G1 radix-`z²` and G2
four-way radix-`z` decompositions. `CurveBoothExecutor.lean` proves the selected
26-window and 13-window constant-time executors and their exact capacity bounds.
`CurveWnafExecutor.lean` proves the source add/subtract-and-shift recurrence,
odd-table schedule, executor recomposition, and exact 129/65-slot exhaustion
using backward minimum-predecessor bounds.

`CurveEndomorphism.lean` resolves the signed-parameter ambiguity directly
against the source constants. The AVX G1 beta is proved to be the square of
scalar BLST's sigma beta, making `(beta*X,-Y,Z)` the negation of `sigma²`;
BLST's `sigma²(P)=[-|z|²]P` relation therefore yields the unsigned
`[|z|²]P` split basis. The G2 component schedule is proved exactly equal to
`-psi`, so BLST's signed `psi(P)=[-|z|]P` relation yields `[|z|]P`. Both maps
are proved curve-preserving over the selected Fp/Fp2 fields, including all
literal coefficient identities.

The first subpacket left the row `in_review` pending concrete subgroup and
projective-group composition. No C++ source, benchmark, native build, or Peach
run was part of that proof-only step.

The next packet adds `CurveSubgroup.lean`. It proves `sigma³ = 1` and
`psi⁶ = -1` from the literal coordinate constants, checks
`z⁶ + 1 = r(z² + 1)`, checks both published cofactors and their required gcds,
and derives relation-implies-`[r]P=0` for G1 and G2 from an explicit full curve
cardinality premise. The converse is proved from an explicit prime-subgroup
generator-span premise and generator relation, yielding conditional iff
theorems without smuggling a point-count or cyclicity claim into the model.
The G1 coordinate map is also proved to commute exactly with the selected
complete-add schedule. The G2 semilinear bridge is factored into a generic
ring-map theorem and a diagonal-scaling theorem; complete addition after `psi`
is one common homogeneous scale of `psi` after complete addition, and the
literal scale is proved nonzero.

`CurveGroupBridge.lean` then starts the concrete projective-group composition.
It instantiates selected G1 as `y²z = x³ + 4z³` in Mathlib over the already
proved-prime `ZMod p`, proves the discriminant nonzero, and maps every nonzero
on-curve homogeneous representative to a nonsingular projective point. On that
domain CatBLST's representative equivalence is proved iff Mathlib projective
point equality, including the point at infinity, and the map is proved to
preserve zero and negation. The all-zero coordinate triple is excluded
explicitly rather than silently admitted by the homogeneous equation. The
next proof establishes the missing selected-Fp2 `Field` instance. It checks
`p mod 4 = 3`, uses Mathlib's quadratic-residue theorem to prove that the norm
`a² + b²` cannot vanish for a nonzero pair, and verifies the explicit
conjugate-over-norm inverse. The generic projective carrier is then instantiated
for the selected G2 twist. Its discriminant is checked nonzero, and the same
nonzero-representative conversion, exact quotient equality, zero, and negation
theorems are proved. `CurveGroupAdd.lean` then proves the complete formula's
bihomogeneous scaling, reduces arbitrary representatives to affine or infinity
normal forms, and identifies every valid selected G1/G2 result with Mathlib
group addition. The proof separately covers infinity, doubling, inverse, and
generic affine inputs; it also proves the selected curves have no affine
2-torsion using finite-field cube criteria, so the RCB exceptional divisor
cannot occur.

Proof commit `f1cfde9` completes the selected multiplication composition.
`CurveBinaryMul.lean` defines a logarithmic proof-only evaluator from the exact
complete-add owner and proves it equal to Mathlib projective scalar
multiplication. `CurveGenerators.lean` checks the standard G1/G2 coordinates,
proves both generators are killed by `r`, and proves the literal AVX beta and
minus-psi maps produce the required radix-`z²` and radix-`z` bases on those
generators. `CurveDispatch.lean` models public bit masking and canonical
modulo-`r` preparation, composes the exact selected Booth and wNAF exponents,
and proves G1/G2 `mult` and `mult_vartime` over both `r`-torsion carriers and
the conventional cyclic subgroups generated by those concrete points.

The cyclic subgroup is the G1/G2 API domain, so a full ambient-curve point
count is not a premise of the multiplication theorem. The separate assertion
that scalar BLST's relation-based `in_g1`/`in_g2` predicates characterize
exactly those cyclic subgroups is now bounded directly in
`CurveMembership.lean`. Proof commit `24de7b4` proves the literal relation maps
are quotient-safe and additive, proves the infinity/x=0 small kernels, and
establishes `standard cyclic subgroup ⊆ accepted relation ⊆ r-torsion` for
both curves without the old full-cardinality premise. Commit `1fdc201` further
turns those exact accepted relations into `ZMod r` subgroup carriers, proves
the selected unsigned GLV/GLS bases on them, and instantiates all four selected
CT/VT multiplication theorems. Commit `65b54f6` proves `r` prime, proves the
standard generators have exact order `r`, and proves span equality is
equivalent to `card accepted ≤ r`. Commit `f6c251b` proves the G1 upper bound
and exact span equality using the generic `2p+1` affine-point bound and the
checked inequality `r² > 2p+1`. Commit `475221a` narrows G2 further:
`r³ > 2p²+1` forces the accepted kernel to have cardinality exactly `r` or
`r²`. El Housni, Guillevic, and Piellard, ePrint 2022/352, Proposition 5,
exclude the rank-two case for the BLS12 family. Commit `a8075e7` certifies the
exact BLS12-381 seed, cofactor polynomial, and coprimality hypothesis in Lean
and cites the corrected proof at the source owner. `AUD-231` is validated;
the published curve-family theorem and twist order remain explicit external
mathematics rather than a claimed foundational Lean point count. `AVX-P4-003` is accepted
for selected arithmetic, while compiler/ISA lowering remains the same
external P0 boundary. The new packet passes focused module and transitive root
type-checks, a trust-word scan, and an axiom inventory. No C++
build or benchmark ran for this proof-only completion.

### GMP/MPIR validation cleanup (`AUD-232`)

Commit `9573d58` removes the optional GMP/MPIR build dependency and its old
runtime-constant generator, `BigInt` overloads, validation/reference tests,
benchmark harnesses, and remote comparison scripts. The selected MatrixNoK
engine remains fixed-width and unchanged: its constants continue to be
regenerated independently with Python, bound into Lean, compared against the
scalar backend, and checked by the exact correspondence ledger.

The GMP-free MSVC portable-Debug configuration rebuilt the production library
and the single test executable. All 30 focused `AvxMatrixNoKPrimitives` and
`AvxCApiContract` tests passed. Structural search finds no live GMP/MPIR,
`BigInt`, or removed CMake-option references outside historical documentation.
The audit manifest accounts for all 784 retained files. Removing disabled
branches changed three full-file hashes, so six accepted correspondence rows
were rebound to `9573d58`; all 20 rows and 59 mapped hashes validate. No
benchmark, portable Release build, native run, CUDA work, or BTx audit ran.
