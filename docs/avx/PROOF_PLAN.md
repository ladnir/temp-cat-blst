# AVX-512 Arithmetic Proof Plan

## Claim boundary

The proof establishes mathematical correctness and range safety for the exact
BLS12-381 MatrixNoK arithmetic selected by CatBLST. It does not claim that Lean
proves the C++ compiler, the x86 ISA, microarchitectural constant-time behavior,
or automatic equivalence between arbitrary future C++ edits and the model.

Those boundaries are connected by explicit assumptions, exact source/constant
hashes, code-generation inspection, native tests, and differential tests.

`MATRIXNOK_CONSTRUCTION.md` is the readable specification of the selected
construction, its representation invariants, exact bound ledger, call graph,
and deliberately duplicated performance schedules. A correspondence theorem
must agree with both this specification and the mapped executable source.

## Tooling

Use Lean 4 with a pinned Mathlib toolchain. The proof tree will live under
`proof/avx/` and must build without network access after dependency setup.
Generated proof inputs must be deterministic and reviewed; production C++
never depends on Lean or a runtime proof artifact.

## Layers

### P0: Lane semantics

Model one eight-lane vector as a fixed function from lane index to integer.
Specify and prove the exact operations used by the backend:

- add, subtract, bitwise masks, shifts, lane permutations, and comparisons;
- `madd52lo` and `madd52hi` as low/high 52-bit products accumulated modulo
  `2^64`; and
- any masked or shuffled composition used by change-base and multiplication.

The ISA semantics of the corresponding `_mm512_*` intrinsics are axioms at
this layer. CatBLST's `_cc512_*` implementation is the executable differential
oracle for those axioms.

### P1: Fixed constants and bounds

For the checked-in BLS12-381 configuration, prove:

- modulus dimensions and bit widths;
- pairwise coprimality and required modular inverses;
- fixed permutation tables are in range and bijective where required;
- correction, rotation, and Montgomery constants satisfy their defining
  congruences; and
- every intermediate respects its declared signed/unsigned bound and the
  52-bit/64-bit lane assumptions.

### P2: RNS and MatrixNoK

Prove conversion to RNS, RNS Montgomery multiplication, reduction, and
conversion back represent the same integer modulo the BLS12-381 prime.

For MatrixNoK specifically, prove that the fixed cyclic permutation and QR
decomposition implement the required r1 change-base map without the omitted
`k` correction under the checked constants and bounds.

### P3: Field tower

Prove the retained Fp operations against arithmetic modulo p, then prove Fp2
and Fp12 formulas against their defining extension polynomials. Include exact
alias-independent input/output specifications and zero/identity cases.

### P4: Curves and pairing

Prove the retained curve formulas under their documented preconditions, then
Miller-step and final-exponentiation algebra. Public validation and subgroup
contracts remain the API audit's responsibility; the proof states them as
preconditions rather than silently widening operations.

## C++ correspondence

Each proved kernel receives a correspondence record containing:

- C++ path and exact commit;
- Lean declaration and theorem;
- hash of the fixed constants consumed by both sides;
- operation-by-operation mapping from vector statements to lane semantics;
- focused portable and native differential test names; and
- code-generation evidence for the selected native build.

A C++ change that affects a mapped expression, constant, bound, or operation
reopens the corresponding theorem record. Passing tests alone cannot preserve a
formal claim after the mapped source changes.

## Constant-time evidence

The functional proof does not establish constant-time execution. Secret-bearing
kernels additionally require:

- source-level fixed control flow and address review;
- optimized native assembly inspection;
- structural checks against scalar-indexed loads and secret-dependent branches;
  and
- later statistical leakage testing as supporting evidence.

## Initial theorem milestone

The first accepted milestone is deliberately narrow:

1. exact semantics of the used IFMA lane operations;
2. all BLS12-381 MatrixNoK constants and permutation invariants;
3. one complete base-field multiplication/reduction round trip; and
4. correspondence to the side-by-side C++ candidate with local emulated and
   Peach native differential tests.

Only after that slice is stable should the proof expand into the tower and
pairing stack.

Current status:

- items 1 and 2 are machine checked;
- the P0/P1 correspondence rows are accepted at the exact hashes in
  `proof/avx/CORRESPONDENCE.csv`; `P0_P1_REVIEW.md` records the manual source
  mapping, focused native comparisons, constant regeneration, and remaining
  compiler/ISA trust boundary;
- item 3 is machine checked by `selected_base_field_pipeline`,
  `CanonicalImport`, and `CanonicalExport`: repaired import produces a
  representative in `[0,2p)`, selected multiplication preserves the internal
  `M1` scaling, and generated ConvertFrom removes that scaling to return the
  canonical field product;
- item 4 is accepted for the base-field P2 rows at the exact hashes in
  `proof/avx/CORRESPONDENCE.csv`. The Lean model reconstructs the actual split
  IFMA streams for R1, shifted R2, ConvertTo, and ConvertFrom; validates the
  selected reducers; and proves the generated quotient and all mapped
  accumulators fit their machine words before recombination. P0's compiler and
  x86-semantic boundaries remain explicit.

The negative-R2 analysis remains checked as a boundary result: the inherited
R2 schedule is not correct for every abstract representative in `(-p,0)`.
`AUD-221` repaired the reachable source of such values by changing only
ConvertTo to floor estimates. The generic estimate is exact or one low, and
the former exact witness now selects `p-t`, lies below `2p`, and preserves its
canonical value. Commits `e01ee59` and `f91268d` close the reachable
base-field P2 composition and its C++ correspondence without changing C++;
P3 is now in progress. The source-shaped selected Fp2 formulas and the exact
generic delayed reducer branch are machine checked, including the tight
`8m²` coordinate and `24m²` combined bounds as well as C++'s conservative
`832m²` type allowance, `<209m` intermediate, and `<2m` output. The
variadic strided layout and pointwise batched R1 recurrence are checked. An
abstract uncompleted Fp2 expression can select `-1`, and Lean proves that R2
preserves that boundary value. The selected C++ path does not feed this value
directly to R2: its negative static range causes `complete()` to add the
generated wide offset first. The semantic completion is proved to produce a
representative in `[0, 22p)`, with exact R2 correction over that whole range.
Commits `fa2af08` and `aa206e8` close the former source-composition boundary:
both generated offset bases, raw multiply/square producers,
M1 reduction, R1 change-base, M2 reduction, and R2 correction are composed in
one source-shaped theorem packet. `AVX-P3-001` is accepted; Fp12 and later
tower layers remain pending. No C++ arithmetic defect was found.

`AVX-P3-002` is the general-Fp12 layer and is intentionally narrower than the
whole `fp12.h` file. Its first source-shaped packet proves the exact
six-consecutive-Fp2 layout, the `v^3 = 1 + u` / `w^2 = v` tower,
general multiplication and optimized squaring. The semantic 22-product sum is
88 units. Direct unsigned subtraction contributes at most 40 units and is
dominated by the `[80m^2,81m^2)` generated offset. General multiplication is
below `169m^2`; square coordinates that double the offset-bearing expression
are below `242m^2`. The actual general-field boundary is MatrixNoK
`MAX_ADD=800`; `ready<200>` belongs only to a cyclotomic helper. All sixteen
generated offset lane encodings are checked.
M1 reduction, R1 change-base, M2 reduction, the `<18p` completed range, and
R2 correction are composed without assumed intermediate residues. The twelve
literal scalar multiply and square schedules now mechanically reconstruct the
Fp2-level formulas. Their exact coordinate-wise zero/one/two-offset
decomposition, unsigned lower bounds, and `4608p^2`/`5760p^2` RNS ceilings
compose with the generic completed-representative theorem, so `AVX-P3-002` is
accepted at the recorded source hashes. Frobenius maps form the next packet;
cyclotomic square and GT exponentiation remain a distinct narrow-domain packet.

`AVX-P3-003` now accepts Frobenius powers one through three at the literal
selected source schedules. A checked Lucas factor tree establishes primality of
the target modulus without an imported primality axiom. The map-1 theorem is
universal over the concrete Fp12 tower; coefficient-parametric source templates
and finite basis certificates prove map 2 and map 3 are its second and third
iterates, hence `x^(p^2)` and `x^(p^3)`. The generated-constant layer checks all
ten canonical values, all 160 physical lanes, exact `to_mont_exact` rotation,
standard P2 encodings, and the canonical-to-MatrixNoK representative relation.
The source-specific physical ledger is `4m^2`, well inside `MAX_ADD=800`, and
the proof includes the literal `standard_negate`/generated-`40p`-offset path and
specializes the accepted completed reducer.

`AVX-P3-004` accepts the policy-selected split cyclotomic square at proof commit
`256b980`. The exact scalar source predicate is converted to an equal flat-array
relation without projecting a normalized Fp12 product. Twelve independent
polynomial theorems prove that the selected result differs from the literal
general square by a fixed scaling of that relation's residual. This proves the
shortcut for precisely the source cyclotomic relation, including zero, without
claiming general Fp12 or narrowing the theorem to GT. The generated `1/3` and
`-2/3` values, all 32 physical lanes, the offset's semantic zero, and both split
bound ledgers are checked. The pre-reduction maxima are 28/888 for the first
M1/RNS split and 101/2041 for the second, below 800/5760; post-reduction factors
three and six fit `CommonElement<6>`. GT exponentiation and the relevant final-
exponentiation formulas were the remaining P3 work and must use their exact GT,
cyclotomic, or Miller-result preconditions.

`AVX-P3-005` now accepts the final-exponentiation formula at proof commit
`21887f9`. The concrete general-Fp12 conjugation is proved to equal Frobenius
six. A literal group model then preserves the repeated-square loop, all
`mul_n_sqr` constants, signed negative-`z` convention, easy part, and every
hard-part local. The fixed BLS12-381 parameter identities prove that the result
is exactly the repository's documented cubed pairing exponent
`3 * ((p^12 - 1) / r)`. The exact source ledger is 35 multiplications, 315
cyclotomic squares, nine conjugations, Frobenius counts 1/2/1, and one inversion;
no unreduced expression crosses an accepted component boundary. The theorem's
domain is a valid nonzero Miller-result unit with a correct supplied inverse.
At that checkpoint the concrete inverse decomposition remained an explicit
residual boundary, not a weakened final-exponentiation contract. It is now
accepted separately by `AVX-P3-008`; the variable-time base reciprocal provider
and P0 lowering remain external.

`AVX-P3-006` accepts the public GT binary and unsigned GLV4 paths at proof
commit `82f26c1`. The exact source-shaped split proves canonical reduction,
radix-`z` recomposition, four `<z` digits, and the odd-Frobenius conjugation
basis. A literal sixteen-entry table theorem feeds the 64-round CT scan and the
VT skip state machine; both compute the same reduced scalar exponent. The
binary CT and basic VT loops are also proved, and the public dispatch theorem
models the preceding little-endian `nbits` mask and the 129-through-256 GLV4
threshold. The physical row fixes table size, operation counts, selection
volume, and `CommonElement<6>`. The remaining P3 work is now narrowly the
balanced signed `avx_glv4_compile` bytecode plus its trusted executor, including
the schedule-length bound used by FFT. It is not covered by the unsigned public
path theorem merely because differential tests compare them.

Proof commits `51462b9` and `effdb3d` complete `AVX-P3-007`. The first proves
the four source-exact rounded divisions, bounded tail, `z^4` carry fold,
Frobenius-basis reconstruction, and signed 64-bit digit bounds. The second
proves exact sign-mask and magnitude encoding, highest-bit discovery, opcode
nibbles, fifteen-square chunking, trailing runs, compiled-executor semantics,
and the canonical FFT-scalar corollary. Valid 256-bit inputs require no more
than 63 of the 70 schedule bytes, 62 executor squares, 62 executor products,
and 74 products including table construction and normalization. The C++
Montgomery-to-canonical `Fr` bridge and common P0 compiler/ISA boundary remain
recorded assumptions rather than silently enlarged Lean claims.

`AVX-P3-008` accepts the nonzero Fp12 inverse decomposition at proof commit
`5042c1b`. The source-shaped theorem takes the quadratic norm to Fp6, forms the
cubic cofactors and Fp2 norm, takes the Fp2 norm to the base field, consumes one
supplied correct base reciprocal, and proves the lifted value multiplies the
input to one. It assumes exactly the selected inverter's mathematical contract;
zero and the provider's variable-time implementation remain outside the claim.

P3 is complete: `AVX-P3-001` through `AVX-P3-008` are accepted at exact source
and generated-constant hashes. P4 composes this accepted P2/P3 tower rather
than reconstructing it.

P4 is accepted through proof commit `4dbb3cd`. Direct and prepared doubling
and mixed addition now share proved point/line formulas. The literal sparse
Fp12 line multiplier is equivalent to accepted dense multiplication and has a
checked 72-product, 2880-unit envelope. The fixed
`1 + (add + n*dbl)` schedule `[2,3,9,32,16]` is exactly the 68-operation
binary schedule for `|x| = 0xd201000000010000`; its point invariant and final
Frobenius-six conjugation are composed.

The physical packet proves the required four-product `16m^2` bound, composes
the exact generated completion and both reducers at that bound, checks every
direct/prepared reduction input and the strided batch reconstruction, and
closes exact loop ledgers. The final theorem starts from the literal flat
square and sparse-multiply schedules—including the omitted first identity
square—and reaches the abstract binary Miller loop. Focused Lean and three
portable-Debug scalar differentials pass. P0 compiler/ISA lowering remains an
explicit external boundary rather than a Lean claim.

Production integration at `5c4dea1` created `AVX-P4-002`; Peach native
integration at `6ce0d44` accepts its exact-hash structural binding. The
accepted P4-001 theorem is the retained Miller construction; it must not be
generalized into a proof of every public point operation. In-review
`AVX-P4-003` proves the source-shaped complete-add/double field formulas and
curve preservation, exact GLV/GLS scalar splits, and the selected constant-time
Booth and variable-time wNAF executors. The wNAF proof includes backward
minimum-predecessor bounds that discharge the exact 129/65-slot capacities.
`CurveEndomorphism.lean` additionally proves the selected G1/G2 constants,
the exact beta/minus-psi maps, curve preservation, `sigma³ = 1`, `psi⁶ = -1`,
G1 complete-add compatibility, and the unsigned-basis sign bridges.
`CurveSubgroup.lean` proves the abstract subgroup equivalences while exposing
full curve cardinality and prime-subgroup generator span as named premises;
the BLS polynomial, both cofactors, and the decisive gcds are checked in Lean.
The G2 semilinear projective-additivity bridge is also proved: conjugation is
packaged as a ring homomorphism, diagonal scaling is proved against the source
complete-add formula, and their composition differs only by a proved-nonzero
common homogeneous scale. `CurveGroupBridge.lean` supplies the first concrete
group carrier: selected G1 is instantiated as a nonsingular Mathlib
Weierstrass curve over `ZMod p`; nonzero on-curve homogeneous representatives
map to its projective `AddCommGroup`; CatBLST equivalence is proved equivalent
to projective equality; and zero and negation commute with the map. `AUD-227`
also proves the selected Fp2 ring is a field by showing its norm cannot vanish
for a nonzero pair when `p mod 4 = 3`. The shared carrier bridge is instantiated
for the selected G2 twist as well, including projective equality, zero, and
negation. `CurveGroupAdd.lean` now proves bihomogeneous normalization and the
full complete-add/group-add identification for arbitrary valid selected G1/G2
representatives. It covers infinity, doubling, inverse, and generic affine
addition, and discharges the selected curves' no-affine-2-torsion condition by
finite-field cube criteria. `CurveBinaryMul.lean` proves that a logarithmic
evaluator built solely from the exact complete-add owner is Mathlib projective
scalar multiplication. `CurveGenerators.lean` checks the standard G1/G2
coordinates, proves both generators are killed by `r`, and proves the selected
literal endomorphisms give the radix-`z²`/radix-`z` bases. Finally,
`CurveDispatch.lean` composes public masking and modulo-`r` preparation with
the exact Booth and wNAF schedules over both `r`-torsion carriers and the
conventional cyclic G1/G2 subgroups. `AVX-P4-003` is accepted at proof commit
`f1cfde9`.

The multiplication proof does not require a point count: its operand domain
is the conventional cyclic subgroup generated by the standard point. The
different claim that the runtime relation-based `in_g1`/`in_g2` predicates
recognize exactly that domain is now substantially reduced by
`CurveMembership.lean`. The literal coordinate maps are proved well-defined
and additive on the selected projective quotient, and the exact source
relations satisfy
`standard cyclic subgroup ⊆ accepted relation ⊆ r-torsion` without a full
curve-cardinality assumption. The exact accepted relation kernels are now
`ZMod r` modules with the selected unsigned GLV/GLS basis identities, and all
four CT/VT multiplication schedules are proved correct on those carriers.
Proof commit `65b54f6` additionally supplies a Lucas certificate for the
primality of `r`, proves both selected generators have exact order `r`, and
proves each accepted kernel has at least `r` elements. Exact equality with the
generator span is now equivalent to the single upper bound that the accepted
kernel has at most `r` elements. Proof commit `f6c251b` establishes this bound
for G1 from the generic `2p+1` affine-point bound and `r² > 2p+1`, yielding
exact equality between the accepted G1 kernel and the selected generator span.
Proof commit `475221a` proves from `r³ > 2p²+1` that the accepted G2 kernel
has cardinality exactly `r` or `r²`. The rank-two alternative is excluded by
the corrected BLS12 subgroup theorem in El Housni, Guillevic, and Piellard,
ePrint 2022/352, Proposition 5. Proof commit `a8075e7` checks its concrete
BLS12-381 hypotheses: the signed seed and base-field polynomial, the published
G2 cofactor polynomial, and `gcd(chi(u), c2) = 1`. `AUD-231` is validated with
the general curve-family theorem and twist-order identification explicitly
external. Re-proving the BLS12 twist order from foundational elliptic-curve
theory in Lean is outside the accepted proof scope.
The production alias validator closes source selection ambiguity but is not a
mathematical theorem. MatrixNoK-to-scalar observation through byte weaving and
scalar REDC also remains outside Lean when used by equality, infinity
observation, validation, and secret normalization. `AUD-259`/`AUD-263`
separately prove its machine envelope and record branchless MSVC carry lowering;
scalar REDC/IFMA semantics and native Clang lowering remain external.

The `AUD-229` constant-initialization packet at `fe09e26` updates four mapped
owner headers while leaving the arithmetic schedules intact. Peach Clang
21.1.8 Release reacceptance passes 35 focused tests and all five audit gates;
guard/atexit relocations are absent, the selected arithmetic probe is
byte-identical, and hot owner stack frames are unchanged. The five affected
rows are accepted at their new hashes. The independently regenerated
curve-literal header is an additional exact-hash input to `AVX-P4-002`; the
original proved constants artifact is unchanged. `AUD-230` makes the ledger's
51 mapped hashes LF-canonical so the same proof boundary validates on Windows
CRLF and Peach LF checkouts.

`AUD-264` removed the comparison-only direct Fp2 microkernel. The selected
tower and the Fp2 proof now have one C++ owner:
`FP2Ring<FpRingCandidate>`.

## Semi-freeze discipline

Once the preparation packet records exact hashes, changes to mapped executable
expressions, constants, bounds, permutations, or operation schedules reopen
the affected correspondence record. Performance-shaped scalar/batch schedules
remain separate C++ implementations and refine one shared Lean specification;
the proof effort is not a reason to merge hot kernels behind an abstraction.
