# CatBLST AVX Arithmetic Proof

This is the machine-checked proof project for the selected BLS12-381
MatrixNoK arithmetic. It implements the layers and claim boundary in
`docs/avx/PROOF_PLAN.md`; `docs/avx/MATRIXNOK_CONSTRUCTION.md` is the readable
specification of the exact selected construction and bound ledger.

The project pins Lean and Mathlib to `v4.32.1`. Updating either pin is a
reviewed dependency change and does not occur implicitly.

## Local setup

Install Elan, then run these commands from `proof/avx`:

```text
lake update
lake exe cache get
lake build
python ../../scripts/audit/generate_matrixnok_lean.py --check
```

Elan selects the exact version in `lean-toolchain`; do not install or select a
different system Lean. `lake-manifest.json` locks Mathlib's transitive
dependencies, while the downloaded `.lake/` cache remains local and ignored.

## First milestone

The first proof slice is intentionally narrow and follows the implementation
dependency order:

1. define 64-bit word, 52-bit digit, IFMA-low, IFMA-high, accumulation, and
   eight-lane permutation semantics;
2. prove the low/high product split and the no-unintended-wrap bounds used by
   the fixed r1 kernel;
3. import the exact generated BLS12-381 constants through a deterministic
   source generator and prove their defining modular relations;
4. prove the fixed cyclic r1 schedule equals the checked-in matrix map;
5. prove both RNS Montgomery reductions and conversion round trips; and
6. prove `FpRingCandidate::multiply` represents scalar multiplication modulo
   the BLS12-381 prime and returns both residue halves below `2m_i`.

No theorem may be marked complete by an executable test, `native_decide`, or a
generated equality alone when the corresponding general arithmetic or bound
claim remains assumed. Computation can discharge fixed finite facts; reusable
arithmetic properties require explicit theorems.

## Correspondence rule

`CORRESPONDENCE.csv` is the review boundary between Lean and C++. Each accepted
row records the exact implementation commit, source hash, constant hash, Lean
declaration, theorem, and native/portable evidence. A mapped C++ or constant
change reopens the row even if tests still pass.

The proof does not establish compiler correctness, x86 ISA correctness, or
microarchitectural constant-time behavior. Those assumptions are tracked by
the native code-generation and constant-time evidence described in the main
proof plan.

## Current proof status

P0 now has a checked-in, buildable foundation:

- `CatBLST.AVX.Lane` defines unsigned 64-bit lanes, eight-lane vectors,
  wrapping arithmetic, shifts, masking, variable permutation, lane-zero
  extraction, and broadcast;
- `CatBLST.AVX.IFMA` defines low/high IFMA semantics and proves the exact
  radix-`2^52` product split, both digit bounds, modulo-`2^64` accumulation,
  and generic no-wrap accumulation rules; and
- `CatBLST.AVX.PortableIFMA` mirrors the retained 26-bit-half C++ emulation and
  proves its low and high results equal the abstract IFMA operations for every
  pair of 52-bit inputs.

The P0 correspondence rows are accepted at the exact hashes in
`CORRESPONDENCE.csv` after the structured review in
`docs/avx/P0_P1_REVIEW.md`. The arithmetic transcription is machine checked;
reviewed source mapping, C++ unsigned-integer semantics, compiler lowering,
and the x86 IFMA specification remain explicit external boundaries. They are
intentionally not represented as Lean axioms that could make a later theorem
appear stronger than it is.

The first two P1 constant slices are accepted at their exact hashes. The
deterministic `generate_matrixnok_lean.py` snapshot binds every generated C++ literal and
its source hash. `CatBLST.AVX.MatrixNoK.Constants` proves the selected
dimensions, modulus ranges, oddness, pairwise coprimality, exact basis
products, dynamic range, rotation and target inverses, Montgomery factors,
conversion multiplier relations, fixed-table shapes/ranges, and cyclic R1
permutation relation. `CatBLST.AVX.MatrixNoK.Regeneration` then defines the
CRT generator independently of the tables, verifies explicit modular-inverse
witnesses, regenerates all four matrix/correction/quotient families, bounds
the fixed-point floor and ceiling errors, and proves the eight-lane quotient
accumulator cannot wrap 128 bits for 52-bit inputs and 64-bit estimates.

The nonnegative selected base-field P1/P2 slice is machine checked. The proof covers
the exact product, correction-free R1, fused M2, shifted R2, and Montgomery
reduction bounds; the cyclic R1 source schedule; fixed-quotient correctness;
the absence of unintended 64-bit accumulator wrap; and a complete
`selected_base_field_pipeline` theorem. Inputs may use either standard lane
representative below `2m_i`. Both output bases remain below `2m_i`, encode the
same representative, and that representative is the field product modulo p.
The selected raw QR representative is below `18p`, safely inside the C++
base-redundancy allowance of `40p`.

Canonical import exposed `AUD-221`: VROOM's one-high ConvertTo quotient could
select a negative raw lift even though R2 and ConvertFrom require a
nonnegative representative. Commit `5f0399e` changes only ConvertTo to floor
fixed-point estimates. `Conversion` proves the estimate is exact or one low;
`CanonicalImport` proves the repaired representative lies in `[0,2p)` and
preserves the imported field value; and `CanonicalImportRepair` proves the
exact former counterexample now selects its positive complement and satisfies
the R2 range. `NegativeR2Counterexample` retains only the abstract warning that
R2 is not correct for every negative representative; repaired canonical import
no longer reaches that domain.

The base-field P2 correspondence rows are now accepted at exact source hashes.
The review maps the actual split `madd52lo`/`madd52hi` streams for R1, shifted
R2, ConvertTo, and ConvertFrom to exact reconstruction theorems. Canonical
ConvertTo is proved below `16m²`, comfortably inside the unchanged C++ type
allowance `800m²`, and the proof invokes the selected reducer with the actual
`800/200` template bounds rather than silently narrowing the implementation.

`CanonicalExport` now closes the full reachable path without assumed
intermediate residues. The theorem
`actual_canonical_import_selected_multiplication_export` starts with canonical
radix-`2^50` digits, follows the actual ConvertTo reducer and R2-derived M1
path, performs selected multiplication, reconstructs ConvertFrom's eight wide
radix-`2^52` coefficients, and returns the canonical field product. The source
review found no further C++ arithmetic defect. `AUD-263` independently proves
the byte-weaving machine envelope through REDC, and `AUD-259` records clean
MSVC carry lowering after repairing a secret-dependent branch. Those facts are
still outside Lean; scalar REDC semantics, x86 IFMA semantics, and native Clang
lowering remain explicit external boundaries. P3 field-tower constants and
formulas are not claimed by this base-field theorem.

P3 begins with the sole retained selected owner,
`FP2Ring<FpRingCandidate>`. `Fp2.lean` maps the source-shaped `prep_left`,
multiply, square, multiply-accumulate,
multiply-sub-accumulate, nonresidue, and `3b` formulas to
`Fp[u]/(u^2 + 1)`. Standard preparation and negation retain the physical
`<2m` component bound, so each two-product coordinate is below `8m^2` and its
exact R1-augmented value is below `24m^2`. The proof separately preserves the
implementation's conservative type path: C++ recasts to `MAX_ADD=800`, adds
the type-level R1 allowance `32`, obtains `832`, returns `<209m` after
`reduce_wide`, then `reduce_auto<2>` selects mask reduction and returns `<2m`.
The variadic `x...`/`y...` layout and pointwise batched R1 recurrence are also
checked. `Fp2SignedBoundary` shows that the uncompleted abstract signed formula
can select correction-free representative `-1`, and proves the existing R2
schedule preserves that boundary value. It does not claim that selected C++
reaches the value: `FpRingCandidate::complete()` first adds the generated wide
negative-range offset. `Fp2.lean` now proves that the corresponding semantic
completion puts every selected two-product coordinate in `[0, 22p)` and that
R2 is exact throughout that range. Commits `fa2af08` and `aa206e8` then prove
the generated `max_neg_offset` lanes in both bases, compose the physical
`complete()`/M1/R1-change-base/M2/R2 path, and prove the raw real/imaginary
multiply and square lane producers and their tight bounds. `AVX-P3-001` is
accepted at the exact source hashes; Fp12 and later tower layers remain P3
work. The reviewed high/low IFMA split and native/compiler semantics remain
the common P0 boundary rather than an Fp2-specific assumption.

`AVX-P3-002` accepts the selected general-Fp12 multiplication and squaring
packet. `Fp12.lean`
fixes the exact six-Fp2 storage order, defines the extension relations
`v^3 = 1 + u` and `w^2 = v`, and proves that the flattened multiplication and
optimized squaring schedules implement that tower. It also checks the largest
22-product coordinate. The selected general-field path uses
`ring.batch_reduce_expand` with `MAX_ADD=800`, not the cyclotomic helper's
`ready<200>` boundary. A single-offset multiplication coordinate is below
`169m^2`; doubled offset-bearing square coordinates are below `242m^2`.
All sixteen generated offset lanes lie in `[80m^2,81m^2)` and carry the proved
encoding. The selected M1 reduction, R1 change-base, M2
reduction, `<18p` completed representative, and exact R2 correction are now
composed in one theorem. `Fp12Raw.lean` mechanically transcribes all twelve
literal multiply and square polynomials and proves that adjacent scalar outputs
reconstruct the Fp2-level formulas. It also proves the exact zero/one/two-offset
decomposition, unsigned subtraction safety, and the respective `4608p^2` and
`5760p^2` RNS ceilings. These feed the generalized completed-representative
theorem, closing the physical reduction composition. Frobenius and cyclotomic
identities remain separate because they add generated-constant and input-domain
obligations.

`AVX-P3-003` accepts the three selected general-Fp12 Frobenius maps. The proof
constructs the concrete twelve-coordinate tower as a commutative ring, checks
the BLS12-381 base modulus with a fully internal Lucas certificate, and proves
the literal map-1 schedule equals exponentiation by `p` for every Fp12 value.
The literal map-2 and map-3 schedules are coefficient-parametric linear maps;
their checked twelve-basis certificates compose map 1 into universal
`x^(p^2)` and `x^(p^3)` theorems without normalizing the 381-bit coefficients
inside every symbolic goal. All ten canonical Frobenius coefficients and all
160 generated M1/M2 lanes are tied to `to_mont_exact`, the standard-element P2
contracts, and their canonical field meanings. A source-specific physical
ledger proves the largest reduced coordinate is below `4m^2` versus the literal
`MAX_ADD=800` boundary. The unchanged/negated coordinates include an exact
model of `standard_negate`, its generated `40p` zero offset, and the selected
one-subtraction `reduce_auto<2>` path. The packet then specializes the accepted
general-Fp12 completed-reduction chain. Frobenius is valid on general Fp12;
cyclotomic square and GT exponentiation retain their narrower domains.

`AVX-P3-004` accepts the selected split cyclotomic square at proof commit
`256b980`. The proof first shows that the scalar source predicate
`Frob^2(x) = Frob^4(x) * x` is exactly equivalent to a literal twelve-lane
predicate. It then proves an unconditional residual identity: optimized square
minus ordinary square is a fixed diagonal scaling of that predicate residual.
All twelve coordinates therefore reduce to the ordinary Fp12 square under the
exact source predicate. Zero remains in this domain, matching C++; GT is the
strictly narrower nonzero order-r subgroup and is not assumed. Both generated
`1/3` and `-2/3` constants and all 32 M1/M2 lanes are tied to their exact
`to_mont_exact` encodings and canonical field values. The two physical split
ledgers peak at 28 M1 units / 888 RNS units and 101 M1 units / 2041 RNS units,
respectively, below the selected 800/5760 reduction boundaries. GT
exponentiation remained the final P3 packet at that checkpoint.

`AVX-P3-005` accepts the complete final-exponentiation formula at proof commit
`21887f9`. The concrete conjugation schedule is proved to be Frobenius to the
sixth power for every selected Fp12 value; it is treated as inversion only in
the nonzero cyclotomic group used after the easy part. The Lean schedule keeps
the literal repeated-square loop, all five `mul_n_sqr` stages, four
`raise_to_z` calls, the additional `raise_to_z_div_by_2`, every named hard-part
local, and the Frobenius-1/2/3 placements. It proves the exact inherited BLST
convention `f^(3 * ((p^12 - 1) / r))`, rather than silently substituting the
uncubed textbook exponent. The physical ledger records 35 Fp12 multiplies, 315
cyclotomic squares, nine conjugations, four Frobenius maps, and one inversion;
the header serially composes accepted operations and never carries an
unreduced expression across calls. At that checkpoint the remaining explicit
boundary was that a valid nonzero Miller result and its inverter output form
the modeled unit. `AVX-P3-008` now proves the concrete tower decomposition from
exactly one correct supplied base reciprocal; the variable-time reciprocal
provider and P0 lowering remain external.

`AVX-P3-006` accepts the public GT multiplication schedules at proof commit
`82f26c1`. The proof preserves the exact unsigned radix-`z` split performed
after scalar preparation, proves all four digits are below `z < 2^64`, and
checks the source basis `[1,-p,p^2,-p^3]` is `[1,z,z^2,z^3]` modulo `r`.
It constructs the literal sixteen-entry table and proves the descending
64-round CT state machine and the unsigned VT leading/zero-window skip machine
both return the canonical scalar exponent. The binary CT fallback and fixed
256-bit basic VT scan are proved from their literal square/select/skip loops;
the public CT dispatch theorem includes byte masking and the 129-through-256
GLV4 threshold. The exact CT ledger is 11 table multiplications, 64 squares,
64 loop multiplications, one normalization multiplication, 960 table-candidate
comparisons, and 11,520 coefficient selections. These theorems are strictly
GT-domain results composed through accepted Fp12 operations. The balanced
signed schedule compiler/executor used by FFT remains a distinct P3 row rather
than being inferred from the public unsigned executor.

`AVX-P3-007` is accepted by proof commits `51462b9` and `effdb3d`. The first
commit transcribes the source's remainder-centering rule, four signed division
identities, `z^4` carry fold, odd-coordinate sign change, and strict signed
64-bit bounds. The second models the exact sign mask, absolute magnitudes,
high-to-low MSB discovery, high/low opcode nibbles, fifteen-square chunking,
trailing runs, and trusted executor. The composed theorem returns the original
scalar exponent for every input below `2^256`. Physical theorems show valid
inputs emit at most 63 of the owner's 70 bytes, execute at most 62 squares and
62 post-initialization multiplications, and require at most 74 multiplications
including the eleven-entry table construction and final normalization. The
FFT corollary records the canonical `Fr < r` precondition used by prepared
roots and inverse scales. The C++ `Fr::from_mont` representation bridge and P0
compiler/ISA lowering remain explicit non-Lean boundaries.

`AVX-P3-008` is accepted by proof commit `5042c1b`. It models every source stage
from the Fp12 norm through Fp6 cofactors, Fp2 and base norms, and the lift back
through the tower. Assuming exactly that the supplied base reciprocal satisfies
`norm * reciprocal = 1`, the final theorem proves multiplication by the
constructed value yields Fp12 one. Zero remains outside the contract.

All eight P3 correspondence rows are now accepted. P4 composes those accepted
field-tower results into the complete pairing implementation.

`AVX-P4-001` begins at proof commit `8e0dc2d` and is accepted through
`4dbb3cd`.
The exact direct/prepared doubling and mixed-addition formulas preserve the
curve, and their tangent/chord lines contain the required points. The literal
`flat_mulxy00z0` schedule is proved to be dense Fp12 multiplication by the
sparse `(c, x, 0, 0, y, 0)` line value. Its 72-product ledger and 2880-unit
RNS ceiling are checked independently of the dense kernel. The source's
68-operation schedule has exactly 63 doublings and five additions, reaches
`0xd201000000010000`, equals the ordinary binary schedule, preserves the curve
through the fold, and ends with the accepted Frobenius-six conjugation.

The physical packet proves every accumulated Fp2 coordinate contains at most
four base products and fits `16m^2`, generalizes the complete M1/R1/M2/R2
reduction chain to that envelope, checks the strided batch reconstruction, and
records exact direct and prepared ledgers. The literal C++ square/sparse paths
are then composed with the algebraic operations, including the source's
omitted first `one^2`. `selected_literal_cpp_miller_loop_correct` is the final
finite-coordinate literal-to-binary-loop theorem; it does not itself decode
CatBLST's affine `(0,0)` infinity sentinel. The `AUD-277` extension models the
public infinity observation bits separately, proves both identity branches,
and proves that the encoded sparse identity line remains one through the exact
first-line, 67-operation tail, and final conjugation used by prepared
evaluation. Scalar export-based infinity observation and P0 compiler/ISA
lowering remain explicit external boundaries.

The focused root Lean build passed 3,082 targets, the Miller proof trust-word
scan is clean, and three portable-Debug scalar differentials passed:
G2 complete operations, sparse Fp12/final exponentiation, and complete Miller
pairing. The already recorded Peach native differential evidence remains
bound to the unchanged C++ hashes.

Proof SHA-256 values are `978e297d6bf61a6c4070829848d0dfeaa6881129878e951bd3263657b7719038`
(`MillerDouble.lean`), `5ac59d84fe1a6543239b9b334af4e3ce79f2414354ff87101d6911f5f36e083f`
(`MillerSparseMul.lean`), `dfeab10383b959f74b052120ae1fa7e628e56d28a73edc7c6cd768248b17a4f0`
(`MillerAdd.lean`), `a30648e2ee74f3aa4bade730e24967db841429d62adc3d63811285c99f8703e1`
(`MillerSchedule.lean`), `671b32356430c0d472bd42668476483eb23bd980b06c4e221469686a8884b7e7`
(`MillerPhysical.lean`), `4db61a4a23458ca201f8743da42a8f5774f69070a02206664d08cbe9abb24eff`
(`MillerReduction.lean`), and `c20fbb26197715668c644db44f012b593811f692a7f3c19e9dee89cfa960c499`
(`MillerEndToEnd.lean`).

The generated snapshot and Python regeneration remain binding and independent
executable evidence, respectively; Lean checks the corresponding formulas
rather than accepting either tool as the mathematical proof.

## Production integration correspondence

Commit `5c4dea1` selects the exact MatrixNoK theorem-owner types for production
AVX point, GT, pairing, and GT-FFT dispatch. `AVX-P4-002` binds every alias and
production owner at exact hashes and is accepted after the Clang 21.1.8 Peach
native Release build, 34 focused tests, three structural audits, and optimized
code-generation inspection at `6ce0d44`. The structural
`validate_matrixnok_dispatch.py` gate prevents a production owner from silently
falling back to the incumbent ring; it is a source-correspondence check rather
than a mathematical theorem.

The constant-initialization cleanup at `fe09e26` changes exact-hash owner files
without changing their arithmetic expressions. It also adds a separately
regenerated curve-literal artifact so generators and endomorphism owners need
no guarded first-use conversion. Peach Clang 21.1.8 Release reacceptance passes
35 focused tests and all five audit gates. Native guard/atexit relocations are
gone, the selected arithmetic probe is byte-identical, and hot owner stack
frames are unchanged. `AVX-P4-002`, `AVX-P4-001`, `AVX-P3-007`, `AVX-P3-006`,
and `AVX-P3-001` are accepted at their new hashes. `AUD-230` canonicalizes all
51 mapped hashes to LF so the exact ledger is independent of worktree line
endings. Their established Lean theorems remain unchanged.

The accepted P4 Miller theorem covers the retained Miller-specific point
formulas and loop. Accepted `AVX-P4-003` transcribes the public shared complete-add and
double schedules, proves curve preservation and the identity/inverse/doubling
checkpoints for selected Fp and Fp2, proves the exact G1 radix-`z²` and G2
radix-`z` splits, and proves both the fixed Booth and variable-time wNAF
executor schedules, including the exact 129/65-slot exhaustion bounds.
`CurveEndomorphism.lean` now proves the concrete beta and minus-psi constant
identities, exact coordinate maps, curve preservation, `sigma³ = 1`,
`psi⁶ = -1`, the G1 map's exact compatibility with complete addition, and the
sign bridge from BLST's subgroup relations to the unsigned `z²`/`z` executor
bases. `CurveSubgroup.lean` proves that each relation implies `[r]P = 0` from
the explicit full curve cardinality, and proves the converse from the explicit
generator relation and prime-subgroup spanning fact. It checks the BLS
polynomial factorization and both cofactor gcds in Lean rather than importing
them as arithmetic assumptions. The G2 semilinear map is now factored into a
ring-homomorphism step and a diagonal curve-isomorphism step; complete addition
after `psi` is proved to be one explicitly nonzero homogeneous scale of `psi`
after complete addition. `CurveGroupBridge.lean` now instantiates the selected
G1 equation as a nonsingular Mathlib Weierstrass curve over the proved prime
field. It maps every nonzero on-curve homogeneous representative to Mathlib's
projective `AddCommGroup`, proves that CatBLST homogeneous equivalence is exactly
projective equality on that domain, and proves the identity and negation maps
agree. `Fp2Field.lean` proves that the selected `Fp2Value` ring is a field: the
base modulus is `3 mod 4`, so `-1` is not a square, every nonzero pair has
nonzero norm, and conjugate-over-norm supplies an inverse. The same generic
bridge now instantiates the selected G2 twist over that field and proves its
projective equality, identity, and negation correspondence.
`CurveGroupAdd.lean` now proves the complete-add schedule is bihomogeneous,
normalizes arbitrary valid representatives, and identifies its output with
Mathlib group addition across infinity, doubling, inverse, and generic affine
cases. The generic theorem exposes the exact no-affine-2-torsion and
characteristic-not-two hypotheses. Finite-field cube criteria discharge those
hypotheses for both selected BLS12-381 curves, yielding concrete G1 and G2
group-addition theorems. `CurveBinaryMul.lean` proves a logarithmic evaluator
built from the exact complete-add owner equal to Mathlib projective scalar
multiplication. `CurveGenerators.lean` checks the standard BLS12-381 G1/G2
coordinates, proves both generators are killed by `r`, and proves the selected
beta/minus-psi maps give the required radix bases. `CurveDispatch.lean`
composes public masking and reduction with the selected Booth and wNAF
schedules over both `r`-torsion and the conventional cyclic G1/G2 subgroup
carriers. `AUD-227` and `AVX-P4-003` are therefore accepted at proof commit
`f1cfde9`. `CurveMembership.lean` then removes the ambient-cardinality premise
from the hard direction of the relation-based runtime membership proof. It
packages the literal homogeneous `sigma` and conjugating `psi` maps as
quotient-safe additive projective endomorphisms and proves, for both selected
curves, `standard cyclic subgroup ⊆ accepted relation ⊆ r-torsion`. The
relation-to-`r`-torsion step uses the concrete infinity/x=0 3-torsion kernels
and `gcd(z²+1,3)=1`. The exact accepted relations are also instantiated as
additive subgroup kernels and `ZMod r` modules. Their unsigned GLV/GLS basis
identities and all four selected CT/VT multiplication theorems hold for every
validator-accepted point, so the optimized arithmetic no longer depends on
the missing span identification. `ScalarOrderPrime.lean` and
`CurveGeneratorOrder.lean` add a Lucas certificate for `r`, prove the selected
G1/G2 generators have exact order `r`, and show their accepted relation kernels
have at least `r` elements. `CurveMembership.lean` then proves equality with
the standard generator span is equivalent to the upper bound
`card accepted ≤ r`. `CurveMembershipCardinality.lean` closes that bound for
G1 without importing the published exact curve order: every affine x-coordinate
has at most two points, so the ambient projective curve has at most `2p+1`
points, while a larger accepted `ZMod r` kernel would have at least `r² >
2p+1` elements. Over `Fp2`, `r³ > 2p²+1` proves the accepted G2 kernel has
dimension at most two and hence exactly `r` or `r²` elements. The remaining
rank-two exclusion is the published BLS12 criterion in El Housni, Guillevic,
and Piellard, ePrint 2022/352, Proposition 5, which repairs Scott's original
G2 proof. `CurveSubgroup.lean` checks the exact BLS12-381 seed and cofactor
polynomials and `gcd(chi(u), c2) = 1`; the paper supplies the general
curve-family theorem and published twist-order fact. `AUD-231` is therefore
validated with that external mathematical dependency. A foundational Lean
derivation of the BLS12 twist order remains deliberately out of scope.
`CurveBatchValidation.lean` proves that per-point curve validation followed
by the subgroup relation gives the same result as full individual validation.
It also proves that the G2 cofactor is coprime to coefficient differences
1 through 12, so each ideal 13-value row has at most one cancelling choice;
22 independent rows give `13^-22 < 2^-80`. The quotient cardinality is an
explicit premise, and the SHA-256 PRF model and physical AVX bucket-addition
correspondence remain separate proof boundaries.
Canonical export followed by scalar REDC and byte weaving remains an explicit
non-Lean observation boundary for equality, infinity observation, validation,
and secret normalization rather than a theorem silently inferred from
`CanonicalExport`. The independent machine-envelope validator and optimized
MSVC carry inspection recorded in `AUD-259`/`AUD-263` cover the concrete C++
boundary without extending the Lean theorem.

The production GT generator required a new fixed MatrixNoK literal.
`AVX-P1-003` and `GtGeneratorConstant.lean` prove that all 192 RNS lanes are
the exact standard encoding of the twelve canonical Fp coefficients and meet
the selected lane contracts. The statement that those coefficients are
CatBLST's chosen scalar GT generator remains an explicit reviewed constant
boundary plus scalar/AVX differential test. No broader generator-order theorem
is claimed.

## Resource policy

Toolchain installation and Mathlib cache retrieval are explicit setup steps.
Proof builds run serially unless separately approved. They do not build
CatBLST, CUDA, benchmarks, BTx/BTE, or the portable-AVX Release product.
