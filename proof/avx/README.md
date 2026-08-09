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

The P0 correspondence row remains in review. The arithmetic transcription is
machine checked, but the proof still relies on reviewed source mapping, C++
unsigned-integer semantics, compiler lowering, and the x86 IFMA specification.
Those boundaries are intentionally not represented as Lean axioms that could
make a later theorem appear stronger than it is.

The first two P1 constant slices are also machine checked. The deterministic
`generate_matrixnok_lean.py` snapshot binds every generated C++ literal and
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
review found no further C++ arithmetic defect. Compiler lowering, x86 IFMA
semantics, and `from_vec`'s byte-weaving/REDC transcription remain explicit
non-Lean boundaries. P3 field-tower constants and formulas are not claimed by
this base-field theorem.

P3 now begins with the selected generic `FP2Ring<FpRingCandidate>` path, not
the isolated comparison-only `matrixnok::Fp2RingCandidate`. `Fp2.lean` maps
the source-shaped `prep_left`, multiply, square, multiply-accumulate,
multiply-sub-accumulate, nonresidue, and `3b` formulas to
`Fp[u]/(u^2 + 1)`. Standard preparation and negation retain the physical
`<2m` component bound, so each two-product coordinate is below `8m²` and its
exact R1-augmented value is below `24m²`. The proof separately preserves the
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
cyclotomic square and GT exponentiation retain their narrower domains and form
the next P3 packet.

The generated snapshot and Python regeneration remain binding and independent
executable evidence, respectively; Lean checks the corresponding formulas
rather than accepting either tool as the mathematical proof.

## Resource policy

Toolchain installation and Mathlib cache retrieval are explicit setup steps.
Proof builds run serially unless separately approved. They do not build
CatBLST, CUDA, benchmarks, BTx/BTE, or the portable-AVX Release product.
