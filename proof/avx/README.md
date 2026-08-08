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

This reopens the full-reachability part of initial milestone item 3. The P2
correspondence rows are in review: source hashes are unchanged from `adf2ed5`,
the deterministic constant check passes, and the shared scalar/batch schedule
has an explicit Lean refinement theorem. Acceptance still requires the normal
review of the C++ transcription and existing native/portable evidence.
`CanonicalExport` closes the remaining mathematical boundary. It reconstructs
the eight wide radix-`2^52` ConvertFrom coefficients, proves the generated
matrix and correction represent the exact bounded R2 representative, proves
the 128-bit quotient and both 64-bit IFMA accumulators cannot wrap, and proves
both canonical import/export and canonical import/multiply/export. The final
theorem is conditional only on the ordinary reduced-lane postconditions proved
by the reducer modules. P3 field-tower constants and formulas are not claimed
by this base-field theorem.

The generated snapshot and Python regeneration remain binding and independent
executable evidence, respectively; Lean checks the corresponding formulas
rather than accepting either tool as the mathematical proof.

## Resource policy

Toolchain installation and Mathlib cache retrieval are explicit setup steps.
Proof builds run serially unless separately approved. They do not build
CatBLST, CUDA, benchmarks, BTx/BTE, or the portable-AVX Release product.
