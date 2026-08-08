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
- item 3 is machine checked by `selected_base_field_pipeline`,
  `CanonicalImport`, and `CanonicalExport`: repaired import produces a
  representative in `[0,2p)`, selected multiplication preserves the internal
  `M1` scaling, and generated ConvertFrom removes that scaling to return the
  canonical field product;
- item 4 is in review against the unchanged `adf2ed5` C++ source hashes and
  existing portable/native evidence. The Lean export model also proves the
  generated quotient fits 128 bits and the non-shifted ConvertFrom high/low
  accumulators fit 64 bits before recombination.

The negative-R2 analysis remains checked as a boundary result: the inherited
R2 schedule is not correct for every abstract representative in `(-p,0)`.
`AUD-221` repaired the reachable source of such values by changing only
ConvertTo to floor estimates. The generic estimate is exact or one low, and
the former exact witness now selects `p-t`, lies below `2p`, and preserves its
canonical value. Commit `e01ee59` closes the reachable base-field P2
composition without changing C++; the next milestone work is the remaining
C++ correspondence review, followed by P3 field-tower formulas.

The comparison-only direct `matrixnok::Fp2RingCandidate` microkernel is outside
this initial scope. The selected tower uses `FP2Ring<FpRingCandidate>`.

## Semi-freeze discipline

Once the preparation packet records exact hashes, changes to mapped executable
expressions, constants, bounds, permutations, or operation schedules reopen
the affected correspondence record. Performance-shaped scalar/batch schedules
remain separate C++ implementations and refine one shared Lean specification;
the proof effort is not a reason to merge hot kernels behind an abstraction.
