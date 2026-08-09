/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.GtExponentiation
import Mathlib.Tactic

/-!
# GT exponentiation physical ledger

The arithmetic proof works in the GT exponent group.  This companion records
the fixed array sizes, literal loop bounds, and source-level operation counts
that distinguish the constant-time schedule from the public-scalar VT path.
-/

namespace CatBLST.AVX.MatrixNoK

structure GtGlv4CtCost where
  tableEntries : Nat
  tableBuildMultiplications : Nat
  frobeniusMaps : Nat
  conjugations : Nat
  loopIterations : Nat
  cyclotomicSquares : Nat
  loopMultiplications : Nat
  tableCandidatesPerIteration : Nat
  coefficientsPerCandidate : Nat
  normalizationMultiplications : Nat
  deriving DecidableEq

/-- Exact source ledger for `avx_cyclotomic_exp_glv4_ct`. -/
def gtGlv4CtCost : GtGlv4CtCost where
  tableEntries := 16
  tableBuildMultiplications := 11
  frobeniusMaps := 3
  conjugations := 2
  loopIterations := 64
  cyclotomicSquares := 64
  loopMultiplications := 64
  tableCandidatesPerIteration := 15
  coefficientsPerCandidate := 12
  normalizationMultiplications := 1

def GtGlv4CtCost.totalMultiplications (cost : GtGlv4CtCost) : Nat :=
  cost.tableBuildMultiplications + cost.loopMultiplications +
    cost.normalizationMultiplications

def GtGlv4CtCost.totalCoefficientSelections (cost : GtGlv4CtCost) : Nat :=
  cost.loopIterations * cost.tableCandidatesPerIteration *
    cost.coefficientsPerCandidate

theorem gt_glv4_ct_exact_cost :
    gtGlv4CtCost.totalMultiplications = 76 ∧
    gtGlv4CtCost.totalCoefficientSelections = 11520 ∧
    gtGlv4CtCost.cyclotomicSquares = 64 ∧
    gtGlv4CtCost.loopIterations = 64 := by
  decide

theorem gt_glv4_ct_fixed_selection_shape :
    gtGlv4CtCost.tableCandidatesPerIteration = 15 ∧
    gtGlv4CtCost.coefficientsPerCandidate = 12 ∧
    gtGlv4CtCost.tableEntries = 16 := by
  decide

/-- The two production unsigned executors use the same final table. -/
theorem gt_source_table_has_sixteen_entries : Fintype.card (Fin 16) = 16 := by
  decide

/-- `CommonElement<6>` is the largest delayed-reduction accumulator used. -/
def gtGlv4CommonElementBound : Nat := 6

theorem gt_glv4_common_element_bound : gtGlv4CommonElementBound = 6 := rfl

structure GtGlv4VtBounds where
  tableEntries : Nat
  tableBuildMultiplications : Nat
  loopIterations : Nat
  maximumCyclotomicSquares : Nat
  maximumLoopMultiplications : Nat
  normalizationMultiplications : Nat
  deriving DecidableEq

/--
The first nonzero joint bit initializes without a multiply.  At most the
remaining 63 iterations square and multiply; the all-zero input does less.
-/
def gtGlv4UnsignedVtBounds : GtGlv4VtBounds where
  tableEntries := 16
  tableBuildMultiplications := 11
  loopIterations := 64
  maximumCyclotomicSquares := 63
  maximumLoopMultiplications := 63
  normalizationMultiplications := 1

def GtGlv4VtBounds.maximumTotalMultiplications
    (cost : GtGlv4VtBounds) : Nat :=
  cost.tableBuildMultiplications + cost.maximumLoopMultiplications +
    cost.normalizationMultiplications

theorem gt_glv4_unsigned_vt_bounds :
    gtGlv4UnsignedVtBounds.maximumTotalMultiplications = 75 ∧
    gtGlv4UnsignedVtBounds.maximumCyclotomicSquares = 63 ∧
    gtGlv4UnsignedVtBounds.loopIterations = 64 := by
  decide

/-- Binary CT performs one square, one product, and twelve selects per bit. -/
structure GtBinaryCtCost where
  loopIterations : Nat
  cyclotomicSquares : Nat
  candidateMultiplications : Nat
  coefficientSelections : Nat
  normalizationMultiplications : Nat
  deriving DecidableEq

def gtBinaryCtCost (bits : Nat) : GtBinaryCtCost where
  loopIterations := bits
  cyclotomicSquares := bits
  candidateMultiplications := bits
  coefficientSelections := 12 * bits
  normalizationMultiplications := 1

theorem gt_binary_ct_exact_cost (bits : Nat) :
    (gtBinaryCtCost bits).loopIterations = bits ∧
    (gtBinaryCtCost bits).cyclotomicSquares = bits ∧
    (gtBinaryCtCost bits).candidateMultiplications = bits ∧
    (gtBinaryCtCost bits).coefficientSelections = 12 * bits := by
  simp [gtBinaryCtCost]

/-- Source dispatch is GLV4 exactly on the public interval 129 through 256. -/
def gtUsesGlv4 (bits : Nat) : Bool := decide (bits > 128 ∧ bits <= 256)

theorem gt_uses_glv4_iff (bits : Nat) :
    gtUsesGlv4 bits = true ↔ 128 < bits ∧ bits <= 256 := by
  simp [gtUsesGlv4]

end CatBLST.AVX.MatrixNoK
