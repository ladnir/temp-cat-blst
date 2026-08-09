/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Fp12CyclotomicConstants
import CatBLST.AVX.MatrixNoK.Fp12FrobeniusPhysical

/-!
# Physical bounds for the split cyclotomic-square schedule

`ReducedElement<N>` retains the selected M1 element bound while `N` tracks the
RNS range.  Every production call reaches this schedule with `N ≤ 6`.  A
quadratic product therefore occupies at most four M1 element units and 144 RNS
units; a product with a generated standard constant occupies at most four M1
element units and 24 RNS units.

The first split has six quadratic units plus one constant unit in every output.
The second split's offset-bearing outputs preserve the literal positive and
negative product counts.  The remaining four outputs contain two quadratic
units and one constant unit.  These ledgers are the exact source-expression
counts before `batch_reduce`; the subsequent factors three and six both fit the
declared `CommonElement<6>` result.
-/

namespace CatBLST.AVX.MatrixNoK

def cyclotomicMaximumInputRnsBound : Nat := 6
def cyclotomicElementQuadraticUnit : Nat := standardBound * standardBound
def cyclotomicElementConstantUnit : Nat := standardBound * standardBound
def cyclotomicRnsQuadraticUnit : Nat :=
  standardBound * standardBound * cyclotomicMaximumInputRnsBound *
    cyclotomicMaximumInputRnsBound
def cyclotomicRnsConstantUnit : Nat :=
  standardBound * standardBound * cyclotomicMaximumInputRnsBound

def cyclotomicFirstQuadraticUnits (_ : Fin 6) : Nat := 6
def cyclotomicFirstConstantUnits (_ : Fin 6) : Nat := 1

def cyclotomicFirstElementUpper (index : Fin 6) : Nat :=
  cyclotomicFirstQuadraticUnits index * cyclotomicElementQuadraticUnit +
    cyclotomicFirstConstantUnits index * cyclotomicElementConstantUnit

def cyclotomicFirstRnsUpper (index : Fin 6) : Nat :=
  cyclotomicFirstQuadraticUnits index * cyclotomicRnsQuadraticUnit +
    cyclotomicFirstConstantUnits index * cyclotomicRnsConstantUnit

/-- `[poly6, poly7, poly8, poly9, poly10, poly11]`. -/
def cyclotomicSecondOffsetUnits : Fin 6 → Nat := ![1, 1, 0, 0, 0, 0]
def cyclotomicSecondPositiveQuadraticUnits : Fin 6 → Nat :=
  ![2, 4, 2, 2, 2, 2]
def cyclotomicSecondNegativeQuadraticUnits : Fin 6 → Nat :=
  ![2, 0, 0, 0, 0, 0]
def cyclotomicSecondConstantUnits (_ : Fin 6) : Nat := 1

def cyclotomicSecondElementLower (index : Fin 6) : Nat :=
  cyclotomicSecondOffsetUnits index * fp12Version3OffsetElementBound -
    cyclotomicSecondNegativeQuadraticUnits index *
      cyclotomicElementQuadraticUnit

def cyclotomicSecondElementUpper (index : Fin 6) : Nat :=
  cyclotomicSecondOffsetUnits index *
      (fp12Version3OffsetElementBound + 1) +
    cyclotomicSecondPositiveQuadraticUnits index *
      cyclotomicElementQuadraticUnit +
    cyclotomicSecondConstantUnits index * cyclotomicElementConstantUnit

def cyclotomicSecondRnsLower (index : Fin 6) : Nat :=
  cyclotomicSecondOffsetUnits index * fp12Version3OffsetRnsBound -
    cyclotomicSecondNegativeQuadraticUnits index *
      cyclotomicRnsQuadraticUnit

def cyclotomicSecondRnsUpper (index : Fin 6) : Nat :=
  cyclotomicSecondOffsetUnits index *
      (fp12Version3OffsetRnsBound + 1) +
    cyclotomicSecondPositiveQuadraticUnits index *
      cyclotomicRnsQuadraticUnit +
    cyclotomicSecondConstantUnits index * cyclotomicRnsConstantUnit

theorem cyclotomic_product_unit_values :
    cyclotomicElementQuadraticUnit = 4 ∧
      cyclotomicElementConstantUnit = 4 ∧
      cyclotomicRnsQuadraticUnit = 144 ∧
      cyclotomicRnsConstantUnit = 24 := by
  norm_num [cyclotomicElementQuadraticUnit,
    cyclotomicElementConstantUnit, cyclotomicRnsQuadraticUnit,
    cyclotomicRnsConstantUnit, cyclotomicMaximumInputRnsBound,
    standardBound]

theorem cyclotomic_first_split_bounds (index : Fin 6) :
    cyclotomicFirstElementUpper index = 28 ∧
      cyclotomicFirstRnsUpper index = 888 ∧
      cyclotomicFirstElementUpper index ≤ fp12ReadyBound ∧
      cyclotomicFirstRnsUpper index < fp12MaximumCompletedRnsBound := by
  fin_cases index <;>
    norm_num [cyclotomicFirstElementUpper, cyclotomicFirstRnsUpper,
      cyclotomicFirstQuadraticUnits, cyclotomicFirstConstantUnits,
      cyclotomicElementQuadraticUnit, cyclotomicElementConstantUnit,
      cyclotomicRnsQuadraticUnit, cyclotomicRnsConstantUnit,
      cyclotomicMaximumInputRnsBound, standardBound, fp12ReadyBound,
      fp12MaximumCompletedRnsBound]

theorem cyclotomic_second_split_exact_bounds :
    (∀ index : Fin 6,
      cyclotomicSecondElementLower index = ![72, 80, 0, 0, 0, 0] index ∧
      cyclotomicSecondElementUpper index = ![93, 101, 12, 12, 12, 12] index ∧
      cyclotomicSecondRnsLower index = ![1152, 1440, 0, 0, 0, 0] index ∧
      cyclotomicSecondRnsUpper index = ![1753, 2041, 312, 312, 312, 312] index) := by
  intro index
  fin_cases index <;>
    norm_num [cyclotomicSecondElementLower,
      cyclotomicSecondElementUpper, cyclotomicSecondRnsLower,
      cyclotomicSecondRnsUpper, cyclotomicSecondOffsetUnits,
      cyclotomicSecondPositiveQuadraticUnits,
      cyclotomicSecondNegativeQuadraticUnits,
      cyclotomicSecondConstantUnits, cyclotomicElementQuadraticUnit,
      cyclotomicElementConstantUnit, cyclotomicRnsQuadraticUnit,
      cyclotomicRnsConstantUnit, cyclotomicMaximumInputRnsBound,
      fp12Version3OffsetElementBound, fp12Version3OffsetRnsBound,
      standardBound, multiplicationProductBound]

theorem cyclotomic_second_split_fits_reduction (index : Fin 6) :
    cyclotomicSecondElementUpper index ≤ fp12ReadyBound ∧
      cyclotomicSecondRnsUpper index < fp12MaximumCompletedRnsBound := by
  fin_cases index <;>
    norm_num [cyclotomicSecondElementUpper, cyclotomicSecondRnsUpper,
      cyclotomicSecondOffsetUnits,
      cyclotomicSecondPositiveQuadraticUnits,
      cyclotomicSecondConstantUnits, cyclotomicElementQuadraticUnit,
      cyclotomicElementConstantUnit, cyclotomicRnsQuadraticUnit,
      cyclotomicRnsConstantUnit, cyclotomicMaximumInputRnsBound,
      fp12Version3OffsetElementBound, fp12Version3OffsetRnsBound,
      fp12ReadyBound, fp12MaximumCompletedRnsBound, standardBound,
      multiplicationProductBound]

theorem cyclotomic_second_split_subtraction_safe :
    cyclotomicSecondNegativeQuadraticUnits 0 *
        cyclotomicElementQuadraticUnit ≤
      fp12Version3OffsetElementBound ∧
      cyclotomicSecondNegativeQuadraticUnits 0 *
          cyclotomicRnsQuadraticUnit ≤
        fp12Version3OffsetRnsBound := by
  norm_num [cyclotomicSecondNegativeQuadraticUnits,
    cyclotomicElementQuadraticUnit, cyclotomicRnsQuadraticUnit,
    cyclotomicMaximumInputRnsBound, fp12Version3OffsetElementBound,
    fp12Version3OffsetRnsBound, standardBound,
    multiplicationProductBound]

theorem cyclotomic_post_reduction_scales_fit_common_six :
    (3 : Nat) ≤ cyclotomicMaximumInputRnsBound ∧
      (6 : Nat) ≤ cyclotomicMaximumInputRnsBound := by
  norm_num [cyclotomicMaximumInputRnsBound]

theorem cyclotomic_semantic_offset_is_zero :
    (fp12SemanticOffset : SelectedBaseField) = 0 := by
  simp [fp12SemanticOffset]

end CatBLST.AVX.MatrixNoK
