/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.SignedPipeline
import Mathlib.Tactic

/-!
# Canonical radix-2^50 import

The ConvertTo matrix reconstructs eight little-endian radix-2^50 digits
modulo the BLS12-381 field modulus and emits M2 residues. Unlike the inherited
VROOM table, its fixed-point CRT coefficients round downward. The aggregate
estimate is therefore exact or one low, and every canonical import selects a
nonnegative representative in `[0, 2p)`.
-/

namespace CatBLST.AVX.MatrixNoK

set_option exponentiation.threshold 512
set_option maxRecDepth 10000
open scoped BigOperators

def canonicalDigitValue (digits : Fin 8 → Nat) : Nat :=
  ∑ index, digits index * 2 ^ (50 * index.val)

def convertToShiftedFactor (index : Fin 8) : Nat :=
  Generated.convert_to_premult * 2 ^ (50 * index.val) %
    Generated.target_modulus

def convertToCorrectedRepresentative (digits : Fin 8 → Nat) : Nat :=
  floorCorrectedCrtRepresentative digits convertToShiftedFactor
    Generated.target_modulus

set_option maxRecDepth 10000 in
theorem convert_to_factor_formula (index : Fin 8) :
    listAt8
        (shiftedFactors [Generated.target_modulus] [1]
          Generated.convert_to_premult 8 50)
        index.val = convertToShiftedFactor index := by
  have targetPositive : 0 < Generated.target_modulus := by
    norm_num [Generated.target_modulus]
  have targetDiv :
      Generated.target_modulus / Generated.target_modulus = 1 :=
    Nat.div_self targetPositive
  fin_cases index <;>
    simp [shiftedFactors, crtFactors, convertToShiftedFactor, listAt8,
      targetDiv]

theorem generated_convert_to_quotient_factor_formula (index : Fin 8) :
    listAt8 Generated.convert_to_shifted_quotient_estimations index.val =
      (convertToShiftedFactor index * uint64Modulus) /
        Generated.target_modulus := by
  fin_cases index <;> decide

theorem generated_convert_to_matrix_formula (row column : Fin 8) :
    matrixEntry Generated.convert_to_rns_mat row.val column.val =
      fixedMatrixEntry (convertToShiftedFactor row)
        Generated.convert_to_postmult
        (listAt8 Generated.moduli2 column.val) := by
  fin_cases row <;> fin_cases column <;> decide

theorem generated_convert_to_correction_formula (column : Fin 8) :
    listAt8 Generated.convert_to_correction column.val =
      fixedCorrection Generated.target_modulus
        Generated.convert_to_postmult
        (listAt8 Generated.moduli2 column.val) := by
  fin_cases column <;> decide

def generatedConvertToEstimatedQuotient (digits : Fin 8 → Nat) : Nat :=
  (∑ row, digits row *
      listAt8 Generated.convert_to_shifted_quotient_estimations row.val) /
    uint64Modulus

theorem generated_convert_to_estimated_quotient_eq
    (digits : Fin 8 → Nat) :
    generatedConvertToEstimatedQuotient digits =
      floorEstimatedCrtQuotient digits convertToShiftedFactor
        Generated.target_modulus := by
  unfold generatedConvertToEstimatedQuotient floorEstimatedCrtQuotient
  apply congrArg (fun numerator => numerator / uint64Modulus)
  unfold approximateCrtNumeratorFloor
  apply Finset.sum_congr rfl
  intro row _
  rw [generated_convert_to_quotient_factor_formula]

def generatedConvertToLaneValue
    (digits : Fin 8 → Nat) (column : Fin 8) :
    ZMod (listAt8 Generated.moduli2 column.val) :=
  (∑ row, (digits row : ZMod (listAt8 Generated.moduli2 column.val)) *
      matrixEntry Generated.convert_to_rns_mat row.val column.val) +
    (generatedConvertToEstimatedQuotient digits :
      ZMod (listAt8 Generated.moduli2 column.val)) *
      listAt8 Generated.convert_to_correction column.val

theorem generated_convert_to_lane_eq_fixed
    (digits : Fin 8 → Nat) (column : Fin 8) :
    generatedConvertToLaneValue digits column =
      floorFixedMatrixLaneValue digits convertToShiftedFactor
        Generated.target_modulus Generated.convert_to_postmult
        (listAt8 Generated.moduli2 column.val) := by
  unfold generatedConvertToLaneValue floorFixedMatrixLaneValue
  rw [generated_convert_to_estimated_quotient_eq]
  rw [generated_convert_to_correction_formula]
  congr 1
  apply Finset.sum_congr rfl
  intro row _
  rw [generated_convert_to_matrix_formula]

theorem convert_to_digit_sum_bound
    (digits : Fin 8 → Nat)
    (digitsFit : ∀ index, digits index < 2 ^ 50) :
    (∑ index, digits index) < uint64Modulus := by
  have sumBound : (∑ index, digits index) < 8 * 2 ^ 50 :=
    fin8_sum_lt digitsFit
  exact sumBound.trans (by norm_num [uint64Modulus])

theorem convert_to_corrected_range
    (digits : Fin 8 → Nat)
    (digitsFit : ∀ index, digits index < 2 ^ 50) :
    convertToCorrectedRepresentative digits <
      2 * Generated.target_modulus := by
  unfold convertToCorrectedRepresentative
  apply floor_corrected_representative_range
  · norm_num [Generated.target_modulus]
  · exact convert_to_digit_sum_bound digits digitsFit

theorem convert_to_exact_numerator_represents_digits
    (digits : Fin 8 → Nat) :
    (exactCrtNumerator digits convertToShiftedFactor :
        ZMod Generated.target_modulus) =
      Generated.convert_to_premult * canonicalDigitValue digits := by
  unfold exactCrtNumerator convertToShiftedFactor canonicalDigitValue
  push_cast
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro index _
  rw [ZMod.natCast_mod]
  push_cast
  ring

theorem convert_to_corrected_represents_digits
    (digits : Fin 8 → Nat) :
    (convertToCorrectedRepresentative digits :
        ZMod Generated.target_modulus) =
      Generated.convert_to_premult * canonicalDigitValue digits := by
  have targetPositive : 0 < Generated.target_modulus := by
    norm_num [Generated.target_modulus]
  have quotientBound := floor_estimated_quotient_le_exact
    digits convertToShiftedFactor Generated.target_modulus targetPositive
  have productBound := Nat.mul_le_mul_right
    Generated.target_modulus quotientBound
  have correctionBelow :
      floorEstimatedCrtQuotient digits convertToShiftedFactor
          Generated.target_modulus * Generated.target_modulus ≤
        exactCrtNumerator digits convertToShiftedFactor := by
    exact productBound.trans (by
      simpa [exactCrtQuotient, Nat.mul_comm] using
        Nat.div_mul_le_self
          (exactCrtNumerator digits convertToShiftedFactor)
          Generated.target_modulus)
  unfold convertToCorrectedRepresentative floorCorrectedCrtRepresentative
  rw [Nat.cast_sub correctionBelow]
  push_cast
  rw [convert_to_exact_numerator_represents_digits]
  simp

theorem generated_convert_to_lane_represents_import
    (digits : Fin 8 → Nat) (column : Fin 8) :
    generatedConvertToLaneValue digits column =
      (Generated.convert_to_postmult :
        ZMod (listAt8 Generated.moduli2 column.val)) *
        convertToCorrectedRepresentative digits := by
  rw [generated_convert_to_lane_eq_fixed]
  exact floor_fixed_matrix_lane_represents_corrected_crt _ _ _ _ _
    (by fin_cases column <;> norm_num [Generated.moduli2, listAt8])

/- The repaired canonical import composes with the existing nonnegative R2
contract once the per-lane Montgomery reducer has established its ordinary
residue postcondition. -/
theorem canonical_import_r2_expansion_lane
    (digits : Fin 8 → Nat)
    (digitsFit : ∀ index, digits index < 2 ^ 50)
    (m2 : Fin 8 → Nat)
    (m2Fit : ∀ index,
      m2 index < 2 * listAt8 Generated.moduli2 index.val)
    (residuesRepresent : ∀ index,
      Nat.ModEq (listAt8 Generated.moduli2 index.val)
        (m2 index)
        (r2InputValue (convertToCorrectedRepresentative digits)))
    (column : Fin 8) :
    generatedR2LaneValue m2 column =
      (Generated.r2_postmult :
        ZMod (listAt8 Generated.moduli1 column.val)) *
        convertToCorrectedRepresentative digits := by
  let representative := convertToCorrectedRepresentative digits
  have representativeBound : representative <
      18 * Generated.target_modulus := by
    have repairedBound := convert_to_corrected_range digits digitsFit
    dsimp [representative]
    nlinarith
  exact generated_r2_standard_lane_represents_input representative m2
    representativeBound m2Fit residuesRepresent column

end CatBLST.AVX.MatrixNoK
