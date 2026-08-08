/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.SignedPipeline
import Mathlib.Tactic

/-!
# Canonical radix-2^50 import

The convert-to matrix reconstructs eight little-endian radix-2^50 digits
modulo the BLS12-381 field modulus and emits M2 residues.  Its fixed-point CRT
quotient may be exact or one high, so the common mathematical lift is signed
and lies in `[-p, p)`.  The AVX lanes themselves remain unsigned.
-/

namespace CatBLST.AVX.MatrixNoK

set_option exponentiation.threshold 512
open scoped BigOperators

def canonicalDigitValue (digits : Fin 8 → Nat) : Nat :=
  ∑ index, digits index * 2 ^ (50 * index.val)

def convertToShiftedFactor (index : Fin 8) : Nat :=
  Generated.convert_to_premult * 2 ^ (50 * index.val) %
    Generated.target_modulus

def convertToCorrectedRepresentative (digits : Fin 8 → Nat) : Int :=
  correctedCrtRepresentative digits convertToShiftedFactor
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
      ceilDiv (convertToShiftedFactor index * uint64Modulus)
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
      estimatedCrtQuotient digits convertToShiftedFactor
        Generated.target_modulus := by
  unfold generatedConvertToEstimatedQuotient estimatedCrtQuotient
  apply congrArg (fun numerator => numerator / uint64Modulus)
  unfold approximateCrtNumerator
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
      fixedMatrixLaneValue digits convertToShiftedFactor
        Generated.target_modulus Generated.convert_to_postmult
        (listAt8 Generated.moduli2 column.val) := by
  unfold generatedConvertToLaneValue fixedMatrixLaneValue
  rw [generated_convert_to_estimated_quotient_eq]
  rw [generated_convert_to_correction_formula]
  congr 1
  apply Finset.sum_congr rfl
  intro row _
  rw [generated_convert_to_matrix_formula]

theorem convert_to_corrected_range
    (digits : Fin 8 → Nat)
    (digitsFit : ∀ index, digits index < 2 ^ 50) :
    -(Generated.target_modulus : Int) ≤
        convertToCorrectedRepresentative digits ∧
      convertToCorrectedRepresentative digits < Generated.target_modulus := by
  unfold convertToCorrectedRepresentative
  apply corrected_representative_range digits convertToShiftedFactor
    Generated.target_modulus (by norm_num [Generated.target_modulus])
  intro index
  exact (digitsFit index).trans (by norm_num [ifmaRadix])

/-- A negative canonical-import lift can only be the small fixed-point
overshoot.  With eight 50-bit input digits and 64 fractional quotient bits,
its magnitude is strictly below `p / 2^11`. -/
theorem convert_to_negative_magnitude_bound
    (digits : Fin 8 → Nat) (magnitude : Nat)
    (digitsFit : ∀ index, digits index < 2 ^ 50)
    (magnitudePositive : 0 < magnitude)
    (negativeLift :
      convertToCorrectedRepresentative digits = -(magnitude : Int)) :
    2048 * magnitude < Generated.target_modulus := by
  let exact := exactCrtNumerator digits convertToShiftedFactor
  let approximate := approximateCrtNumerator digits convertToShiftedFactor
    Generated.target_modulus
  let exactQuotient := exactCrtQuotient digits convertToShiftedFactor
    Generated.target_modulus
  let estimate := estimatedCrtQuotient digits convertToShiftedFactor
    Generated.target_modulus
  let remainder := exact % Generated.target_modulus
  have targetPositive : 0 < Generated.target_modulus := by
    norm_num [Generated.target_modulus]
  have scalarFits : ∀ index, digits index < ifmaRadix := by
    intro index
    exact (digitsFit index).trans (by norm_num [ifmaRadix])
  have bracket := estimated_quotient_bracket digits convertToShiftedFactor
    Generated.target_modulus targetPositive scalarFits
  have exactSplit := Nat.mod_add_div exact Generated.target_modulus
  have exactDecomposition :
      exact = exactQuotient * Generated.target_modulus + remainder := by
    calc
      exact = exact % Generated.target_modulus +
          Generated.target_modulus * (exact / Generated.target_modulus) :=
        exactSplit.symm
      _ = (exact / Generated.target_modulus) *
          Generated.target_modulus + exact % Generated.target_modulus := by ring
      _ = exactQuotient * Generated.target_modulus + remainder := by rfl
  have remainderBound : remainder < Generated.target_modulus := by
    exact Nat.mod_lt exact targetPositive
  have exactSplitInt :
      (exact : Int) = (remainder : Int) +
        (Generated.target_modulus : Int) * (exactQuotient : Int) := by
    exact_mod_cast exactSplit.symm
  have estimateIsOneHigh : estimate = exactQuotient + 1 := by
    have estimateNotExact : estimate ≠ exactQuotient := by
      intro estimateExact
      have liftNonnegative :
          (0 : Int) ≤ convertToCorrectedRepresentative digits := by
        unfold convertToCorrectedRepresentative correctedCrtRepresentative
        change (0 : Int) ≤ (exact : Int) - (estimate : Int) *
          Generated.target_modulus
        rw [estimateExact, exactSplitInt]
        ring_nf
        exact_mod_cast Nat.zero_le remainder
      rw [negativeLift] at liftNonnegative
      omega
    dsimp [estimate, exactQuotient] at bracket ⊢
    omega
  have magnitudeEq : magnitude = Generated.target_modulus - remainder := by
    have liftFormula :
        convertToCorrectedRepresentative digits =
          (remainder : Int) - Generated.target_modulus := by
      unfold convertToCorrectedRepresentative correctedCrtRepresentative
      change (exact : Int) - (estimate : Int) *
        Generated.target_modulus = _
      rw [estimateIsOneHigh, exactSplitInt]
      push_cast
      ring
    rw [negativeLift] at liftFormula
    have castEquality :
        (magnitude : Int) =
          (Generated.target_modulus - remainder : Nat) := by
      rw [Nat.cast_sub remainderBound.le]
      omega
    exact_mod_cast castEquality
  have errorIdentity := aggregate_ceiling_error_identity
    digits convertToShiftedFactor Generated.target_modulus targetPositive
  have quotientFloorLower := Nat.mul_div_le approximate uint64Modulus
  have approximateLower :
      (exactQuotient + 1) * uint64Modulus ≤ approximate := by
    change approximate / uint64Modulus = exactQuotient + 1 at estimateIsOneHigh
    rw [← estimateIsOneHigh]
    simpa [Nat.mul_comm] using quotientFloorLower
  have errorLower :
      (Generated.target_modulus - remainder) * uint64Modulus ≤
        ∑ index, digits index *
          fixedPointCeilingError (convertToShiftedFactor index)
            Generated.target_modulus := by
    change exact * uint64Modulus +
        ∑ index, digits index *
          fixedPointCeilingError (convertToShiftedFactor index)
            Generated.target_modulus =
      approximate * Generated.target_modulus at errorIdentity
    rw [exactDecomposition] at errorIdentity
    have remainderSplit :
        remainder + (Generated.target_modulus - remainder) =
          Generated.target_modulus :=
      Nat.add_sub_of_le remainderBound.le
    have scaledApproximateLower := Nat.mul_le_mul_right
      Generated.target_modulus approximateLower
    have leftDecomposition :
        (exactQuotient + 1) * uint64Modulus *
            Generated.target_modulus =
          (exactQuotient * Generated.target_modulus + remainder) *
              uint64Modulus +
            (Generated.target_modulus - remainder) * uint64Modulus := by
      calc
        _ = exactQuotient * Generated.target_modulus * uint64Modulus +
            Generated.target_modulus * uint64Modulus := by ring
        _ = exactQuotient * Generated.target_modulus * uint64Modulus +
            (remainder + (Generated.target_modulus - remainder)) *
              uint64Modulus := by rw [remainderSplit]
        _ = _ := by ring
    rw [leftDecomposition, ← errorIdentity] at scaledApproximateLower
    exact Nat.add_le_add_iff_left.mp scaledApproximateLower
  have enclosure := aggregate_quotient_enclosure
    digits convertToShiftedFactor Generated.target_modulus targetPositive
  have errorUpper :
      ∑ index, digits index *
          fixedPointCeilingError (convertToShiftedFactor index)
            Generated.target_modulus ≤
        (Generated.target_modulus - 1) * ∑ index, digits index := by
    have upper := enclosure.2
    rw [← errorIdentity] at upper
    omega
  have digitSumBound : (∑ index, digits index) < 8 * 2 ^ 50 :=
    fin8_sum_lt digitsFit
  have errorEnvelope :
      (Generated.target_modulus - 1) * ∑ index, digits index <
        Generated.target_modulus * 2 ^ 53 := by
    calc
      _ ≤ Generated.target_modulus * ∑ index, digits index :=
        Nat.mul_le_mul_right _ (Nat.sub_le _ _)
      _ < Generated.target_modulus * (8 * 2 ^ 50) :=
        (Nat.mul_lt_mul_left targetPositive).2 digitSumBound
      _ = Generated.target_modulus * 2 ^ 53 := by ring
  rw [magnitudeEq]
  norm_num [uint64Modulus] at errorLower ⊢
  nlinarith [errorLower, errorUpper.trans_lt errorEnvelope]

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
  unfold convertToCorrectedRepresentative correctedCrtRepresentative
  push_cast
  rw [convert_to_exact_numerator_represents_digits]
  simp

theorem generated_convert_to_lane_represents_signed_import
    (digits : Fin 8 → Nat) (column : Fin 8) :
    generatedConvertToLaneValue digits column =
      (Generated.convert_to_postmult :
        ZMod (listAt8 Generated.moduli2 column.val)) *
        convertToCorrectedRepresentative digits := by
  rw [generated_convert_to_lane_eq_fixed]
  exact fixed_matrix_lane_represents_corrected_crt _ _ _ _ _
    (by fin_cases column <;> norm_num [Generated.moduli2, listAt8])

end CatBLST.AVX.MatrixNoK
