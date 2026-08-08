/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.CanonicalImport
import Mathlib.Tactic

/-!
# Canonical radix-2^52 export

`ConvertFromMatrix::reduce_raw` emits eight wide radix-2^52 coefficients.
`from_vec` recombines each low/high pair and then recombines the eight lanes.
The definitions below model that recombined integer directly.  This keeps the
proof aligned with the generated table while avoiding an artificial limb
normalization requirement: carries between the wide coefficients are handled
by ordinary integer recombination.

ConvertTo stores an imported field value with one factor of `M1`.  MatrixNoK
multiplication preserves that representation, and ConvertFrom removes it with
`inv_factor_target`.  The final theorems close both the canonical round trip
and the selected import/multiply/export composition.
-/

namespace CatBLST.AVX.MatrixNoK

set_option exponentiation.threshold 512
set_option maxRecDepth 10000
open scoped BigOperators

def generatedConvertFromRowValue (row : Fin 8) : Nat :=
  ∑ column : Fin 8,
    matrixEntry Generated.convert_from_rns_mat row.val column.val *
      montgomeryRadix ^ column.val

def generatedConvertFromCorrectionValue : Nat :=
  ∑ column : Fin 8,
    listAt8 Generated.convert_from_correction column.val *
      montgomeryRadix ^ column.val

theorem generated_convert_from_row_formula (row : Fin 8) :
    generatedConvertFromRowValue row =
      fixedMatrixEntry (r2ShiftedFactor row)
        Generated.convert_from_postmult Generated.target_modulus := by
  fin_cases row <;> decide

theorem generated_convert_from_correction_formula :
    generatedConvertFromCorrectionValue =
      fixedCorrection Generated.M2 Generated.convert_from_postmult
        Generated.target_modulus := by
  decide

theorem generated_convert_from_quotient_factor_formula (row : Fin 8) :
    listAt8 Generated.convert_from_shifted_quotient_estimations row.val =
      ceilDiv (r2ShiftedFactor row * uint64Modulus) Generated.M2 := by
  fin_cases row <;> decide

def generatedConvertFromEstimatedQuotient
    (residues : Fin 8 → Nat) : Nat :=
  (∑ row, residues row *
      listAt8 Generated.convert_from_shifted_quotient_estimations row.val) /
    uint64Modulus

theorem generated_convert_from_estimated_quotient_eq
    (residues : Fin 8 → Nat) :
    generatedConvertFromEstimatedQuotient residues =
      estimatedCrtQuotient residues r2ShiftedFactor Generated.M2 := by
  unfold generatedConvertFromEstimatedQuotient estimatedCrtQuotient
  apply congrArg (fun numerator => numerator / uint64Modulus)
  unfold approximateCrtNumerator
  apply Finset.sum_congr rfl
  intro row _
  rw [generated_convert_from_quotient_factor_formula row]

def generatedConvertFromValue (residues : Fin 8 → Nat) : Nat :=
  (∑ row, residues row * generatedConvertFromRowValue row) +
    generatedConvertFromEstimatedQuotient residues *
      generatedConvertFromCorrectionValue

def generatedConvertFromWideDigit
    (residues : Fin 8 → Nat) (column : Fin 8) : Nat :=
  (∑ row, residues row *
      matrixEntry Generated.convert_from_rns_mat row.val column.val) +
    generatedConvertFromEstimatedQuotient residues *
      listAt8 Generated.convert_from_correction column.val

def generatedConvertFromRecombinedValue
    (residues : Fin 8 → Nat) : Nat :=
  ∑ column, generatedConvertFromWideDigit residues column *
    montgomeryRadix ^ column.val

/-- Recombining the eight wide output coefficients is exactly the expanded
integer model used below. This is the source-level algebra of `reduce_raw`
followed by `from_vec`; no per-coefficient carry normalization is required. -/
theorem generated_convert_from_recombined_eq
    (residues : Fin 8 → Nat) :
    generatedConvertFromRecombinedValue residues =
      generatedConvertFromValue residues := by
  unfold generatedConvertFromRecombinedValue generatedConvertFromWideDigit
    generatedConvertFromValue generatedConvertFromRowValue
    generatedConvertFromCorrectionValue
  simp_rw [add_mul]
  rw [Finset.sum_add_distrib]
  congr 1
  · simp_rw [Finset.sum_mul]
    rw [Finset.sum_comm]
    apply Finset.sum_congr rfl
    intro row _
    rw [Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro column _
    ring
  · rw [Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro column _
    ring

theorem generated_convert_from_matrix_entry_lt_radix
    (row column : Fin 8) :
    matrixEntry Generated.convert_from_rns_mat row.val column.val <
      ifmaRadix := by
  fin_cases row <;> fin_cases column <;> decide

theorem generated_convert_from_correction_digit_lt_radix
    (column : Fin 8) :
    listAt8 Generated.convert_from_correction column.val < ifmaRadix := by
  fin_cases column <;> decide

theorem generated_convert_from_estimate_lt_word (row : Fin 8) :
    listAt8 Generated.convert_from_shifted_quotient_estimations row.val <
      uint64Modulus := by
  fin_cases row <;> decide

theorem generated_convert_from_quotient_words_fit
    (residues : Fin 8 → Nat)
    (residuesFit : ∀ index,
      residues index < 2 * listAt8 Generated.moduli2 index.val) :
    let quotientNumerator := ∑ index, residues index *
      listAt8 Generated.convert_from_shifted_quotient_estimations index.val
    quotientNumerator < uint128Modulus ∧
      generatedConvertFromEstimatedQuotient residues < 2 ^ 55 ∧
      generatedConvertFromEstimatedQuotient residues / ifmaRadix < 8 := by
  dsimp only
  have scalarFits := r2_standard_residue_lt_ifma residues residuesFit
  have words := quotient_words_fit residues
    (fun index =>
      listAt8 Generated.convert_from_shifted_quotient_estimations index.val)
    scalarFits generated_convert_from_estimate_lt_word
  simpa [quotientNumerator, generatedConvertFromEstimatedQuotient] using words

def generatedConvertFromFinalLow
    (residues : Fin 8 → Nat) (column : Fin 8) : Nat :=
  lowMatrixAccumulator residues Generated.convert_from_rns_mat column +
    lowDigitProduct
      (listAt8 Generated.convert_from_correction column.val)
      (generatedConvertFromEstimatedQuotient residues)

def generatedConvertFromFinalHigh
    (residues : Fin 8 → Nat) (column : Fin 8) : Nat :=
  highMatrixAccumulator residues Generated.convert_from_rns_mat column +
    highDigitProduct
      (listAt8 Generated.convert_from_correction column.val)
      (generatedConvertFromEstimatedQuotient residues) +
    listAt8 Generated.convert_from_correction column.val *
      (generatedConvertFromEstimatedQuotient residues / ifmaRadix)

theorem generated_convert_from_accumulators_fit_word
    (residues : Fin 8 → Nat)
    (residuesFit : ∀ index,
      residues index < 2 * listAt8 Generated.moduli2 index.val)
    (column : Fin 8) :
    generatedConvertFromFinalLow residues column < uint64Modulus ∧
      generatedConvertFromFinalHigh residues column < uint64Modulus := by
  have lowFits : lowMatrixAccumulator residues
      Generated.convert_from_rns_mat column < 8 * ifmaRadix := by
    apply fin8_sum_lt
    intro row
    exact (digit_products_fit (residues row)
      (matrixEntry Generated.convert_from_rns_mat row.val column.val)).1
  have highFits : highMatrixAccumulator residues
      Generated.convert_from_rns_mat column < 8 * ifmaRadix := by
    apply fin8_sum_lt
    intro row
    exact (digit_products_fit (residues row)
      (matrixEntry Generated.convert_from_rns_mat row.val column.val)).2
  have quotientFits :=
    (generated_convert_from_quotient_words_fit residues residuesFit).2.1
  have correctionFits :=
    generated_convert_from_correction_digit_lt_radix column
  simpa [generatedConvertFromFinalLow, generatedConvertFromFinalHigh] using
    fixed_matrix_correction_fits_word
      (highMatrixAccumulator residues Generated.convert_from_rns_mat column)
      (lowMatrixAccumulator residues Generated.convert_from_rns_mat column)
      (listAt8 Generated.convert_from_correction column.val)
      (generatedConvertFromEstimatedQuotient residues)
      highFits lowFits correctionFits quotientFits

theorem digit_product_reconstruction (left right : Nat) :
    lowDigitProduct left right + ifmaRadix * highDigitProduct left right =
      (left % ifmaRadix) * (right % ifmaRadix) := by
  simpa [lowDigitProduct, highDigitProduct, Nat.mul_comm] using
    Nat.mod_add_div
      ((left % ifmaRadix) * (right % ifmaRadix)) ifmaRadix

theorem generated_convert_from_accumulators_reconstruct_wide_digit
    (residues : Fin 8 → Nat)
    (residuesFit : ∀ index,
      residues index < 2 * listAt8 Generated.moduli2 index.val)
    (column : Fin 8) :
    generatedConvertFromFinalLow residues column + ifmaRadix *
        generatedConvertFromFinalHigh residues column =
      generatedConvertFromWideDigit residues column := by
  have scalarFits := r2_standard_residue_lt_ifma residues residuesFit
  have matrixReconstruction :
      lowMatrixAccumulator residues Generated.convert_from_rns_mat column +
          ifmaRadix * highMatrixAccumulator residues
            Generated.convert_from_rns_mat column =
        ∑ row, residues row *
          matrixEntry Generated.convert_from_rns_mat row.val column.val := by
    unfold lowMatrixAccumulator highMatrixAccumulator
    rw [Finset.mul_sum, ← Finset.sum_add_distrib]
    apply Finset.sum_congr rfl
    intro row _
    rw [digit_product_reconstruction]
    rw [Nat.mod_eq_of_lt (scalarFits row)]
    rw [Nat.mod_eq_of_lt
      (generated_convert_from_matrix_entry_lt_radix row column)]
  have correctionReconstruction :
      lowDigitProduct
          (listAt8 Generated.convert_from_correction column.val)
          (generatedConvertFromEstimatedQuotient residues) +
        ifmaRadix *
          (highDigitProduct
              (listAt8 Generated.convert_from_correction column.val)
              (generatedConvertFromEstimatedQuotient residues) +
            listAt8 Generated.convert_from_correction column.val *
              (generatedConvertFromEstimatedQuotient residues / ifmaRadix)) =
        listAt8 Generated.convert_from_correction column.val *
          generatedConvertFromEstimatedQuotient residues := by
    let correction := listAt8 Generated.convert_from_correction column.val
    let quotient := generatedConvertFromEstimatedQuotient residues
    have correctionFits :=
      generated_convert_from_correction_digit_lt_radix column
    have productSplit := digit_product_reconstruction correction quotient
    have quotientSplit := Nat.mod_add_div quotient ifmaRadix
    dsimp [correction, quotient] at productSplit quotientSplit ⊢
    calc
      _ = (lowDigitProduct
            (listAt8 Generated.convert_from_correction column.val)
            (generatedConvertFromEstimatedQuotient residues) +
          ifmaRadix * highDigitProduct
            (listAt8 Generated.convert_from_correction column.val)
            (generatedConvertFromEstimatedQuotient residues)) +
          listAt8 Generated.convert_from_correction column.val *
            (ifmaRadix *
              (generatedConvertFromEstimatedQuotient residues / ifmaRadix)) := by
        ring
      _ = listAt8 Generated.convert_from_correction column.val *
          (generatedConvertFromEstimatedQuotient residues % ifmaRadix) +
          listAt8 Generated.convert_from_correction column.val *
            (ifmaRadix *
              (generatedConvertFromEstimatedQuotient residues / ifmaRadix)) := by
        rw [productSplit, Nat.mod_eq_of_lt correctionFits]
      _ = listAt8 Generated.convert_from_correction column.val *
          (generatedConvertFromEstimatedQuotient residues % ifmaRadix +
            ifmaRadix *
              (generatedConvertFromEstimatedQuotient residues / ifmaRadix)) := by
        ring
      _ = listAt8 Generated.convert_from_correction column.val *
          generatedConvertFromEstimatedQuotient residues := by
        rw [quotientSplit]
  unfold generatedConvertFromFinalLow generatedConvertFromFinalHigh
    generatedConvertFromWideDigit
  calc
    _ = (lowMatrixAccumulator residues Generated.convert_from_rns_mat column +
          ifmaRadix * highMatrixAccumulator residues
            Generated.convert_from_rns_mat column) +
        (lowDigitProduct
            (listAt8 Generated.convert_from_correction column.val)
            (generatedConvertFromEstimatedQuotient residues) +
          ifmaRadix *
            (highDigitProduct
                (listAt8 Generated.convert_from_correction column.val)
                (generatedConvertFromEstimatedQuotient residues) +
              listAt8 Generated.convert_from_correction column.val *
                (generatedConvertFromEstimatedQuotient residues / ifmaRadix))) := by
      ring
    _ = (∑ row, residues row *
          matrixEntry Generated.convert_from_rns_mat row.val column.val) +
        listAt8 Generated.convert_from_correction column.val *
          generatedConvertFromEstimatedQuotient residues := by
      rw [matrixReconstruction, correctionReconstruction]
    _ = _ := by ring

theorem generated_convert_from_value_eq_fixed
    (residues : Fin 8 → Nat) :
    (generatedConvertFromValue residues : ZMod Generated.target_modulus) =
      fixedMatrixLaneValue residues r2ShiftedFactor Generated.M2
        Generated.convert_from_postmult Generated.target_modulus := by
  unfold generatedConvertFromValue fixedMatrixLaneValue
  rw [generated_convert_from_estimated_quotient_eq]
  rw [generated_convert_from_correction_formula]
  simp only [Nat.cast_add, Nat.cast_mul, Nat.cast_sum]
  simp_rw [generated_convert_from_row_formula]

theorem generated_convert_from_represents_representative
    (representative : Nat) (residues : Fin 8 → Nat)
    (representativeBound : representative < 18 * Generated.target_modulus)
    (residuesFit : ∀ index,
      residues index < 2 * listAt8 Generated.moduli2 index.val)
    (residuesRepresent : ∀ index,
      Nat.ModEq (listAt8 Generated.moduli2 index.val)
        (residues index) (r2InputValue representative)) :
    (generatedConvertFromValue residues : ZMod Generated.target_modulus) =
      Generated.convert_from_postmult * representative := by
  rw [generated_convert_from_value_eq_fixed]
  rw [fixed_matrix_lane_represents_corrected_crt]
  · rw [r2_standard_corrected_representative_is_exact representative residues
      representativeBound residuesFit residuesRepresent]
    norm_cast
  · norm_num [Generated.target_modulus]

theorem target_inverse_factor_zmod :
    (Generated.M1 : ZMod Generated.target_modulus) *
        Generated.convert_from_postmult = 1 := by
  have relation := target_inverse_factor
  have castRelation := congrArg
    (fun value : Nat => (value : ZMod Generated.target_modulus)) relation
  simpa [Generated.convert_from_postmult, Generated.inv_factor_target] using
    castRelation

theorem export_scaling_composes
    (leftInternal rightInternal productRepresentative leftValue rightValue :
      ZMod Generated.target_modulus)
    (leftRepresents : leftInternal = Generated.M1 * leftValue)
    (rightRepresents : rightInternal = Generated.M1 * rightValue)
    (productRepresents : productRepresentative * Generated.M1 =
      leftInternal * rightInternal) :
    Generated.convert_from_postmult * productRepresentative =
      leftValue * rightValue := by
  calc
    Generated.convert_from_postmult * productRepresentative =
        (Generated.M1 * Generated.convert_from_postmult) *
          (Generated.convert_from_postmult * productRepresentative) := by
      rw [target_inverse_factor_zmod]
      simp
    _ = (Generated.convert_from_postmult *
          Generated.convert_from_postmult) *
        (productRepresentative * Generated.M1) := by ring
    _ = (Generated.convert_from_postmult *
          Generated.convert_from_postmult) *
        (leftInternal * rightInternal) := by rw [productRepresents]
    _ = (Generated.convert_from_postmult *
          Generated.convert_from_postmult) *
        ((Generated.M1 * leftValue) *
          (Generated.M1 * rightValue)) := by
      rw [leftRepresents, rightRepresents]
    _ = (Generated.M1 * Generated.convert_from_postmult) *
        (Generated.M1 * Generated.convert_from_postmult) *
          (leftValue * rightValue) := by ring
    _ = leftValue * rightValue := by
      rw [target_inverse_factor_zmod]
      simp

theorem canonical_import_export_round_trip
    (digits : Fin 8 → Nat)
    (digitsFit : ∀ index, digits index < 2 ^ 50)
    (m2 : Fin 8 → Nat)
    (m2Fit : ∀ index,
      m2 index < 2 * listAt8 Generated.moduli2 index.val)
    (m2Represent : ∀ index,
      Nat.ModEq (listAt8 Generated.moduli2 index.val)
        (m2 index)
        (r2InputValue (convertToCorrectedRepresentative digits))) :
    (generatedConvertFromValue m2 : ZMod Generated.target_modulus) =
      canonicalDigitValue digits := by
  let representative := convertToCorrectedRepresentative digits
  have representativeBound : representative <
      18 * Generated.target_modulus := by
    have repairedBound := convert_to_corrected_range digits digitsFit
    dsimp [representative]
    nlinarith
  have exported := generated_convert_from_represents_representative
    representative m2 representativeBound m2Fit m2Represent
  have imported := convert_to_corrected_represents_digits digits
  have convertToRelation := conversion_multiplier_relations.1
  calc
    (generatedConvertFromValue m2 : ZMod Generated.target_modulus) =
        Generated.convert_from_postmult * representative := exported
    _ = Generated.convert_from_postmult *
        (Generated.convert_to_premult * canonicalDigitValue digits) := by
      rw [imported]
    _ = (Generated.M1 * Generated.convert_from_postmult) *
        canonicalDigitValue digits := by
      rw [convertToRelation]
      ring
    _ = canonicalDigitValue digits := by
      rw [target_inverse_factor_zmod]
      simp

/-- The selected generated tables implement canonical field multiplication:
canonical radix-2^50 import, one MatrixNoK multiplication, and canonical
radix-2^52 export compose to multiplication modulo the target field modulus. -/
theorem canonical_import_selected_multiplication_export
    (leftDigits rightDigits : Fin 8 → Nat)
    (leftM1 rightM1 leftM2 rightM2 : Fin 8 → Nat)
    (leftDigitsFit : ∀ index, leftDigits index < 2 ^ 50)
    (rightDigitsFit : ∀ index, rightDigits index < 2 ^ 50)
    (leftM1Fit : ∀ index, leftM1 index < 2 * modulus1At index)
    (rightM1Fit : ∀ index, rightM1 index < 2 * modulus1At index)
    (leftM2Fit : ∀ index, leftM2 index < 2 * modulus2At index)
    (rightM2Fit : ∀ index, rightM2 index < 2 * modulus2At index)
    (leftM1Represent : ∀ index,
      Nat.ModEq (modulus1At index) (leftM1 index)
        (convertToCorrectedRepresentative leftDigits * ifmaRadix *
          Generated.qr_premult))
    (rightM1Represent : ∀ index,
      Nat.ModEq (modulus1At index) (rightM1 index)
        (convertToCorrectedRepresentative rightDigits * ifmaRadix *
          Generated.qr_premult))
    (leftM2Represent : ∀ index,
      Nat.ModEq (modulus2At index) (leftM2 index)
        (convertToCorrectedRepresentative leftDigits * ifmaRadix *
          crossBasisInverse))
    (rightM2Represent : ∀ index,
      Nat.ModEq (modulus2At index) (rightM2 index)
        (convertToCorrectedRepresentative rightDigits * ifmaRadix *
          crossBasisInverse)) :
    let intermediateM1 := m1ProductResidues leftM1 rightM1
    let resultM2 := m2CombinedResidues leftM2 rightM2 intermediateM1
    (generatedConvertFromValue resultM2 :
        ZMod Generated.target_modulus) =
      canonicalDigitValue leftDigits * canonicalDigitValue rightDigits := by
  dsimp only
  let leftValue := convertToCorrectedRepresentative leftDigits
  let rightValue := convertToCorrectedRepresentative rightDigits
  let intermediateM1 := m1ProductResidues leftM1 rightM1
  let representative := r1ProductFromResidues
    leftValue rightValue intermediateM1
  let resultM2 := m2CombinedResidues leftM2 rightM2 intermediateM1
  have leftValueBound : leftValue < 40 * Generated.target_modulus := by
    have bound := convert_to_corrected_range leftDigits leftDigitsFit
    dsimp [leftValue]
    nlinarith
  have rightValueBound : rightValue < 40 * Generated.target_modulus := by
    have bound := convert_to_corrected_range rightDigits rightDigitsFit
    dsimp [rightValue]
    nlinarith
  have pipeline := selected_base_field_pipeline leftValue rightValue
    leftM1 rightM1 leftM2 rightM2 leftValueBound rightValueBound
    leftM1Fit rightM1Fit leftM2Fit rightM2Fit
    leftM1Represent rightM1Represent leftM2Represent rightM2Represent
  have representativeBound : representative <
      18 * Generated.target_modulus := pipeline.1
  have productRelation := pipeline.2.1
  have resultM2Facts := pipeline.2.2.2
  have resultM2Fit : ∀ index,
      resultM2 index < 2 * listAt8 Generated.moduli2 index.val := by
    intro index
    simpa [modulus2At] using (resultM2Facts index).1
  have resultM2Represent : ∀ index,
      Nat.ModEq (listAt8 Generated.moduli2 index.val)
        (resultM2 index) (r2InputValue representative) := by
    intro index
    have represented := (resultM2Facts index).2
    have inputRelation := r2_input_represents_m2_encoding representative index
    rw [modulus2At] at represented inputRelation
    exact represented.trans inputRelation.symm
  have exported := generated_convert_from_represents_representative
    representative resultM2 representativeBound resultM2Fit resultM2Represent
  have productModEq : Nat.ModEq Generated.target_modulus
      (representative * Generated.M1) (leftValue * rightValue) :=
    productRelation
  have productCast :=
    (ZMod.natCast_eq_natCast_iff
      (representative * Generated.M1) (leftValue * rightValue)
      Generated.target_modulus).2 productModEq
  have leftImported := convert_to_corrected_represents_digits leftDigits
  have rightImported := convert_to_corrected_represents_digits rightDigits
  have convertToRelation := conversion_multiplier_relations.1
  have leftRepresents :
      (leftValue : ZMod Generated.target_modulus) =
        Generated.M1 * canonicalDigitValue leftDigits := by
    dsimp [leftValue]
    rw [leftImported, convertToRelation]
  have rightRepresents :
      (rightValue : ZMod Generated.target_modulus) =
        Generated.M1 * canonicalDigitValue rightDigits := by
    dsimp [rightValue]
    rw [rightImported, convertToRelation]
  calc
    (generatedConvertFromValue resultM2 :
        ZMod Generated.target_modulus) =
      Generated.convert_from_postmult * representative := exported
    _ = canonicalDigitValue leftDigits * canonicalDigitValue rightDigits :=
      export_scaling_composes
        (leftValue : ZMod Generated.target_modulus)
        (rightValue : ZMod Generated.target_modulus)
        (representative : ZMod Generated.target_modulus)
        (canonicalDigitValue leftDigits : ZMod Generated.target_modulus)
        (canonicalDigitValue rightDigits : ZMod Generated.target_modulus)
        leftRepresents rightRepresents (by simpa using productCast)

end CatBLST.AVX.MatrixNoK
