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

theorem generated_convert_to_estimate_lt_word (row : Fin 8) :
    listAt8 Generated.convert_to_shifted_quotient_estimations row.val <
      uint64Modulus := by
  fin_cases row <;> decide

theorem generated_convert_to_quotient_lt_2pow55
    (digits : Fin 8 → Nat)
    (digitsFit : ∀ index, digits index < 2 ^ 50) :
    generatedConvertToEstimatedQuotient digits < 2 ^ 55 := by
  have digitsFitRadix : ∀ index, digits index < ifmaRadix := by
    intro index
    exact (digitsFit index).trans (by norm_num [ifmaRadix])
  have words := quotient_words_fit digits
    (fun row =>
      listAt8 Generated.convert_to_shifted_quotient_estimations row.val)
    digitsFitRadix generated_convert_to_estimate_lt_word
  simpa [generatedConvertToEstimatedQuotient, quotientNumerator] using
    words.2.1

theorem generated_convert_to_matrix_entry_lt_radix (row column : Fin 8) :
    matrixEntry Generated.convert_to_rns_mat row.val column.val <
      ifmaRadix := by
  fin_cases row <;> fin_cases column <;> decide

theorem generated_convert_to_matrix_entry_lt_modulus (row column : Fin 8) :
    matrixEntry Generated.convert_to_rns_mat row.val column.val <
      modulus2At column := by
  fin_cases row <;> fin_cases column <;> decide

def convertToCorrectionShift (column : Fin 8) : Nat :=
  listAt8 Generated.convert_to_correction column.val * ifmaRadix %
    modulus2At column

def generatedConvertToShiftedWideValue
    (digits : Fin 8 → Nat) (column : Fin 8) : Nat :=
  let quotient := generatedConvertToEstimatedQuotient digits
  matrixProductValue digits Generated.convert_to_rns_mat column +
    listAt8 Generated.convert_to_correction column.val *
      (quotient % ifmaRadix) +
    convertToCorrectionShift column * (quotient / ifmaRadix)

def generatedConvertToFinalLow
    (digits : Fin 8 → Nat) (column : Fin 8) : Nat :=
  let quotient := generatedConvertToEstimatedQuotient digits
  lowMatrixAccumulator digits Generated.convert_to_rns_mat column +
    lowDigitProduct
      (listAt8 Generated.convert_to_correction column.val) quotient +
    lowDigitProduct (convertToCorrectionShift column)
      (quotient / ifmaRadix)

def generatedConvertToFinalHigh
    (digits : Fin 8 → Nat) (column : Fin 8) : Nat :=
  let quotient := generatedConvertToEstimatedQuotient digits
  highMatrixAccumulator digits Generated.convert_to_rns_mat column +
    highDigitProduct
      (listAt8 Generated.convert_to_correction column.val) quotient +
    highDigitProduct (convertToCorrectionShift column)
      (quotient / ifmaRadix)

/-- The split IFMA streams in the actual ConvertTo shifted schedule
reconstruct its mathematical wide value exactly. -/
theorem generated_convert_to_accumulators_reconstruct_wide
    (digits : Fin 8 → Nat)
    (digitsFit : ∀ index, digits index < 2 ^ 50)
    (column : Fin 8) :
    generatedConvertToFinalLow digits column + ifmaRadix *
        generatedConvertToFinalHigh digits column =
      generatedConvertToShiftedWideValue digits column := by
  have digitsFitRadix : ∀ index, digits index < ifmaRadix := by
    intro index
    exact (digitsFit index).trans (by norm_num [ifmaRadix])
  have matrixReconstruction := matrix_accumulators_reconstruct_product
    digits Generated.convert_to_rns_mat column digitsFitRadix
    (fun row => generated_convert_to_matrix_entry_lt_radix row column)
  have correctionFits :
      listAt8 Generated.convert_to_correction column.val < ifmaRadix := by
    fin_cases column <;>
      norm_num [Generated.convert_to_correction, listAt8, ifmaRadix]
  have shiftFits : convertToCorrectionShift column < ifmaRadix := by
    have modulusFits : modulus2At column < ifmaRadix := by
      fin_cases column <;>
        norm_num [modulus2At, Generated.moduli2, listAt8, ifmaRadix]
    exact (Nat.mod_lt _ (by
      fin_cases column <;>
        norm_num [modulus2At, Generated.moduli2, listAt8])).trans modulusFits
  have quotientBound :=
    generated_convert_to_quotient_lt_2pow55 digits digitsFit
  have quotientHighFits :
      generatedConvertToEstimatedQuotient digits / ifmaRadix < ifmaRadix := by
    apply (Nat.div_lt_iff_lt_mul (by norm_num [ifmaRadix])).2
    exact quotientBound.trans (by norm_num [ifmaRadix])
  have firstProduct := digit_product_reconstruction
    (listAt8 Generated.convert_to_correction column.val)
    (generatedConvertToEstimatedQuotient digits)
  rw [Nat.mod_eq_of_lt correctionFits] at firstProduct
  have secondProduct := digit_product_reconstruction
    (convertToCorrectionShift column)
    (generatedConvertToEstimatedQuotient digits / ifmaRadix)
  rw [Nat.mod_eq_of_lt shiftFits,
    Nat.mod_eq_of_lt quotientHighFits] at secondProduct
  unfold generatedConvertToFinalLow generatedConvertToFinalHigh
    generatedConvertToShiftedWideValue
  dsimp only
  nlinarith

theorem generated_convert_to_shifted_wide_represents_lane
    (digits : Fin 8 → Nat) (column : Fin 8) :
    (generatedConvertToShiftedWideValue digits column :
        ZMod (listAt8 Generated.moduli2 column.val)) =
      generatedConvertToLaneValue digits column := by
  unfold generatedConvertToShiftedWideValue generatedConvertToLaneValue
    matrixProductValue
  dsimp only
  push_cast
  have shiftCast :
      (convertToCorrectionShift column :
          ZMod (listAt8 Generated.moduli2 column.val)) =
        listAt8 Generated.convert_to_correction column.val * ifmaRadix := by
    unfold convertToCorrectionShift modulus2At
    rw [ZMod.natCast_mod]
    push_cast
    rfl
  have quotientSplit := Nat.mod_add_div
    (generatedConvertToEstimatedQuotient digits) ifmaRadix
  have quotientSplitCast :
      (generatedConvertToEstimatedQuotient digits :
          ZMod (listAt8 Generated.moduli2 column.val)) =
        (generatedConvertToEstimatedQuotient digits % ifmaRadix :
          ZMod (listAt8 Generated.moduli2 column.val)) +
        ifmaRadix *
          (generatedConvertToEstimatedQuotient digits / ifmaRadix :
            ZMod (listAt8 Generated.moduli2 column.val)) := by
    have castEquality := congrArg
      (fun value : Nat =>
        (value : ZMod (listAt8 Generated.moduli2 column.val)))
      quotientSplit.symm
    push_cast at castEquality
    exact castEquality
  rw [shiftCast, quotientSplitCast]
  ring

theorem convert_to_correction_and_shift_lt_modulus (column : Fin 8) :
    listAt8 Generated.convert_to_correction column.val < modulus2At column ∧
      convertToCorrectionShift column < modulus2At column := by
  constructor
  · fin_cases column <;>
      norm_num [Generated.convert_to_correction, modulus2At,
        Generated.moduli2, listAt8]
  · unfold convertToCorrectionShift
    exact Nat.mod_lt _ (by
      fin_cases column <;>
        norm_num [modulus2At, Generated.moduli2, listAt8])

theorem generated_convert_to_ifma_accumulators_fit_word
    (digits : Fin 8 → Nat) (column : Fin 8) :
    generatedConvertToFinalLow digits column < uint64Modulus ∧
      generatedConvertToFinalHigh digits column < uint64Modulus := by
  have lowFits :
      lowMatrixAccumulator digits Generated.convert_to_rns_mat column <
        8 * ifmaRadix := by
    unfold lowMatrixAccumulator
    exact fin8_sum_lt (fun index =>
      (digit_products_fit (digits index)
        (matrixEntry Generated.convert_to_rns_mat index.val column.val)).1)
  have highFits :
      highMatrixAccumulator digits Generated.convert_to_rns_mat column <
        8 * ifmaRadix := by
    unfold highMatrixAccumulator
    exact fin8_sum_lt (fun index =>
      (digit_products_fit (digits index)
        (matrixEntry Generated.convert_to_rns_mat index.val column.val)).2)
  simpa [generatedConvertToFinalLow, generatedConvertToFinalHigh] using
    shifted_matrix_correction_fits_word
      (highMatrixAccumulator digits Generated.convert_to_rns_mat column)
      (lowMatrixAccumulator digits Generated.convert_to_rns_mat column)
      (listAt8 Generated.convert_to_correction column.val)
      (convertToCorrectionShift column)
      (generatedConvertToEstimatedQuotient digits)
      highFits lowFits

theorem convert_to_matrix_product_lt_9_square
    (digits : Fin 8 → Nat)
    (digitsFit : ∀ index, digits index < 2 ^ 50)
    (column : Fin 8) :
    matrixProductValue digits Generated.convert_to_rns_mat column <
      9 * (modulus2At column) ^ 2 := by
  have eachFits : ∀ row : Fin 8,
      digits row *
          matrixEntry Generated.convert_to_rns_mat row.val column.val <
        2 ^ 50 * modulus2At column := by
    intro row
    have entryBound :=
      generated_convert_to_matrix_entry_lt_modulus row column
    have modulusPositive : 0 < modulus2At column := by
      fin_cases column <;>
        norm_num [modulus2At, Generated.moduli2, listAt8]
    have scaledDigitBound :
        digits row * modulus2At column < 2 ^ 50 * modulus2At column :=
      (Nat.mul_lt_mul_right modulusPositive).2 (digitsFit row)
    calc
      _ ≤ digits row * modulus2At column :=
        Nat.mul_le_mul_left _ (le_of_lt entryBound)
      _ < 2 ^ 50 * modulus2At column := scaledDigitBound
  have sumFits := fin8_sum_lt eachFits
  unfold matrixProductValue
  calc
    _ < 8 * (2 ^ 50 * modulus2At column) := sumFits
    _ < 9 * (modulus2At column) ^ 2 := by
      fin_cases column <;>
        norm_num [modulus2At, Generated.moduli2, listAt8]

theorem convert_to_shifted_correction_lt_5_square
    (digits : Fin 8 → Nat)
    (digitsFit : ∀ index, digits index < 2 ^ 50)
    (column : Fin 8) :
    let quotient := generatedConvertToEstimatedQuotient digits
    listAt8 Generated.convert_to_correction column.val *
          (quotient % ifmaRadix) +
        convertToCorrectionShift column * (quotient / ifmaRadix) <
      5 * (modulus2At column) ^ 2 := by
  dsimp only
  have quotientBound :=
    generated_convert_to_quotient_lt_2pow55 digits digitsFit
  have quotientLow :
      generatedConvertToEstimatedQuotient digits % ifmaRadix < ifmaRadix :=
    Nat.mod_lt _ (by norm_num [ifmaRadix])
  have quotientHigh :
      generatedConvertToEstimatedQuotient digits / ifmaRadix < 8 := by
    apply (Nat.div_lt_iff_lt_mul (by norm_num [ifmaRadix])).2
    norm_num [ifmaRadix]
    exact quotientBound
  have correctionBounds := convert_to_correction_and_shift_lt_modulus column
  have modulusWindow : ifmaRadix + 8 < 5 * modulus2At column := by
    fin_cases column <;>
      norm_num [modulus2At, Generated.moduli2, listAt8, ifmaRadix]
  nlinarith

theorem generated_convert_to_shifted_wide_lt_16_square
    (digits : Fin 8 → Nat)
    (digitsFit : ∀ index, digits index < 2 ^ 50)
    (column : Fin 8) :
    generatedConvertToShiftedWideValue digits column <
      16 * (modulus2At column) ^ 2 := by
  have matrixBound := convert_to_matrix_product_lt_9_square
    digits digitsFit column
  have correctionBound := convert_to_shifted_correction_lt_5_square
    digits digitsFit column
  unfold generatedConvertToShiftedWideValue
  nlinarith

def canonicalImportedM2Residue
    (digits : Fin 8 → Nat) (column : Fin 8) : Nat :=
  let wide := generatedConvertToShiftedWideValue digits column
  selectedReduceMasked (wide / ifmaRadix) (wide % ifmaRadix)
    (modulus2At column) (reducer2FactorAt column)

/--
The actual ConvertTo shifted accumulators followed by the selected `800`
masked Montgomery reducer establish the standard M2 residue postcondition.
The tighter canonical-input proof bound is retained as evidence of margin;
the reducer parameters remain those instantiated by the C++ source.
-/
theorem canonical_imported_m2_residue_contract
    (digits : Fin 8 → Nat)
    (digitsFit : ∀ index, digits index < 2 ^ 50) :
    ∀ column,
      (canonicalImportedM2Residue digits column :
          ZMod (modulus2At column)) * ifmaRadix =
        generatedConvertToShiftedWideValue digits column ∧
      canonicalImportedM2Residue digits column <
        2 * modulus2At column := by
  intro column
  let wide := generatedConvertToShiftedWideValue digits column
  have tightWideBound : wide < 16 * (modulus2At column) ^ 2 := by
    simpa [wide] using
      generated_convert_to_shifted_wide_lt_16_square digits digitsFit column
  have wideBound : wide ≤ 800 * modulus2At column * modulus2At column := by
    have modulusPositive : 0 < modulus2At column := by
      fin_cases column <;>
        norm_num [modulus2At, Generated.moduli2, listAt8]
    nlinarith [tightWideBound]
  have split := Nat.mod_add_div wide ifmaRadix
  have splitReordered :
      (wide / ifmaRadix) * ifmaRadix + wide % ifmaRadix = wide := by
    simpa [Nat.add_comm, Nat.mul_comm] using split
  have reducer := selected_masked_reducer_contract
    (wide / ifmaRadix) (wide % ifmaRadix)
    (modulus2At column) (reducer2FactorAt column) 800 200
    (modulus2_member column) (reducer2_factor_relation column)
    (by norm_num) (by
      rw [splitReordered]
      simpa [Nat.mul_assoc] using wideBound)
    (by norm_num)
    (by norm_num [elementRadix, uint64Modulus])
  have splitCast := congrArg
    (fun value : Nat => (value : ZMod (modulus2At column))) splitReordered
  push_cast at splitCast
  unfold canonicalImportedM2Residue
  dsimp only
  constructor
  · rw [reducer.1]
    simpa [wide] using splitCast
  · exact reducer.2

theorem canonical_imported_m2_residue_represents
    (digits : Fin 8 → Nat)
    (digitsFit : ∀ index, digits index < 2 ^ 50)
    (column : Fin 8) :
    Nat.ModEq (modulus2At column)
      (canonicalImportedM2Residue digits column)
      (r2InputValue (convertToCorrectedRepresentative digits)) := by
  let representative := convertToCorrectedRepresentative digits
  have reducedContract :=
    (canonical_imported_m2_residue_contract digits digitsFit column).1
  change
    (canonicalImportedM2Residue digits column :
        ZMod (listAt8 Generated.moduli2 column.val)) * ifmaRadix =
      generatedConvertToShiftedWideValue digits column at reducedContract
  have reducedContract' :
      (canonicalImportedM2Residue digits column :
          ZMod (listAt8 Generated.moduli2 column.val)) * ifmaRadix =
        generatedConvertToShiftedWideValue digits column := by
    exact reducedContract
  have wideCast :
      (generatedConvertToShiftedWideValue digits column :
          ZMod (listAt8 Generated.moduli2 column.val)) =
        generatedConvertToLaneValue digits column := by
    exact generated_convert_to_shifted_wide_represents_lane digits column
  have laneCast :
      generatedConvertToLaneValue digits column =
        (Generated.convert_to_postmult :
          ZMod (listAt8 Generated.moduli2 column.val)) *
          representative := by
    simpa only [representative] using
      generated_convert_to_lane_represents_import digits column
  have outputTimesRadix :
      (canonicalImportedM2Residue digits column :
          ZMod (listAt8 Generated.moduli2 column.val)) *
          ifmaRadix =
        (Generated.convert_to_postmult :
          ZMod (listAt8 Generated.moduli2 column.val)) *
          representative := by
    rw [reducedContract', wideCast, laneCast]
  have targetModEq : Nat.ModEq
      (listAt8 Generated.moduli2 column.val)
      (r2InputValue representative)
      (representative * ifmaRadix * crossBasisInverse) := by
    unfold r2InputValue
    change
      ((representative * montgomeryRadix * crossBasisInverse) %
          Generated.M2) % listAt8 Generated.moduli2 column.val =
        (representative * ifmaRadix * crossBasisInverse) %
          listAt8 Generated.moduli2 column.val
    rw [Nat.mod_mod_of_dvd _ (by
      simpa only [modulus2At] using modulus2_dvd_product column)]
    simp only [montgomeryRadix, ifmaRadix]
  have targetCast :=
    (ZMod.natCast_eq_natCast_iff _ _
      (listAt8 Generated.moduli2 column.val)).2 targetModEq
  have postmultRelation := conversion_multiplier_relations.2.1
  have targetTimesRadix :
      (r2InputValue representative :
          ZMod (listAt8 Generated.moduli2 column.val)) * ifmaRadix =
        (Generated.convert_to_postmult :
          ZMod (listAt8 Generated.moduli2 column.val)) *
          representative := by
    rw [targetCast, postmultRelation]
    simp only [montgomeryRadix, ifmaRadix]
    push_cast
    ring
  have equalTimesRadix :
      (canonicalImportedM2Residue digits column :
          ZMod (listAt8 Generated.moduli2 column.val)) *
          ifmaRadix =
        (r2InputValue representative :
          ZMod (listAt8 Generated.moduli2 column.val)) *
          ifmaRadix := by
    exact outputTimesRadix.trans targetTimesRadix.symm
  have radixUnit := ifma_radix_is_unit_modulus2 column
  change IsUnit
    (ifmaRadix : ZMod (listAt8 Generated.moduli2 column.val)) at radixUnit
  have equalCast :=
    radixUnit.mul_right_cancel equalTimesRadix
  have result := (ZMod.natCast_eq_natCast_iff _ _
    (listAt8 Generated.moduli2 column.val)).1 equalCast
  simpa only [modulus2At] using result

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

/-- Canonical ConvertTo, its selected M2 reducer, and the checked-in R2
expansion table compose without an assumed intermediate residue contract. -/
theorem canonical_import_actual_r2_expansion_lane
    (digits : Fin 8 → Nat)
    (digitsFit : ∀ index, digits index < 2 ^ 50)
    (column : Fin 8) :
    generatedR2LaneValue (canonicalImportedM2Residue digits) column =
      (Generated.r2_postmult :
        ZMod (listAt8 Generated.moduli1 column.val)) *
        convertToCorrectedRepresentative digits := by
  exact canonical_import_r2_expansion_lane digits digitsFit
    (canonicalImportedM2Residue digits)
    (fun index => by
      simpa only [modulus2At] using
        (canonical_imported_m2_residue_contract digits digitsFit index).2)
    (fun index => by
      simpa only [modulus2At] using
        canonical_imported_m2_residue_represents digits digitsFit index)
    column

def canonicalImportedM1Residue
    (digits : Fin 8 → Nat) (column : Fin 8) : Nat :=
  m1ExpandedResidues (canonicalImportedM2Residue digits) column

/-- The complete canonical import produces both ordinary MatrixNoK residue
bases with the bounds and encodings required by multiplication. -/
theorem canonical_imported_m1_residue_contract
    (digits : Fin 8 → Nat)
    (digitsFit : ∀ index, digits index < 2 ^ 50) :
    ∀ column,
      canonicalImportedM1Residue digits column < 2 * modulus1At column ∧
      Nat.ModEq (modulus1At column)
        (canonicalImportedM1Residue digits column)
        (convertToCorrectedRepresentative digits * ifmaRadix *
          Generated.qr_premult) := by
  intro column
  have m2Fit : ∀ index,
      canonicalImportedM2Residue digits index < 2 * modulus2At index :=
    fun index =>
      (canonical_imported_m2_residue_contract digits digitsFit index).2
  have m2Represent : ∀ index,
      Nat.ModEq (modulus2At index)
        (canonicalImportedM2Residue digits index)
        (r2InputValue (convertToCorrectedRepresentative digits)) :=
    fun index => canonical_imported_m2_residue_represents
      digits digitsFit index
  have representativeBound :
      convertToCorrectedRepresentative digits <
        18 * Generated.target_modulus := by
    have tightBound := convert_to_corrected_range digits digitsFit
    nlinarith
  constructor
  · unfold canonicalImportedM1Residue
    exact (m1_expanded_residues_contract
      (canonicalImportedM2Residue digits) m2Fit column).2
  · unfold canonicalImportedM1Residue
    exact m1_expanded_residues_represent
      (convertToCorrectedRepresentative digits)
      (canonicalImportedM2Residue digits) representativeBound
      m2Fit m2Represent column

end CatBLST.AVX.MatrixNoK
