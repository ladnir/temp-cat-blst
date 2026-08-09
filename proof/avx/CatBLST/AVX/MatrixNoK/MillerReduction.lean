/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.MillerPhysical
import Mathlib.Tactic

/-!
# Four-product Miller completion and reduction

This composes the accepted MatrixNoK completion and two reducers with the
`16m^2` raw bound required by Miller `mul_accumulate` expressions.  It changes
no executable code and retains the exact generated offsets and matrices.
-/

namespace CatBLST.AVX.MatrixNoK

theorem miller_fp2_completed_m1_wide_fits
    (rawWide high low : Fin 8 → Nat)
    (rawFit : ∀ index,
      rawWide index ≤ millerFp2AccumulatorBound * modulus1At index *
        modulus1At index)
    (wideSplit : ∀ index, high index * ifmaRadix + low index =
      rawWide index + listAt8 Generated.max_neg_offset_m1 index.val) :
    ∀ index, high index * ifmaRadix + low index ≤
      delayedCombinedTypeBound * modulus1At index * modulus1At index := by
  intro index
  have offsetBound := fp2_max_neg_offset_m1_lane_bound index
  have modulusPositive : 0 < modulus1At index := by
    fin_cases index <;> norm_num [modulus1At, listAt8, Generated.moduli1]
  have offsetSquare :
      listAt8 Generated.max_neg_offset_m1 index.val ≤
        modulus1At index * modulus1At index := by
    nlinarith
  have rawBound := rawFit index
  have typeBound := miller_fp2_completed_m1_fits_selected_type
  rw [wideSplit index]
  nlinarith

theorem miller_fp2_completed_m1_residue_represents
    (value : Int) (rawWide high low : Fin 8 → Nat)
    (rawRepresents : ∀ index,
      (rawWide index : ZMod (modulus1At index)) =
        value * ifmaRadix * ifmaRadix * Generated.qr_premult *
          Generated.qr_premult)
    (rawFit : ∀ index,
      rawWide index ≤ millerFp2AccumulatorBound * modulus1At index *
        modulus1At index)
    (wideSplit : ∀ index, high index * ifmaRadix + low index =
      rawWide index + listAt8 Generated.max_neg_offset_m1 index.val) :
    ∀ index, Nat.ModEq (modulus1At index)
      (cppFp2CompletedM1Residue high low index)
      (signedR1ReducedValue (cppFp2CompletedExpression value)) := by
  intro index
  unfold cppFp2CompletedM1Residue
  exact fp2_completed_m1_reducer_represents value (rawWide index)
    (high index) (low index) index (rawRepresents index) (wideSplit index)
    (miller_fp2_completed_m1_wide_fits rawWide high low rawFit wideSplit index)

theorem miller_fp2_completed_m1_residue_fit
    (rawWide high low : Fin 8 → Nat)
    (rawFit : ∀ index,
      rawWide index ≤ millerFp2AccumulatorBound * modulus1At index *
        modulus1At index)
    (wideSplit : ∀ index, high index * ifmaRadix + low index =
      rawWide index + listAt8 Generated.max_neg_offset_m1 index.val) :
    ∀ index,
      cppFp2CompletedM1Residue high low index < 2 * modulus1At index := by
  intro index
  have reducer := selected_delayed_reducer_contract
    (high index) (low index) (modulus1At index) (reducer1FactorAt index)
    (modulus1_member index) (reducer1_factor_relation index)
    (miller_fp2_completed_m1_wide_fits rawWide high low rawFit wideSplit index)
  exact reducer.2.2

theorem miller_fp2_completed_m2_combined_wide_lt_33_square
    (rawWide m1Residues : Fin 8 → Nat)
    (rawFit : ∀ index,
      rawWide index ≤ millerFp2AccumulatorBound * modulus2At index *
        modulus2At index)
    (m1Fit : ∀ index, m1Residues index < 2 * modulus1At index)
    (column : Fin 8) :
    cppFp2CompletedM2CombinedWideValue rawWide m1Residues column <
      33 * modulus2At column * modulus2At column := by
  have matrixBound := r1_matrix_product_lt_16_square
    m1Residues (by simpa [modulus1At] using m1Fit) column
  have matrixBound' :
      matrixProductValue m1Residues Generated.r1_rns_mat column <
        16 * modulus2At column * modulus2At column := by
    simpa [modulus2At, Nat.pow_two, Nat.mul_assoc] using matrixBound
  have offsetBound := fp2_max_neg_offset_m2_lane_bound column
  have modulusPositive : 0 < modulus2At column := by
    fin_cases column <;> norm_num [modulus2At, listAt8, Generated.moduli2]
  have offsetSquare :
      listAt8 Generated.max_neg_offset_m2 column.val ≤
        modulus2At column * modulus2At column := by
    nlinarith
  have rawBound := rawFit column
  rw [selected_miller_fp2_accumulator_bound] at rawBound
  unfold cppFp2CompletedM2CombinedWideValue
  nlinarith

theorem miller_fp2_completed_m2_residue_represents
    (value : Int) (rawWide m1Residues : Fin 8 → Nat)
    (rawRepresents : ∀ index,
      (rawWide index : ZMod (modulus2At index)) =
        value * ifmaRadix * ifmaRadix * crossBasisInverse *
          crossBasisInverse)
    (rawFit : ∀ index,
      rawWide index ≤ millerFp2AccumulatorBound * modulus2At index *
        modulus2At index)
    (m1Fit : ∀ index, m1Residues index < 2 * modulus1At index)
    (m1Represent : ∀ index,
      Nat.ModEq (modulus1At index) (m1Residues index)
        (signedR1ReducedValue (cppFp2CompletedExpression value))) :
    ∀ column, Nat.ModEq (modulus2At column)
      (cppFp2CompletedM2Residue rawWide m1Residues column)
      (signedR2InputValue (r1SignedValueFromResidues
        (cppFp2CompletedExpression value) m1Residues)) := by
  intro column
  let modulus := modulus2At column
  let completed := cppFp2CompletedExpression value
  let representative := r1SignedValueFromResidues completed m1Residues
  let wide := cppFp2CompletedM2CombinedWideValue rawWide m1Residues column
  have wideBound : wide ≤ delayedCombinedTypeBound * modulus * modulus := by
    have tight := miller_fp2_completed_m2_combined_wide_lt_33_square
      rawWide m1Residues rawFit m1Fit column
    have typeBound : 33 ≤ delayedCombinedTypeBound := by
      norm_num [delayedCombinedTypeBound, delayedInputBound,
        delayedR1TypeBound, r1ExactAccumulationBound,
        multiplicationProductBound, standardBound]
    dsimp only [wide, modulus]
    nlinarith
  have split := Nat.mod_add_div wide ifmaRadix
  have splitReordered :
      (wide / ifmaRadix) * ifmaRadix + wide % ifmaRadix = wide := by
    simpa [Nat.add_comm, Nat.mul_comm] using split
  have reducer := selected_delayed_reducer_contract
    (wide / ifmaRadix) (wide % ifmaRadix) modulus
    (reducer2FactorAt column) (modulus2_member column)
    (reducer2_factor_relation column) (by
      rw [splitReordered]
      exact wideBound)
  have completedCast := fp2_completed_m2_wide_lane_encoding
    value (rawWide column) column (rawRepresents column)
  have matrixCast :
      (matrixProductValue m1Residues Generated.r1_rns_mat column :
          ZMod modulus) = generatedR1LaneValue m1Residues column := by
    unfold matrixProductValue generatedR1LaneValue modulus modulus2At
    push_cast
    rfl
  have rawMatrix := generated_r1_lane_represents_raw_crt m1Residues column
  have combinedCast :
      (wide : ZMod modulus) =
        (completed + Generated.target_modulus *
          r1RawLiftFromResidues m1Residues) * ifmaRadix * ifmaRadix *
            crossBasisInverse * crossBasisInverse := by
    unfold wide cppFp2CompletedM2CombinedWideValue
    push_cast
    rw [completedCast, matrixCast, rawMatrix]
    unfold completed
    rw [r1_postmult_relation]
    unfold r1RawLiftFromResidues
    push_cast
    simp only [montgomeryRadix, ifmaRadix]
    ring_nf
    ac_rfl
  have reconstruction := r1_signed_value_reconstruction
    completed m1Residues (by simpa [completed] using m1Represent)
  have reconstructionCast :
      (representative : ZMod modulus) * Generated.M1 =
        completed + Generated.target_modulus *
          r1RawLiftFromResidues m1Residues := by
    have castEquality := congrArg
      (fun result : Int => (result : ZMod modulus)) reconstruction
    unfold r1SignedValueNumerator at castEquality
    push_cast at castEquality
    simpa [representative] using castEquality
  have inverseCast := fp2_cross_basis_inverse_lane column
  have targetCast := signed_r2_input_lane representative column
  have splitCast := congrArg
    (fun result : Nat => (result : ZMod modulus)) splitReordered
  push_cast at splitCast
  have outputTimesRadix :
      (cppFp2CompletedM2Residue rawWide m1Residues column : ZMod modulus) *
          ifmaRadix =
        (signedR2InputValue representative : ZMod modulus) * ifmaRadix := by
    unfold cppFp2CompletedM2Residue
    dsimp only
    rw [reducer.1, splitCast, combinedCast, targetCast]
    rw [← reconstructionCast]
    change
      (representative : ZMod modulus) * Generated.M1 * ifmaRadix *
          ifmaRadix * crossBasisInverse * crossBasisInverse =
        representative * ifmaRadix * crossBasisInverse * ifmaRadix
    have inverseCast' :
        (Generated.M1 : ZMod modulus) * crossBasisInverse = 1 := by
      simpa only [Nat.cast_mul, Nat.cast_one] using inverseCast
    rw [show
      (representative : ZMod modulus) * Generated.M1 * ifmaRadix *
          ifmaRadix * crossBasisInverse * crossBasisInverse =
        representative * ifmaRadix * ifmaRadix * crossBasisInverse *
          (Generated.M1 * crossBasisInverse) by ring,
      inverseCast']
    ring
  have radixUnit : IsUnit (ifmaRadix : ZMod modulus) := by
    simpa [modulus] using ifma_radix_is_unit_modulus2 column
  have outputCast := radixUnit.mul_right_cancel outputTimesRadix
  exact (ZMod.natCast_eq_natCast_iff _ _ modulus).1 outputCast

theorem miller_fp2_completed_m2_residue_fit
    (rawWide m1Residues : Fin 8 → Nat)
    (rawFit : ∀ index,
      rawWide index ≤ millerFp2AccumulatorBound * modulus2At index *
        modulus2At index)
    (m1Fit : ∀ index, m1Residues index < 2 * modulus1At index) :
    ∀ column,
      cppFp2CompletedM2Residue rawWide m1Residues column <
        2 * modulus2At column := by
  intro column
  let modulus := modulus2At column
  let wide := cppFp2CompletedM2CombinedWideValue rawWide m1Residues column
  have tight := miller_fp2_completed_m2_combined_wide_lt_33_square
    rawWide m1Residues rawFit m1Fit column
  have wideBound : wide ≤ delayedCombinedTypeBound * modulus * modulus := by
    have typeBound : 33 ≤ delayedCombinedTypeBound := by
      norm_num [delayedCombinedTypeBound, delayedInputBound,
        delayedR1TypeBound, r1ExactAccumulationBound,
        multiplicationProductBound, standardBound]
    dsimp only [wide, modulus]
    nlinarith
  have split := Nat.mod_add_div wide ifmaRadix
  have splitReordered :
      (wide / ifmaRadix) * ifmaRadix + wide % ifmaRadix = wide := by
    simpa [Nat.add_comm, Nat.mul_comm] using split
  have reducer := selected_delayed_reducer_contract
    (wide / ifmaRadix) (wide % ifmaRadix) modulus
    (reducer2FactorAt column) (modulus2_member column)
    (reducer2_factor_relation column) (by
      rw [splitReordered]
      exact wideBound)
  simpa [cppFp2CompletedM2Residue, wide, modulus] using reducer.2.2

theorem miller_fp2_completed_r2_representative_is_exact
    (value : Int) (m1Residues m2Residues : Fin 8 → Nat)
    (valueBound : |value| <
      4 * (fp2BaseRedundantMultiple * Generated.target_modulus) ^ 2)
    (m1ResiduesFit : ∀ index,
      m1Residues index < 2 * listAt8 Generated.moduli1 index.val)
    (m1ResiduesRepresent : ∀ index,
      Nat.ModEq (modulus1At index) (m1Residues index)
        (signedR1ReducedValue (cppFp2CompletedExpression value)))
    (m2ResiduesFit : ∀ index,
      m2Residues index < 2 * listAt8 Generated.moduli2 index.val)
    (m2ResiduesRepresent : ∀ index,
      Nat.ModEq (modulus2At index) (m2Residues index)
        (signedR2InputValue (r1SignedValueFromResidues
          (cppFp2CompletedExpression value) m1Residues))) :
    correctedCrtRepresentative m2Residues r2ShiftedFactor Generated.M2 =
      r1SignedValueFromResidues
        (cppFp2CompletedExpression value) m1Residues := by
  let representative := r1SignedValueFromResidues
    (cppFp2CompletedExpression value) m1Residues
  have representativeRange := miller_fp2_completed_r1_representative_range
    value m1Residues valueBound m1ResiduesFit m1ResiduesRepresent
  have representativeNonnegative : (0 : Int) ≤ representative := by
    simpa [representative] using representativeRange.1
  have representativeBound :
      representative < 22 * Generated.target_modulus := by
    simpa [representative] using representativeRange.2
  have representativeNatBound :
      representative.toNat < 22 * Generated.target_modulus := by
    rw [← Int.ofNat_lt]
    simpa [Int.toNat_of_nonneg representativeNonnegative] using
      representativeBound
  have representedNat : ∀ index,
      Nat.ModEq (listAt8 Generated.moduli2 index.val)
        (m2Residues index) (r2InputValue representative.toNat) := by
    intro index
    simpa [representative, modulus2At,
      signed_r2_input_of_nonnegative representative
        representativeNonnegative] using m2ResiduesRepresent index
  have exact := r2_standard_corrected_representative_is_exact_22
    representative.toNat m2Residues representativeNatBound
      m2ResiduesFit representedNat
  simpa [representative, Int.toNat_of_nonneg representativeNonnegative] using exact

theorem miller_fp2_selected_completed_reduction_is_exact
    (value : Int)
    (rawM1 m1High m1Low rawM2 : Fin 8 → Nat)
    (valueBound : |value| <
      4 * (fp2BaseRedundantMultiple * Generated.target_modulus) ^ 2)
    (rawM1Represents : ∀ index,
      (rawM1 index : ZMod (modulus1At index)) =
        value * ifmaRadix * ifmaRadix * Generated.qr_premult *
          Generated.qr_premult)
    (rawM1Fit : ∀ index,
      rawM1 index ≤ millerFp2AccumulatorBound * modulus1At index *
        modulus1At index)
    (m1WideSplit : ∀ index,
      m1High index * ifmaRadix + m1Low index =
        rawM1 index + listAt8 Generated.max_neg_offset_m1 index.val)
    (rawM2Represents : ∀ index,
      (rawM2 index : ZMod (modulus2At index)) =
        value * ifmaRadix * ifmaRadix * crossBasisInverse *
          crossBasisInverse)
    (rawM2Fit : ∀ index,
      rawM2 index ≤ millerFp2AccumulatorBound * modulus2At index *
        modulus2At index) :
    let m1Residues := cppFp2CompletedM1Residue m1High m1Low
    correctedCrtRepresentative
        (cppFp2CompletedM2Residue rawM2 m1Residues)
        r2ShiftedFactor Generated.M2 =
      r1SignedValueFromResidues
        (cppFp2CompletedExpression value) m1Residues := by
  dsimp only
  let m1Residues := cppFp2CompletedM1Residue m1High m1Low
  have m1Fit := miller_fp2_completed_m1_residue_fit
    rawM1 m1High m1Low rawM1Fit m1WideSplit
  have m1Represent := miller_fp2_completed_m1_residue_represents value rawM1
    m1High m1Low rawM1Represents rawM1Fit m1WideSplit
  apply miller_fp2_completed_r2_representative_is_exact value m1Residues
    (cppFp2CompletedM2Residue rawM2 m1Residues) valueBound
    (by simpa [m1Residues, modulus1At] using m1Fit)
    (by simpa [m1Residues] using m1Represent)
  · exact miller_fp2_completed_m2_residue_fit rawM2 m1Residues rawM2Fit
      (by simpa [m1Residues] using m1Fit)
  · exact miller_fp2_completed_m2_residue_represents value rawM2 m1Residues
      rawM2Represents rawM2Fit
      (by simpa [m1Residues] using m1Fit)
      (by simpa [m1Residues] using m1Represent)

end CatBLST.AVX.MatrixNoK
