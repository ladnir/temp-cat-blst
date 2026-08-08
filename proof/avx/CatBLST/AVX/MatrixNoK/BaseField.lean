/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.CRT
import Mathlib.Tactic

/-!
# MatrixNoK base-field multiplication semantics

This file closes the positive-representative R2 expansion used by the selected
base-field multiplication.  The correction-free R1 stage first produces a
nonnegative representative below `18p`; the fixed-point R2 correction is then
proved exact for that range, so expansion reconstructs the same representative
in M1.
-/

namespace CatBLST.AVX.MatrixNoK

open scoped BigOperators

def r2InputValue (representative : Nat) : Nat :=
  representative * montgomeryRadix * crossBasisInverse % Generated.M2

def r2Residues (representative : Nat) (index : Fin 8) : Nat :=
  r2InputValue representative % listAt8 Generated.moduli2 index.val

def r2ShiftedFactor (index : Fin 8) : Nat :=
  Generated.r2_premult * listAt8 crtFactors2 index.val % Generated.M2

def r2RawLift (representative : Nat) : Nat :=
  scaledCrtReconstruction8 Generated.moduli2 crtFactors2
    Generated.r2_premult (r2InputValue representative)

theorem r2_raw_lift_eq_exact_numerator (representative : Nat) :
    r2RawLift representative =
      exactCrtNumerator (r2Residues representative) r2ShiftedFactor := by
  unfold r2RawLift r2Residues r2ShiftedFactor
  rw [M2_eq_product]
  simp [scaledCrtReconstruction8, exactCrtNumerator, Fin.sum_univ_succ]
  ring

theorem r2_raw_lift_mod_input (representative : Nat) :
    r2RawLift representative % Generated.M2 =
      (Generated.r2_premult * r2InputValue representative) % Generated.M2 := by
  have scaled := scaled_crt_reconstruction_mod_product
    Generated.moduli2 crtFactors2 Generated.r2_premult
      (r2InputValue representative)
  have reconstruction := crt_reconstruction2_mod_product
    (r2InputValue representative)
  have multiplied :=
    (Nat.ModEq.refl (n := Generated.M2) Generated.r2_premult).mul
      reconstruction
  rw [M2_eq_product] at multiplied
  unfold r2RawLift
  rw [M2_eq_product]
  exact scaled.trans multiplied

theorem r2_input_scale_relation (representative : Nat) :
    (Generated.r2_premult * r2InputValue representative) % Generated.M2 =
      representative % Generated.M2 := by
  change Nat.ModEq Generated.M2
    (Generated.r2_premult * r2InputValue representative) representative
  have inputRelation : Nat.ModEq Generated.M2
      (r2InputValue representative)
      (representative * montgomeryRadix * crossBasisInverse) := by
    unfold r2InputValue
    change
      (representative * montgomeryRadix * crossBasisInverse % Generated.M2) %
          Generated.M2 =
        (representative * montgomeryRadix * crossBasisInverse) % Generated.M2
    rw [Nat.mod_mod]
  have premultRelation : Nat.ModEq Generated.M2
      (Generated.r2_premult * montgomeryRadix) Generated.M1 :=
    r2_premult_relation
  have crossRelation : Nat.ModEq Generated.M2
      (Generated.M1 * crossBasisInverse) 1 := by
    exact cross_basis_inverse_relation
  have scaledInput :=
    (Nat.ModEq.refl (n := Generated.M2) Generated.r2_premult).mul
      inputRelation
  have rotated :=
    ((Nat.ModEq.refl (n := Generated.M2) representative).mul
      premultRelation).mul
        (Nat.ModEq.refl (n := Generated.M2) crossBasisInverse)
  have unrotated :=
    (Nat.ModEq.refl (n := Generated.M2) representative).mul crossRelation
  have scaledOrder :
      Generated.r2_premult *
          (representative * montgomeryRadix * crossBasisInverse) =
        representative * (Generated.r2_premult * montgomeryRadix) *
          crossBasisInverse := by ring
  rw [scaledOrder] at scaledInput
  have targetOrder : representative * Generated.M1 * crossBasisInverse =
      representative * (Generated.M1 * crossBasisInverse) := by ring
  rw [targetOrder] at rotated
  rw [Nat.mul_one] at unrotated
  exact scaledInput.trans (rotated.trans unrotated)

theorem r2_raw_lift_remainder
    (representative : Nat) (representativeBound : representative < Generated.M2) :
    r2RawLift representative % Generated.M2 = representative := by
  rw [r2_raw_lift_mod_input, r2_input_scale_relation]
  exact Nat.mod_eq_of_lt representativeBound

theorem r2_residue_lt_element_radix
    (representative : Nat) (index : Fin 8) :
    r2Residues representative index < elementRadix := by
  fin_cases index
  all_goals
    unfold r2Residues listAt8
    exact (Nat.mod_lt _ (by norm_num [Generated.moduli2])).trans
      (by norm_num [Generated.moduli2, elementRadix])

theorem r2_scalar_sum_bound (representative : Nat) :
    (∑ index, r2Residues representative index) < 8 * elementRadix :=
  fin8_sum_lt (r2_residue_lt_element_radix representative)

theorem r2_positive_margin
    (representative : Nat)
    (representativeBound : representative < 18 * Generated.target_modulus) :
    representative * uint64Modulus +
        (Generated.M2 - 1) * ∑ index, r2Residues representative index <
      Generated.M2 * uint64Modulus := by
  have sumBound := r2_scalar_sum_bound representative
  have envelope :
      (18 * Generated.target_modulus) * uint64Modulus +
          (Generated.M2 - 1) * (8 * elementRadix) <
        Generated.M2 * uint64Modulus := by
    norm_num [Generated.target_modulus, Generated.M2,
      uint64Modulus, elementRadix]
  nlinarith

theorem r2_estimated_quotient_is_exact
    (representative : Nat)
    (representativeBound : representative < 18 * Generated.target_modulus) :
    estimatedCrtQuotient (r2Residues representative) r2ShiftedFactor
        Generated.M2 =
      exactCrtQuotient (r2Residues representative) r2ShiftedFactor
        Generated.M2 := by
  have belowM2 : representative < Generated.M2 :=
    representativeBound.trans (by
      norm_num [Generated.target_modulus, Generated.M2])
  have remainder := r2_raw_lift_remainder representative belowM2
  have exactNumerator := r2_raw_lift_eq_exact_numerator representative
  have scalarFits : ∀ index, r2Residues representative index < ifmaRadix := by
    intro index
    exact (r2_residue_lt_element_radix representative index).trans
      (by norm_num [elementRadix, ifmaRadix])
  apply estimated_quotient_exact_of_margin
    (r2Residues representative) r2ShiftedFactor Generated.M2
    (by norm_num [Generated.M2]) scalarFits
  rw [← exactNumerator, remainder]
  exact r2_positive_margin representative representativeBound

theorem r2_corrected_representative_is_exact
    (representative : Nat)
    (representativeBound : representative < 18 * Generated.target_modulus) :
    correctedCrtRepresentative (r2Residues representative) r2ShiftedFactor
        Generated.M2 = representative := by
  have belowM2 : representative < Generated.M2 :=
    representativeBound.trans (by
      norm_num [Generated.target_modulus, Generated.M2])
  have exactNumerator := r2_raw_lift_eq_exact_numerator representative
  have remainder := r2_raw_lift_remainder representative belowM2
  have split := Nat.mod_add_div (r2RawLift representative) Generated.M2
  unfold correctedCrtRepresentative
  rw [r2_estimated_quotient_is_exact representative representativeBound]
  unfold exactCrtQuotient
  rw [← exactNumerator]
  have reconstructionNat :
      r2RawLift representative =
        representative +
          (r2RawLift representative / Generated.M2) * Generated.M2 := by
    calc
      r2RawLift representative =
          r2RawLift representative % Generated.M2 +
            Generated.M2 * (r2RawLift representative / Generated.M2) :=
        split.symm
      _ = representative +
          Generated.M2 * (r2RawLift representative / Generated.M2) := by
        rw [remainder]
      _ = representative +
          (r2RawLift representative / Generated.M2) * Generated.M2 := by
        rw [Nat.mul_comm]
  have reconstruction :
      (r2RawLift representative : Int) =
        (representative : Int) +
          ((r2RawLift representative / Generated.M2 : Nat) : Int) *
            (Generated.M2 : Int) := by
    exact_mod_cast reconstructionNat
  rw [reconstruction]
  ring

theorem r2_expansion_lane_represents_input
    (representative : Nat)
    (representativeBound : representative < 18 * Generated.target_modulus)
    (outputModulus : Nat) (outputPositive : 0 < outputModulus) :
    fixedMatrixLaneValue (r2Residues representative) r2ShiftedFactor
        Generated.M2 Generated.r2_postmult outputModulus =
      (Generated.r2_postmult : ZMod outputModulus) * representative := by
  rw [fixed_matrix_lane_represents_corrected_crt _ _ _ _ _ outputPositive]
  rw [r2_corrected_representative_is_exact representative representativeBound]
  norm_cast

theorem generated_r2_entry_formula (row column : Fin 8) :
    matrixEntry Generated.r2_rns_mat row.val column.val =
      fixedMatrixEntry (r2ShiftedFactor row) Generated.r2_postmult
        (listAt8 Generated.moduli1 column.val) := by
  fin_cases row <;> fin_cases column <;> decide

theorem generated_r2_correction_formula (column : Fin 8) :
    listAt8 Generated.r2_correction column.val =
      fixedCorrection Generated.M2 Generated.r2_postmult
        (listAt8 Generated.moduli1 column.val) := by
  fin_cases column <;> decide

theorem generated_r2_quotient_factor_formula (row : Fin 8) :
    listAt8 Generated.r2_shifted_quotient_estimations row.val =
      ceilDiv (r2ShiftedFactor row * uint64Modulus) Generated.M2 := by
  fin_cases row <;> decide

def generatedR2EstimatedQuotient (residues : Fin 8 → Nat) : Nat :=
  (∑ row, residues row *
      listAt8 Generated.r2_shifted_quotient_estimations row.val) /
    uint64Modulus

theorem generated_r2_estimated_quotient_eq
    (residues : Fin 8 → Nat) :
    generatedR2EstimatedQuotient residues =
      estimatedCrtQuotient residues r2ShiftedFactor Generated.M2 := by
  unfold generatedR2EstimatedQuotient estimatedCrtQuotient
  apply congrArg (fun numerator => numerator / uint64Modulus)
  unfold approximateCrtNumerator
  apply Finset.sum_congr rfl
  intro row _
  rw [generated_r2_quotient_factor_formula row]

def generatedR2LaneValue
    (residues : Fin 8 → Nat) (column : Fin 8) :
    ZMod (listAt8 Generated.moduli1 column.val) :=
  (∑ row, (residues row :
      ZMod (listAt8 Generated.moduli1 column.val)) *
      matrixEntry Generated.r2_rns_mat row.val column.val) +
    (generatedR2EstimatedQuotient residues :
      ZMod (listAt8 Generated.moduli1 column.val)) *
      listAt8 Generated.r2_correction column.val

theorem generated_r2_lane_eq_fixed
    (residues : Fin 8 → Nat) (column : Fin 8) :
    generatedR2LaneValue residues column =
      fixedMatrixLaneValue residues r2ShiftedFactor Generated.M2
        Generated.r2_postmult (listAt8 Generated.moduli1 column.val) := by
  unfold generatedR2LaneValue fixedMatrixLaneValue
  rw [generated_r2_estimated_quotient_eq]
  rw [generated_r2_correction_formula]
  congr 1
  apply Finset.sum_congr rfl
  intro row _
  rw [generated_r2_entry_formula row column]

theorem generated_r2_lane_represents_input
    (representative : Nat)
    (representativeBound : representative < 18 * Generated.target_modulus)
    (column : Fin 8) :
    generatedR2LaneValue (r2Residues representative) column =
      (Generated.r2_postmult :
        ZMod (listAt8 Generated.moduli1 column.val)) * representative := by
  rw [generated_r2_lane_eq_fixed]
  apply r2_expansion_lane_represents_input representative representativeBound
  fin_cases column <;> norm_num [Generated.moduli1, listAt8]

/-- One selected MatrixNoK multiplication round: correction-free R1 produces
an `< 18p` representative of the field product, and generated R2 reconstructs
that exact representative in every M1 lane. -/
theorem selected_base_field_round_trip
    (left right : Nat)
    (leftBound : left < 40 * Generated.target_modulus)
    (rightBound : right < 40 * Generated.target_modulus) :
    let representative := r1ActualProductRepresentative left right
    representative < 18 * Generated.target_modulus ∧
      (representative * Generated.M1) % Generated.target_modulus =
        left * right % Generated.target_modulus ∧
      ∀ column : Fin 8,
        generatedR2LaneValue (r2Residues representative) column =
          (Generated.r2_postmult :
            ZMod (listAt8 Generated.moduli1 column.val)) * representative := by
  dsimp only
  exact ⟨r1_actual_product_lt_18p left right leftBound rightBound,
    r1_actual_product_is_field_multiplication left right,
    fun column => generated_r2_lane_represents_input
      (r1ActualProductRepresentative left right)
      (r1_actual_product_lt_18p left right leftBound rightBound) column⟩

end CatBLST.AVX.MatrixNoK
