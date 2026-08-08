/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.NoCorrection
import Mathlib.Data.List.GetD
import Mathlib.Data.Nat.ChineseRemainder
import Mathlib.Tactic

/-!
# Fixed eight-lane CRT reconstruction

This file bridges the checked inverse witnesses to a function-level
reconstruction theorem for the exact M1 and M2 bases.
-/

namespace CatBLST.AVX.MatrixNoK

open scoped BigOperators

def listAt8 (values : List Nat) (index : Nat) : Nat :=
  values.getD index 0

def crtReconstruction8
    (moduli factors : List Nat) (value : Nat) : Nat :=
  let term (index : Nat) :=
    (value % listAt8 moduli index) * listAt8 factors index
  term 0 + term 1 + term 2 + term 3 + term 4 + term 5 + term 6 + term 7

def crtFactors1 : List Nat := crtFactors Generated.moduli1 crtInverses1
def crtFactors2 : List Nat := crtFactors Generated.moduli2 crtInverses2

def scaledCrtReconstruction8
    (moduli factors : List Nat) (premult value : Nat) : Nat :=
  let term (index : Nat) :=
    (value % listAt8 moduli index) *
      (premult * listAt8 factors index % moduli.prod)
  term 0 + term 1 + term 2 + term 3 + term 4 + term 5 + term 6 + term 7

theorem scaled_crt_reconstruction_mod_product
    (moduli factors : List Nat) (premult value : Nat) :
    Nat.ModEq moduli.prod
      (scaledCrtReconstruction8 moduli factors premult value)
      (premult * crtReconstruction8 moduli factors value) := by
  rw [← ZMod.natCast_eq_natCast_iff]
  simp only [scaledCrtReconstruction8, crtReconstruction8,
    Nat.cast_add, Nat.cast_mul, ZMod.natCast_mod]
  ring

theorem crt_factor1_orthogonal
    (row column : Nat) (rowLt : row < 8) (columnLt : column < 8) :
    listAt8 crtFactors1 row % listAt8 Generated.moduli1 column =
      if row = column then 1 else 0 := by
  interval_cases row <;> interval_cases column <;> decide

theorem crt_factor2_orthogonal
    (row column : Nat) (rowLt : row < 8) (columnLt : column < 8) :
    listAt8 crtFactors2 row % listAt8 Generated.moduli2 column =
      if row = column then 1 else 0 := by
  interval_cases row <;> interval_cases column <;> decide

theorem crt_reconstruction_mod_each_generic
    (moduli factors : List Nat) (value index : Nat)
    (indexLt : index < 8)
    (orthogonal : ∀ row, row < 8 →
      listAt8 factors row % listAt8 moduli index =
        if row = index then 1 else 0) :
    Nat.ModEq (listAt8 moduli index)
      (crtReconstruction8 moduli factors value) value := by
  rw [← ZMod.natCast_eq_natCast_iff]
  simp only [crtReconstruction8, Nat.cast_add, Nat.cast_mul]
  have factorCast : ∀ row, row < 8 →
      (listAt8 factors row : ZMod (listAt8 moduli index)) =
        if row = index then 1 else 0 := by
    intro row rowLt
    rw [← ZMod.natCast_mod]
    rw [orthogonal row rowLt]
    split <;> simp
  rw [factorCast 0 (by omega), factorCast 1 (by omega),
    factorCast 2 (by omega), factorCast 3 (by omega),
    factorCast 4 (by omega), factorCast 5 (by omega),
    factorCast 6 (by omega), factorCast 7 (by omega)]
  interval_cases index <;> simp [ZMod.natCast_mod]

theorem crt_reconstruction1_mod_each
    (value index : Nat) (indexLt : index < 8) :
    Nat.ModEq (listAt8 Generated.moduli1 index)
      (crtReconstruction8 Generated.moduli1 crtFactors1 value) value := by
  apply crt_reconstruction_mod_each_generic
    Generated.moduli1 crtFactors1 value index indexLt
  intro row rowLt
  exact crt_factor1_orthogonal row index rowLt indexLt

theorem crt_reconstruction2_mod_each
    (value index : Nat) (indexLt : index < 8) :
    Nat.ModEq (listAt8 Generated.moduli2 index)
      (crtReconstruction8 Generated.moduli2 crtFactors2 value) value := by
  apply crt_reconstruction_mod_each_generic
    Generated.moduli2 crtFactors2 value index indexLt
  intro row rowLt
  exact crt_factor2_orthogonal row index rowLt indexLt

theorem moduli1_pairwise_coprime :
    Generated.moduli1.Pairwise Nat.Coprime := by
  decide

theorem moduli2_pairwise_coprime :
    Generated.moduli2.Pairwise Nat.Coprime := by
  decide

theorem crt_reconstruction1_mod_product (value : Nat) :
    Nat.ModEq Generated.M1
      (crtReconstruction8 Generated.moduli1 crtFactors1 value) value := by
  rw [M1_eq_product]
  apply (Nat.modEq_list_prod_iff moduli1_pairwise_coprime).2
  intro index
  have indexLt : index.val < 8 := by
    have result := index.isLt
    norm_num [Generated.moduli1] at result
    exact result
  have result := crt_reconstruction1_mod_each value index.val indexLt
  simpa [listAt8, List.getD_eq_get] using result

theorem crt_reconstruction2_mod_product (value : Nat) :
    Nat.ModEq Generated.M2
      (crtReconstruction8 Generated.moduli2 crtFactors2 value) value := by
  rw [M2_eq_product]
  apply (Nat.modEq_list_prod_iff moduli2_pairwise_coprime).2
  intro index
  have indexLt : index.val < 8 := by
    have result := index.isLt
    norm_num [Generated.moduli2] at result
    exact result
  have result := crt_reconstruction2_mod_each value index.val indexLt
  simpa [listAt8, List.getD_eq_get] using result

def r1RawLift (value : Nat) : Nat :=
  scaledCrtReconstruction8 Generated.moduli1 crtFactors1
    Generated.r1_premult value

def r1Residues (value : Nat) (index : Fin 8) : Nat :=
  value % listAt8 Generated.moduli1 index.val

def r1ShiftedFactor (index : Fin 8) : Nat :=
  Generated.r1_premult * listAt8 crtFactors1 index.val % Generated.M1

theorem r1_raw_lift_eq_exact_numerator (value : Nat) :
    r1RawLift value = exactCrtNumerator (r1Residues value) r1ShiftedFactor := by
  unfold r1RawLift r1Residues r1ShiftedFactor
  rw [M1_eq_product]
  simp [scaledCrtReconstruction8, exactCrtNumerator, Fin.sum_univ_succ]
  ring

theorem generated_r1_entry_formula (row column : Fin 8) :
    matrixEntry Generated.r1_rns_mat row.val column.val =
      fixedMatrixEntry (r1ShiftedFactor row) Generated.r1_postmult
        (listAt8 Generated.moduli2 column.val) := by
  fin_cases row <;> fin_cases column <;> decide

def generatedR1LaneValue
    (residues : Fin 8 → Nat) (column : Fin 8) :
    ZMod (listAt8 Generated.moduli2 column.val) :=
  ∑ row, (residues row : ZMod (listAt8 Generated.moduli2 column.val)) *
    matrixEntry Generated.r1_rns_mat row.val column.val

theorem generated_r1_lane_represents_raw_crt
    (residues : Fin 8 → Nat) (column : Fin 8) :
    generatedR1LaneValue residues column =
      (Generated.r1_postmult :
        ZMod (listAt8 Generated.moduli2 column.val)) *
        exactCrtNumerator residues r1ShiftedFactor := by
  unfold generatedR1LaneValue
  have raw := raw_fixed_matrix_lane_represents_crt
    residues r1ShiftedFactor Generated.r1_postmult
      (listAt8 Generated.moduli2 column.val)
  rw [← raw]
  unfold rawFixedMatrixLaneValue
  apply Finset.sum_congr rfl
  intro row _
  rw [generated_r1_entry_formula row column]

theorem r1_quotient_estimate_is_fixed (index : Fin 8) :
    ceilDiv (r1ShiftedFactor index * uint64Modulus) Generated.M1 = 16385 := by
  fin_cases index <;> decide

theorem r1_residue_lt_element_radix (value : Nat) (index : Fin 8) :
    r1Residues value index < elementRadix := by
  fin_cases index
  all_goals
    unfold r1Residues listAt8
    exact (Nat.mod_lt value (by norm_num [Generated.moduli1])).trans
      (by norm_num [Generated.moduli1, elementRadix])

theorem r1_approximate_numerator_formula (value : Nat) :
    approximateCrtNumerator (r1Residues value) r1ShiftedFactor Generated.M1 =
      (∑ index, r1Residues value index) * 16385 := by
  unfold approximateCrtNumerator
  rw [Finset.sum_mul]
  apply Finset.sum_congr rfl
  intro index _
  rw [r1_quotient_estimate_is_fixed]

theorem r1_estimated_quotient_lt_17 (value : Nat) :
    estimatedCrtQuotient (r1Residues value) r1ShiftedFactor Generated.M1 < 17 := by
  have sumBound : (∑ index, r1Residues value index) < 8 * elementRadix :=
    fin8_sum_lt (r1_residue_lt_element_radix value)
  have approximateBound :
      approximateCrtNumerator (r1Residues value) r1ShiftedFactor Generated.M1 <
        17 * uint64Modulus := by
    rw [r1_approximate_numerator_formula]
    calc
      (∑ index, r1Residues value index) * 16385 <
          (8 * elementRadix) * 16385 :=
        by nlinarith
      _ < 17 * uint64Modulus := by
        norm_num [elementRadix, uint64Modulus]
  unfold estimatedCrtQuotient
  exact (Nat.div_lt_iff_lt_mul (by norm_num [uint64Modulus])).2
    (by simpa [Nat.mul_comm] using approximateBound)

theorem r1_exact_quotient_lt_17 (value : Nat) :
    exactCrtQuotient (r1Residues value) r1ShiftedFactor Generated.M1 < 17 := by
  have scalarFits : ∀ index, r1Residues value index < ifmaRadix := by
    intro index
    exact (r1_residue_lt_element_radix value index).trans
      (by norm_num [elementRadix, ifmaRadix])
  have bracket := estimated_quotient_bracket
    (r1Residues value) r1ShiftedFactor Generated.M1
    (by norm_num [Generated.M1]) scalarFits
  exact bracket.1.trans_lt (r1_estimated_quotient_lt_17 value)

theorem r1_raw_lift_lt_17_product (value : Nat) :
    r1RawLift value < 17 * Generated.M1 := by
  rw [r1_raw_lift_eq_exact_numerator]
  have split := Nat.mod_add_div
    (exactCrtNumerator (r1Residues value) r1ShiftedFactor) Generated.M1
  have modulusPositive : 0 < Generated.M1 := by
    norm_num [Generated.M1]
  have remainder := Nat.mod_lt
    (exactCrtNumerator (r1Residues value) r1ShiftedFactor)
    modulusPositive
  have quotient := r1_exact_quotient_lt_17 value
  unfold exactCrtQuotient at quotient
  nlinarith

theorem r1_raw_lift_mod_product (value : Nat) :
    r1RawLift value % Generated.M1 =
      (Generated.r1_premult * value) % Generated.M1 := by
  have scaled := scaled_crt_reconstruction_mod_product
    Generated.moduli1 crtFactors1 Generated.r1_premult value
  have reconstruction := crt_reconstruction1_mod_product value
  have multiplied :=
    (Nat.ModEq.refl (n := Generated.M1) Generated.r1_premult).mul
      reconstruction
  rw [M1_eq_product] at multiplied
  unfold r1RawLift
  rw [M1_eq_product]
  exact scaled.trans multiplied

theorem r1_raw_lift_matches_canonical (left right : Nat) :
    r1RawLift
        (left * right * montgomeryRadix * Generated.qr_premult *
          Generated.qr_premult % Generated.M1) % Generated.M1 =
      r1CanonicalLift left right := by
  rw [r1_raw_lift_mod_product]
  rfl

def r1ReducedProduct (left right : Nat) : Nat :=
  left * right * montgomeryRadix * Generated.qr_premult *
    Generated.qr_premult % Generated.M1

def r1ActualProductRepresentative (left right : Nat) : Nat :=
  (left * right +
      Generated.target_modulus * r1RawLift (r1ReducedProduct left right)) /
    Generated.M1

theorem r1_actual_numerator_divisible (left right : Nat) :
    Generated.M1 ∣
      left * right +
        Generated.target_modulus * r1RawLift (r1ReducedProduct left right) := by
  apply r1_raw_numerator_divisible
  unfold r1ReducedProduct
  exact r1_raw_lift_matches_canonical left right

theorem r1_actual_product_reconstruction (left right : Nat) :
    r1ActualProductRepresentative left right * Generated.M1 =
      left * right +
        Generated.target_modulus * r1RawLift (r1ReducedProduct left right) := by
  unfold r1ActualProductRepresentative
  exact Nat.div_mul_cancel (r1_actual_numerator_divisible left right)

theorem r1_actual_product_lt_18p
    (left right : Nat)
    (leftBound : left < 40 * Generated.target_modulus)
    (rightBound : right < 40 * Generated.target_modulus) :
    r1ActualProductRepresentative left right <
      18 * Generated.target_modulus := by
  have targetPositive : 0 < Generated.target_modulus := by
    norm_num [Generated.target_modulus]
  have modulusPositive : 0 < Generated.M1 := by
    norm_num [Generated.M1]
  have dynamicRange :
      1600 * Generated.target_modulus < Generated.M1 := by
    norm_num [Generated.target_modulus, Generated.M1]
  have productBound :
      left * right < Generated.target_modulus * Generated.M1 := by
    nlinarith
  have liftBound := r1_raw_lift_lt_17_product (r1ReducedProduct left right)
  have numeratorBound :
      left * right +
          Generated.target_modulus * r1RawLift (r1ReducedProduct left right) <
        (18 * Generated.target_modulus) * Generated.M1 := by
    nlinarith
  unfold r1ActualProductRepresentative
  exact (Nat.div_lt_iff_lt_mul modulusPositive).2
    (by simpa [Nat.mul_assoc, Nat.mul_comm, Nat.mul_left_comm] using numeratorBound)

theorem r1_actual_product_is_field_multiplication (left right : Nat) :
    (r1ActualProductRepresentative left right * Generated.M1) %
        Generated.target_modulus =
      left * right % Generated.target_modulus := by
  rw [r1_actual_product_reconstruction]
  simp [Nat.add_mod]

theorem r1_actual_product_cross_basis_output (left right : Nat) :
    ((left * right + Generated.target_modulus *
        r1RawLift (r1ReducedProduct left right)) *
      montgomeryRadix * crossBasisInverse * crossBasisInverse) %
        Generated.M2 =
      (r1ActualProductRepresentative left right * montgomeryRadix *
        crossBasisInverse) % Generated.M2 := by
  exact cross_basis_output_generic
    Generated.M1 Generated.M2 crossBasisInverse montgomeryRadix
    (left * right + Generated.target_modulus *
      r1RawLift (r1ReducedProduct left right))
    (r1ActualProductRepresentative left right)
    (by norm_num [Generated.M2])
    (r1_actual_product_reconstruction left right)
    cross_basis_inverse_relation

end CatBLST.AVX.MatrixNoK
