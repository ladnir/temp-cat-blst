/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.MatrixBounds
import Mathlib.Tactic

/-!
# Selected MatrixNoK base-field pipeline

This file lifts the lane and matrix theorems to the actual standard-residue
contract.  In particular, R1 consumes residues below `2m_i`; it does not assume
that the preceding Montgomery reducer happened to choose the canonical
representative below `m_i`.
-/

namespace CatBLST.AVX.MatrixNoK

open scoped BigOperators

def crtFactor1At (index : Fin 8) : Nat :=
  listAt8 crtFactors1 index.val

def r1RawLiftFromResidues (residues : Fin 8 → Nat) : Nat :=
  exactCrtNumerator residues r1ShiftedFactor

theorem exact_crt1_residues_mod_each
    (residues : Fin 8 → Nat) (column : Fin 8) :
    Nat.ModEq (listAt8 Generated.moduli1 column.val)
      (exactCrtNumerator residues crtFactor1At) (residues column) := by
  rw [← ZMod.natCast_eq_natCast_iff]
  unfold exactCrtNumerator crtFactor1At
  push_cast
  have factorCast : ∀ row : Fin 8,
      (listAt8 crtFactors1 row.val :
          ZMod (listAt8 Generated.moduli1 column.val)) =
        if row = column then 1 else 0 := by
    intro row
    rw [← ZMod.natCast_mod]
    rw [crt_factor1_orthogonal row.val column.val row.isLt column.isLt]
    by_cases equal : row = column
    · subst row
      simp
    · have valuesDiffer : row.val ≠ column.val := by
        intro valuesEqual
        exact equal (Fin.ext valuesEqual)
      simp [equal, valuesDiffer]
  simp_rw [factorCast]
  simp

theorem exact_crt1_residues_mod_product
    (residues : Fin 8 → Nat) (value : Nat)
    (residuesRepresent : ∀ index,
      Nat.ModEq (listAt8 Generated.moduli1 index.val)
        (residues index) value) :
    Nat.ModEq Generated.M1
      (exactCrtNumerator residues crtFactor1At) value := by
  rw [M1_eq_product]
  apply (Nat.modEq_list_prod_iff moduli1_pairwise_coprime).2
  intro index
  have indexLt : index.val < 8 := by
    have result := index.isLt
    norm_num [Generated.moduli1] at result
    exact result
  let lane : Fin 8 := ⟨index.val, indexLt⟩
  have reconstructed := exact_crt1_residues_mod_each residues lane
  have represented := residuesRepresent lane
  simpa [lane, listAt8, List.getD_eq_get] using
    reconstructed.trans represented

theorem r1_raw_lift_from_residues_scaled
    (residues : Fin 8 → Nat) :
    Nat.ModEq Generated.M1 (r1RawLiftFromResidues residues)
      (Generated.r1_premult *
        exactCrtNumerator residues crtFactor1At) := by
  rw [← ZMod.natCast_eq_natCast_iff]
  unfold r1RawLiftFromResidues exactCrtNumerator
  push_cast
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro index _
  unfold r1ShiftedFactor crtFactor1At
  rw [ZMod.natCast_mod]
  push_cast
  ring

theorem r1_raw_lift_from_residues_relation
    (residues : Fin 8 → Nat) (value : Nat)
    (residuesRepresent : ∀ index,
      Nat.ModEq (listAt8 Generated.moduli1 index.val)
        (residues index) value) :
    r1RawLiftFromResidues residues % Generated.M1 =
      (Generated.r1_premult * value) % Generated.M1 := by
  have reconstructed := exact_crt1_residues_mod_product
    residues value residuesRepresent
  have scaled :=
    (Nat.ModEq.refl (n := Generated.M1) Generated.r1_premult).mul
      reconstructed
  exact (r1_raw_lift_from_residues_scaled residues).trans scaled

theorem r1_standard_residue_lt_ifma
    (residues : Fin 8 → Nat)
    (residuesFit : ∀ index,
      residues index < 2 * listAt8 Generated.moduli1 index.val) :
    ∀ index, residues index < ifmaRadix := by
  intro index
  have modulusFit :
      2 * listAt8 Generated.moduli1 index.val < ifmaRadix := by
    fin_cases index <;> norm_num [Generated.moduli1, listAt8, ifmaRadix]
  exact (residuesFit index).trans modulusFit

theorem r1_standard_approximate_numerator_formula
    (residues : Fin 8 → Nat) :
    approximateCrtNumerator residues r1ShiftedFactor Generated.M1 =
      (∑ index, residues index) * 16385 := by
  unfold approximateCrtNumerator
  rw [Finset.sum_mul]
  apply Finset.sum_congr rfl
  intro index _
  rw [r1_quotient_estimate_is_fixed]

theorem r1_standard_estimated_quotient_lt_17
    (residues : Fin 8 → Nat)
    (residuesFit : ∀ index,
      residues index < 2 * listAt8 Generated.moduli1 index.val) :
    estimatedCrtQuotient residues r1ShiftedFactor Generated.M1 < 17 := by
  have eachFits : ∀ index, residues index < 2 * elementRadix := by
    intro index
    have modulusBelow :
        listAt8 Generated.moduli1 index.val < elementRadix := by
      fin_cases index <;>
        norm_num [Generated.moduli1, listAt8, elementRadix]
    have doubledBelow :
        2 * listAt8 Generated.moduli1 index.val < 2 * elementRadix := by
      nlinarith
    exact (residuesFit index).trans doubledBelow
  have sumBound := fin8_sum_lt eachFits
  have approximateBound :
      approximateCrtNumerator residues r1ShiftedFactor Generated.M1 <
        17 * uint64Modulus := by
    rw [r1_standard_approximate_numerator_formula]
    calc
      (∑ index, residues index) * 16385 <
          (8 * (2 * elementRadix)) * 16385 :=
        (Nat.mul_lt_mul_right (by norm_num)).2 sumBound
      _ < 17 * uint64Modulus := by
        norm_num [elementRadix, uint64Modulus]
  unfold estimatedCrtQuotient
  exact (Nat.div_lt_iff_lt_mul (by norm_num [uint64Modulus])).2
    (by simpa [Nat.mul_comm] using approximateBound)

theorem r1_standard_raw_lift_lt_17_product
    (residues : Fin 8 → Nat)
    (residuesFit : ∀ index,
      residues index < 2 * listAt8 Generated.moduli1 index.val) :
    r1RawLiftFromResidues residues < 17 * Generated.M1 := by
  have scalarFits := r1_standard_residue_lt_ifma residues residuesFit
  have bracket := estimated_quotient_bracket
    residues r1ShiftedFactor Generated.M1
    (by norm_num [Generated.M1]) scalarFits
  have exactQuotient :
      exactCrtQuotient residues r1ShiftedFactor Generated.M1 < 17 :=
    bracket.1.trans_lt
      (r1_standard_estimated_quotient_lt_17 residues residuesFit)
  have split := Nat.mod_add_div
    (exactCrtNumerator residues r1ShiftedFactor) Generated.M1
  have modulusPositive : 0 < Generated.M1 := by norm_num [Generated.M1]
  have remainder := Nat.mod_lt
    (exactCrtNumerator residues r1ShiftedFactor)
    modulusPositive
  unfold exactCrtQuotient at exactQuotient
  unfold r1RawLiftFromResidues
  nlinarith

def r1ProductFromResidues
    (left right : Nat) (residues : Fin 8 → Nat) : Nat :=
  (left * right +
      Generated.target_modulus * r1RawLiftFromResidues residues) /
    Generated.M1

theorem r1_product_from_standard_residues
    (left right : Nat) (residues : Fin 8 → Nat)
    (residuesRepresent : ∀ index,
      Nat.ModEq (listAt8 Generated.moduli1 index.val)
        (residues index) (r1ReducedProduct left right)) :
    let representative := r1ProductFromResidues left right residues
    representative * Generated.M1 =
        left * right + Generated.target_modulus *
          r1RawLiftFromResidues residues ∧
      (representative * Generated.M1) % Generated.target_modulus =
        left * right % Generated.target_modulus := by
  dsimp only
  have rawRelation := r1_raw_lift_from_residues_relation
    residues (r1ReducedProduct left right) residuesRepresent
  have canonicalRelation :
      r1RawLiftFromResidues residues % Generated.M1 =
        r1CanonicalLift left right := by
    rw [rawRelation]
    unfold r1ReducedProduct r1CanonicalLift
    rw [Nat.mul_mod]
  have divisible := r1_raw_numerator_divisible left right
    (r1RawLiftFromResidues residues) canonicalRelation
  have reconstruction :
      r1ProductFromResidues left right residues * Generated.M1 =
        left * right + Generated.target_modulus *
          r1RawLiftFromResidues residues := by
    unfold r1ProductFromResidues
    exact Nat.div_mul_cancel divisible
  refine ⟨reconstruction, ?_⟩
  rw [reconstruction]
  simp [Nat.add_mod]

theorem r1_product_from_standard_residues_lt_18p
    (left right : Nat) (residues : Fin 8 → Nat)
    (leftBound : left < 40 * Generated.target_modulus)
    (rightBound : right < 40 * Generated.target_modulus)
    (residuesFit : ∀ index,
      residues index < 2 * listAt8 Generated.moduli1 index.val) :
    r1ProductFromResidues left right residues <
      18 * Generated.target_modulus := by
  have targetPositive : 0 < Generated.target_modulus := by
    norm_num [Generated.target_modulus]
  have modulusPositive : 0 < Generated.M1 := by norm_num [Generated.M1]
  have dynamicRange : 1600 * Generated.target_modulus < Generated.M1 := by
    norm_num [Generated.target_modulus, Generated.M1]
  have productBound :
      left * right < Generated.target_modulus * Generated.M1 := by
    nlinarith
  have liftBound := r1_standard_raw_lift_lt_17_product residues residuesFit
  unfold r1ProductFromResidues
  apply (Nat.div_lt_iff_lt_mul modulusPositive).2
  nlinarith

def crtFactor2At (index : Fin 8) : Nat :=
  listAt8 crtFactors2 index.val

def r2RawLiftFromResidues (residues : Fin 8 → Nat) : Nat :=
  exactCrtNumerator residues r2ShiftedFactor

theorem exact_crt2_residues_mod_each
    (residues : Fin 8 → Nat) (column : Fin 8) :
    Nat.ModEq (listAt8 Generated.moduli2 column.val)
      (exactCrtNumerator residues crtFactor2At) (residues column) := by
  rw [← ZMod.natCast_eq_natCast_iff]
  unfold exactCrtNumerator crtFactor2At
  push_cast
  have factorCast : ∀ row : Fin 8,
      (listAt8 crtFactors2 row.val :
          ZMod (listAt8 Generated.moduli2 column.val)) =
        if row = column then 1 else 0 := by
    intro row
    rw [← ZMod.natCast_mod]
    rw [crt_factor2_orthogonal row.val column.val row.isLt column.isLt]
    by_cases equal : row = column
    · subst row
      simp
    · have valuesDiffer : row.val ≠ column.val := by
        intro valuesEqual
        exact equal (Fin.ext valuesEqual)
      simp [equal, valuesDiffer]
  simp_rw [factorCast]
  simp

theorem exact_crt2_residues_mod_product
    (residues : Fin 8 → Nat) (value : Nat)
    (residuesRepresent : ∀ index,
      Nat.ModEq (listAt8 Generated.moduli2 index.val)
        (residues index) value) :
    Nat.ModEq Generated.M2
      (exactCrtNumerator residues crtFactor2At) value := by
  rw [M2_eq_product]
  apply (Nat.modEq_list_prod_iff moduli2_pairwise_coprime).2
  intro index
  have indexLt : index.val < 8 := by
    have result := index.isLt
    norm_num [Generated.moduli2] at result
    exact result
  let lane : Fin 8 := ⟨index.val, indexLt⟩
  have reconstructed := exact_crt2_residues_mod_each residues lane
  have represented := residuesRepresent lane
  simpa [lane, listAt8, List.getD_eq_get] using
    reconstructed.trans represented

theorem r2_raw_lift_from_residues_scaled
    (residues : Fin 8 → Nat) :
    Nat.ModEq Generated.M2 (r2RawLiftFromResidues residues)
      (Generated.r2_premult *
        exactCrtNumerator residues crtFactor2At) := by
  rw [← ZMod.natCast_eq_natCast_iff]
  unfold r2RawLiftFromResidues exactCrtNumerator
  push_cast
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro index _
  unfold r2ShiftedFactor crtFactor2At
  rw [ZMod.natCast_mod]
  push_cast
  ring

theorem r2_raw_lift_from_residues_remainder
    (representative : Nat) (residues : Fin 8 → Nat)
    (representativeBound : representative < Generated.M2)
    (residuesRepresent : ∀ index,
      Nat.ModEq (listAt8 Generated.moduli2 index.val)
        (residues index) (r2InputValue representative)) :
    r2RawLiftFromResidues residues % Generated.M2 = representative := by
  have reconstructed := exact_crt2_residues_mod_product residues
    (r2InputValue representative) residuesRepresent
  have scaled :=
    (Nat.ModEq.refl (n := Generated.M2) Generated.r2_premult).mul
      reconstructed
  have rawRelation :=
    (r2_raw_lift_from_residues_scaled residues).trans scaled
  have scaleRelation : Nat.ModEq Generated.M2
      (Generated.r2_premult * r2InputValue representative) representative :=
    r2_input_scale_relation representative
  simpa [Nat.ModEq, Nat.mod_eq_of_lt representativeBound] using
    rawRelation.trans scaleRelation

theorem r2_standard_residue_lt_ifma
    (residues : Fin 8 → Nat)
    (residuesFit : ∀ index,
      residues index < 2 * listAt8 Generated.moduli2 index.val) :
    ∀ index, residues index < ifmaRadix := by
  intro index
  have modulusFit :
      2 * listAt8 Generated.moduli2 index.val < ifmaRadix := by
    fin_cases index <;> norm_num [Generated.moduli2, listAt8, ifmaRadix]
  exact (residuesFit index).trans modulusFit

theorem r2_standard_scalar_sum_bound
    (residues : Fin 8 → Nat)
    (residuesFit : ∀ index,
      residues index < 2 * listAt8 Generated.moduli2 index.val) :
    (∑ index, residues index) < 16 * elementRadix := by
  have eachFits : ∀ index, residues index < 2 * elementRadix := by
    intro index
    have modulusBelow :
        listAt8 Generated.moduli2 index.val < elementRadix := by
      fin_cases index <;>
        norm_num [Generated.moduli2, listAt8, elementRadix]
    have doubledBelow :
        2 * listAt8 Generated.moduli2 index.val < 2 * elementRadix := by
      nlinarith
    exact (residuesFit index).trans doubledBelow
  simpa [Nat.mul_assoc, Nat.mul_comm, Nat.mul_left_comm] using
    (fin8_sum_lt eachFits)

theorem r2_standard_positive_margin
    (representative : Nat) (residues : Fin 8 → Nat)
    (representativeBound : representative < 18 * Generated.target_modulus)
    (residuesFit : ∀ index,
      residues index < 2 * listAt8 Generated.moduli2 index.val) :
    representative * uint64Modulus +
        (Generated.M2 - 1) * ∑ index, residues index <
      Generated.M2 * uint64Modulus := by
  have sumBound := r2_standard_scalar_sum_bound residues residuesFit
  have envelope :
      (18 * Generated.target_modulus) * uint64Modulus +
          (Generated.M2 - 1) * (16 * elementRadix) <
        Generated.M2 * uint64Modulus := by
    norm_num [Generated.target_modulus, Generated.M2,
      uint64Modulus, elementRadix]
  nlinarith

theorem r2_standard_estimated_quotient_is_exact
    (representative : Nat) (residues : Fin 8 → Nat)
    (representativeBound : representative < 18 * Generated.target_modulus)
    (residuesFit : ∀ index,
      residues index < 2 * listAt8 Generated.moduli2 index.val)
    (residuesRepresent : ∀ index,
      Nat.ModEq (listAt8 Generated.moduli2 index.val)
        (residues index) (r2InputValue representative)) :
    estimatedCrtQuotient residues r2ShiftedFactor Generated.M2 =
      exactCrtQuotient residues r2ShiftedFactor Generated.M2 := by
  have belowM2 : representative < Generated.M2 :=
    representativeBound.trans (by
      norm_num [Generated.target_modulus, Generated.M2])
  have remainder := r2_raw_lift_from_residues_remainder representative
    residues belowM2 residuesRepresent
  have scalarFits := r2_standard_residue_lt_ifma residues residuesFit
  apply estimated_quotient_exact_of_margin residues r2ShiftedFactor
    Generated.M2 (by norm_num [Generated.M2]) scalarFits
  unfold r2RawLiftFromResidues at remainder
  rw [remainder]
  exact r2_standard_positive_margin representative residues
    representativeBound residuesFit

theorem r2_standard_corrected_representative_is_exact
    (representative : Nat) (residues : Fin 8 → Nat)
    (representativeBound : representative < 18 * Generated.target_modulus)
    (residuesFit : ∀ index,
      residues index < 2 * listAt8 Generated.moduli2 index.val)
    (residuesRepresent : ∀ index,
      Nat.ModEq (listAt8 Generated.moduli2 index.val)
        (residues index) (r2InputValue representative)) :
    correctedCrtRepresentative residues r2ShiftedFactor Generated.M2 =
      representative := by
  have belowM2 : representative < Generated.M2 :=
    representativeBound.trans (by
      norm_num [Generated.target_modulus, Generated.M2])
  have remainder := r2_raw_lift_from_residues_remainder representative
    residues belowM2 residuesRepresent
  have split := Nat.mod_add_div
    (r2RawLiftFromResidues residues) Generated.M2
  unfold correctedCrtRepresentative
  rw [r2_standard_estimated_quotient_is_exact representative residues
    representativeBound residuesFit residuesRepresent]
  unfold exactCrtQuotient r2RawLiftFromResidues at *
  have reconstructionNat :
      exactCrtNumerator residues r2ShiftedFactor =
        representative +
          (exactCrtNumerator residues r2ShiftedFactor / Generated.M2) *
            Generated.M2 := by
    calc
      exactCrtNumerator residues r2ShiftedFactor =
          exactCrtNumerator residues r2ShiftedFactor % Generated.M2 +
            Generated.M2 *
              (exactCrtNumerator residues r2ShiftedFactor / Generated.M2) :=
        split.symm
      _ = representative + Generated.M2 *
          (exactCrtNumerator residues r2ShiftedFactor / Generated.M2) := by
        rw [remainder]
      _ = representative +
          (exactCrtNumerator residues r2ShiftedFactor / Generated.M2) *
            Generated.M2 := by rw [Nat.mul_comm]
  have reconstructionInt :
      (exactCrtNumerator residues r2ShiftedFactor : Int) =
        (representative : Int) +
          ((exactCrtNumerator residues r2ShiftedFactor / Generated.M2 : Nat) : Int) *
            (Generated.M2 : Int) := by
    exact_mod_cast reconstructionNat
  rw [reconstructionInt]
  ring

theorem generated_r2_standard_lane_represents_input
    (representative : Nat) (residues : Fin 8 → Nat)
    (representativeBound : representative < 18 * Generated.target_modulus)
    (residuesFit : ∀ index,
      residues index < 2 * listAt8 Generated.moduli2 index.val)
    (residuesRepresent : ∀ index,
      Nat.ModEq (listAt8 Generated.moduli2 index.val)
        (residues index) (r2InputValue representative))
    (column : Fin 8) :
    generatedR2LaneValue residues column =
      (Generated.r2_postmult :
        ZMod (listAt8 Generated.moduli1 column.val)) * representative := by
  rw [generated_r2_lane_eq_fixed]
  rw [fixed_matrix_lane_represents_corrected_crt]
  · rw [r2_standard_corrected_representative_is_exact representative residues
      representativeBound residuesFit residuesRepresent]
    norm_cast
  · fin_cases column <;> norm_num [Generated.moduli1, listAt8]

def modulus1At (index : Fin 8) : Nat :=
  listAt8 Generated.moduli1 index.val

def modulus2At (index : Fin 8) : Nat :=
  listAt8 Generated.moduli2 index.val

def reducer1FactorAt (index : Fin 8) : Nat :=
  listAt8 Generated.reducer1_mont_factor index.val

def reducer2FactorAt (index : Fin 8) : Nat :=
  listAt8 Generated.reducer2_mont_factor index.val

theorem modulus1_member (index : Fin 8) : modulus1At index ∈ allModuli := by
  fin_cases index <;>
    simp [modulus1At, listAt8, allModuli, Generated.moduli1,
      Generated.moduli2]

theorem modulus2_member (index : Fin 8) : modulus2At index ∈ allModuli := by
  fin_cases index <;>
    simp [modulus2At, listAt8, allModuli, Generated.moduli1,
      Generated.moduli2]

theorem reducer1_factor_relation (index : Fin 8) :
    modulus1At index * reducer1FactorAt index % ifmaRadix = 1 := by
  fin_cases index <;> decide

theorem reducer2_factor_relation (index : Fin 8) :
    modulus2At index * reducer2FactorAt index % ifmaRadix = 1 := by
  fin_cases index <;> decide

theorem modulus1_dvd_product (index : Fin 8) :
    modulus1At index ∣ Generated.M1 := by
  fin_cases index <;>
    norm_num [modulus1At, listAt8, Generated.moduli1, Generated.M1]

theorem modulus2_dvd_product (index : Fin 8) :
    modulus2At index ∣ Generated.M2 := by
  fin_cases index <;>
    norm_num [modulus2At, listAt8, Generated.moduli2, Generated.M2]

theorem ifma_radix_is_unit_modulus1 (index : Fin 8) :
    IsUnit (ifmaRadix : ZMod (modulus1At index)) := by
  rw [ZMod.isUnit_iff_coprime]
  unfold ifmaRadix
  rw [Nat.coprime_pow_left_iff (by norm_num)]
  exact (Nat.coprime_two_left).2
    (all_moduli_odd (modulus1At index) (modulus1_member index))

theorem ifma_radix_is_unit_modulus2 (index : Fin 8) :
    IsUnit (ifmaRadix : ZMod (modulus2At index)) := by
  rw [ZMod.isUnit_iff_coprime]
  unfold ifmaRadix
  rw [Nat.coprime_pow_left_iff (by norm_num)]
  exact (Nat.coprime_two_left).2
    (all_moduli_odd (modulus2At index) (modulus2_member index))

def productReducedLane
    (left right modulus factor : Nat) : Nat :=
  let product := left * right
  selectedMontgomeryStep (product / ifmaRadix) (product % ifmaRadix)
    modulus factor

theorem product_reduced_lane_contract
    (left right modulus factor : Nat)
    (modulusMember : modulus ∈ allModuli)
    (inverseRelation : modulus * factor % ifmaRadix = 1)
    (leftFit : left < 2 * modulus) (rightFit : right < 2 * modulus) :
    (productReducedLane left right modulus factor : ZMod modulus) *
        ifmaRadix = left * right ∧
      productReducedLane left right modulus factor < 2 * modulus := by
  have productBound : left * right ≤ 4 * modulus * modulus := by
    nlinarith
  have split := Nat.mod_add_div (left * right) ifmaRadix
  have splitReordered :
      (left * right / ifmaRadix) * ifmaRadix +
          left * right % ifmaRadix = left * right := by
    simpa [Nat.add_comm, Nat.mul_comm] using split
  have wideFits :
      (left * right / ifmaRadix) * ifmaRadix +
          left * right % ifmaRadix ≤ 4 * modulus * modulus := by
    rw [splitReordered]
    exact productBound
  have congruent := selected_step_congruent_wide
    (left * right / ifmaRadix) (left * right % ifmaRadix)
    modulus factor
    (selected_modulus_quarter_window modulus modulusMember).1 inverseRelation
  have congruentCast :
      (productReducedLane left right modulus factor : ZMod modulus) *
          ifmaRadix = left * right := by
    unfold productReducedLane
    have castEquality :=
      (ZMod.natCast_eq_natCast_iff
        (selectedMontgomeryStep (left * right / ifmaRadix)
          (left * right % ifmaRadix) modulus factor * ifmaRadix)
        ((left * right / ifmaRadix) * ifmaRadix +
          left * right % ifmaRadix) modulus).2 congruent
    push_cast at castEquality
    have splitCast := congrArg
      (fun value : Nat => (value : ZMod modulus)) splitReordered
    push_cast at splitCast
    rw [splitCast] at castEquality
    simpa [productReducedLane] using castEquality
  refine ⟨congruentCast, ?_⟩
  unfold productReducedLane
  exact selected_product_step_bound
    (left * right / ifmaRadix) (left * right % ifmaRadix)
    modulus factor modulusMember wideFits

def m1ProductResidues
    (left right : Fin 8 → Nat) (index : Fin 8) : Nat :=
  productReducedLane (left index) (right index)
    (modulus1At index) (reducer1FactorAt index)

theorem m1_product_residues_fit
    (left right : Fin 8 → Nat)
    (leftFit : ∀ index, left index < 2 * modulus1At index)
    (rightFit : ∀ index, right index < 2 * modulus1At index) :
    ∀ index, m1ProductResidues left right index < 2 * modulus1At index := by
  intro index
  exact (product_reduced_lane_contract
    (left index) (right index) (modulus1At index) (reducer1FactorAt index)
    (modulus1_member index) (reducer1_factor_relation index)
    (leftFit index) (rightFit index)).2

theorem m1_product_residues_represent
    (leftValue rightValue : Nat) (left right : Fin 8 → Nat)
    (leftRepresent : ∀ index,
      Nat.ModEq (modulus1At index) (left index)
        (leftValue * ifmaRadix * Generated.qr_premult))
    (rightRepresent : ∀ index,
      Nat.ModEq (modulus1At index) (right index)
        (rightValue * ifmaRadix * Generated.qr_premult))
    (leftFit : ∀ index, left index < 2 * modulus1At index)
    (rightFit : ∀ index, right index < 2 * modulus1At index) :
    ∀ index, Nat.ModEq (modulus1At index)
      (m1ProductResidues left right index)
      (r1ReducedProduct leftValue rightValue) := by
  intro index
  let modulus := modulus1At index
  have reducedContract := (product_reduced_lane_contract
    (left index) (right index) modulus (reducer1FactorAt index)
    (modulus1_member index) (reducer1_factor_relation index)
    (leftFit index) (rightFit index)).1
  have leftCast := (ZMod.natCast_eq_natCast_iff
    (left index) (leftValue * ifmaRadix * Generated.qr_premult) modulus).2
      (leftRepresent index)
  have rightCast := (ZMod.natCast_eq_natCast_iff
    (right index) (rightValue * ifmaRadix * Generated.qr_premult) modulus).2
      (rightRepresent index)
  have reducedTarget : Nat.ModEq modulus
      (r1ReducedProduct leftValue rightValue)
      (leftValue * rightValue * ifmaRadix * Generated.qr_premult *
        Generated.qr_premult) := by
    unfold r1ReducedProduct
    change
      (leftValue * rightValue * ifmaRadix * Generated.qr_premult *
          Generated.qr_premult % Generated.M1) % modulus =
        (leftValue * rightValue * ifmaRadix * Generated.qr_premult *
          Generated.qr_premult) % modulus
    exact Nat.mod_mod_of_dvd _ (modulus1_dvd_product index)
  have targetCast := (ZMod.natCast_eq_natCast_iff _ _ modulus).2 reducedTarget
  have targetTimesRadix :
      (r1ReducedProduct leftValue rightValue : ZMod modulus) * ifmaRadix =
        leftValue * rightValue * ifmaRadix * ifmaRadix *
          Generated.qr_premult * Generated.qr_premult := by
    rw [targetCast]
    push_cast
    ring_nf
  have equalTimesRadix :
      (m1ProductResidues left right index : ZMod modulus) * ifmaRadix =
        (r1ReducedProduct leftValue rightValue : ZMod modulus) * ifmaRadix := by
    unfold m1ProductResidues
    rw [reducedContract, leftCast, rightCast, targetTimesRadix]
    push_cast
    ring
  have equalCast :=
    (ifma_radix_is_unit_modulus1 index).mul_right_cancel equalTimesRadix
  exact (ZMod.natCast_eq_natCast_iff _ _ modulus).1 equalCast

def combinedM2WideValue
    (left right m1Residues : Fin 8 → Nat) (column : Fin 8) : Nat :=
  left column * right column +
    matrixProductValue m1Residues Generated.r1_rns_mat column

def m2CombinedResidues
    (left right m1Residues : Fin 8 → Nat) (column : Fin 8) : Nat :=
  let wide := combinedM2WideValue left right m1Residues column
  selectedReduceMasked (wide / ifmaRadix) (wide % ifmaRadix)
    (modulus2At column) (reducer2FactorAt column)

theorem combined_m2_wide_lt_20_square
    (left right m1Residues : Fin 8 → Nat)
    (leftFit : ∀ index, left index < 2 * modulus2At index)
    (rightFit : ∀ index, right index < 2 * modulus2At index)
    (m1Fit : ∀ index, m1Residues index < 2 * modulus1At index)
    (column : Fin 8) :
    combinedM2WideValue left right m1Residues column <
      20 * (modulus2At column) ^ 2 := by
  have productBound :
      left column * right column < 4 * (modulus2At column) ^ 2 := by
    have leftBound := leftFit column
    have rightBound := rightFit column
    nlinarith
  have matrixBound := r1_matrix_product_lt_16_square
    m1Residues (by simpa [modulus1At] using m1Fit) column
  have matrixBound' :
      matrixProductValue m1Residues Generated.r1_rns_mat column <
        16 * (modulus2At column) ^ 2 := by
    simpa [modulus2At] using matrixBound
  unfold combinedM2WideValue
  nlinarith

theorem m2_combined_residues_contract
    (left right m1Residues : Fin 8 → Nat)
    (leftFit : ∀ index, left index < 2 * modulus2At index)
    (rightFit : ∀ index, right index < 2 * modulus2At index)
    (m1Fit : ∀ index, m1Residues index < 2 * modulus1At index) :
    ∀ column,
      (m2CombinedResidues left right m1Residues column :
          ZMod (modulus2At column)) * ifmaRadix =
        combinedM2WideValue left right m1Residues column ∧
      m2CombinedResidues left right m1Residues column <
        2 * modulus2At column := by
  intro column
  let wide := combinedM2WideValue left right m1Residues column
  have wideBound : wide ≤ 20 * (modulus2At column) ^ 2 :=
    (combined_m2_wide_lt_20_square left right m1Residues
      leftFit rightFit m1Fit column).le
  have split := Nat.mod_add_div wide ifmaRadix
  have splitReordered :
      (wide / ifmaRadix) * ifmaRadix + wide % ifmaRadix = wide := by
    simpa [Nat.add_comm, Nat.mul_comm] using split
  have reducer := selected_combined_reducer_contract
    (wide / ifmaRadix) (wide % ifmaRadix)
    (modulus2At column) (reducer2FactorAt column)
    (modulus2_member column) (reducer2_factor_relation column)
    (by
      rw [splitReordered]
      simpa [Nat.pow_two, Nat.mul_assoc] using wideBound)
  have splitCast := congrArg
    (fun value : Nat => (value : ZMod (modulus2At column))) splitReordered
  push_cast at splitCast
  unfold m2CombinedResidues
  dsimp only
  constructor
  · rw [reducer.1]
    simpa [wide] using splitCast
  · exact reducer.2

theorem r1_product_from_residues_cross_basis_output
    (left right : Nat) (residues : Fin 8 → Nat)
    (reconstruction :
      r1ProductFromResidues left right residues * Generated.M1 =
        left * right + Generated.target_modulus *
          r1RawLiftFromResidues residues) :
    ((left * right + Generated.target_modulus *
        r1RawLiftFromResidues residues) * ifmaRadix *
          crossBasisInverse * crossBasisInverse) % Generated.M2 =
      (r1ProductFromResidues left right residues * ifmaRadix *
        crossBasisInverse) % Generated.M2 := by
  exact cross_basis_output_generic
    Generated.M1 Generated.M2 crossBasisInverse ifmaRadix
    (left * right + Generated.target_modulus *
      r1RawLiftFromResidues residues)
    (r1ProductFromResidues left right residues)
    (by norm_num [Generated.M2]) reconstruction cross_basis_inverse_relation

set_option maxHeartbeats 600000 in
theorem m2_combined_residues_represent
    (leftValue rightValue : Nat)
    (left right m1Residues : Fin 8 → Nat)
    (leftRepresent : ∀ index,
      Nat.ModEq (modulus2At index) (left index)
        (leftValue * ifmaRadix * crossBasisInverse))
    (rightRepresent : ∀ index,
      Nat.ModEq (modulus2At index) (right index)
        (rightValue * ifmaRadix * crossBasisInverse))
    (leftFit : ∀ index, left index < 2 * modulus2At index)
    (rightFit : ∀ index, right index < 2 * modulus2At index)
    (m1Fit : ∀ index, m1Residues index < 2 * modulus1At index)
    (m1Represent : ∀ index,
      Nat.ModEq (modulus1At index) (m1Residues index)
        (r1ReducedProduct leftValue rightValue)) :
    ∀ column, Nat.ModEq (modulus2At column)
      (m2CombinedResidues left right m1Residues column)
      (r2InputValue
        (r1ProductFromResidues leftValue rightValue m1Residues)) := by
  intro column
  let modulus := modulus2At column
  let representative := r1ProductFromResidues
    leftValue rightValue m1Residues
  have productFacts := r1_product_from_standard_residues
    leftValue rightValue m1Residues m1Represent
  have reducedContract := (m2_combined_residues_contract
    left right m1Residues leftFit rightFit m1Fit column).1
  have leftCast := (ZMod.natCast_eq_natCast_iff _ _ modulus).2
    (leftRepresent column)
  have rightCast := (ZMod.natCast_eq_natCast_iff _ _ modulus).2
    (rightRepresent column)
  have matrixCast :
      (matrixProductValue m1Residues Generated.r1_rns_mat column :
          ZMod modulus) = generatedR1LaneValue m1Residues column := by
    unfold matrixProductValue generatedR1LaneValue modulus modulus2At
    push_cast
    rfl
  have rawMatrix := generated_r1_lane_represents_raw_crt
    m1Residues column
  have combinedCast :
      (combinedM2WideValue left right m1Residues column : ZMod modulus) =
        (leftValue * rightValue + Generated.target_modulus *
          r1RawLiftFromResidues m1Residues) * ifmaRadix * ifmaRadix *
            crossBasisInverse * crossBasisInverse := by
    unfold combinedM2WideValue
    push_cast
    rw [leftCast, rightCast, matrixCast, rawMatrix]
    unfold r1RawLiftFromResidues
    rw [r1_postmult_relation]
    push_cast
    simp only [montgomeryRadix, ifmaRadix]
    ring_nf
    ac_rfl
  have crossGlobal := r1_product_from_residues_cross_basis_output
    leftValue rightValue m1Residues productFacts.1
  have crossSmall :
      ((leftValue * rightValue + Generated.target_modulus *
          r1RawLiftFromResidues m1Residues) * ifmaRadix *
            crossBasisInverse * crossBasisInverse : ZMod modulus) =
        representative * ifmaRadix * crossBasisInverse := by
    have reducedEquality :
        (((leftValue * rightValue + Generated.target_modulus *
            r1RawLiftFromResidues m1Residues) * ifmaRadix *
              crossBasisInverse * crossBasisInverse) % Generated.M2) % modulus =
          ((representative * ifmaRadix * crossBasisInverse) % Generated.M2) %
            modulus := congrArg (fun value => value % modulus) crossGlobal
    rw [Nat.mod_mod_of_dvd _ (modulus2_dvd_product column),
      Nat.mod_mod_of_dvd _ (modulus2_dvd_product column)] at reducedEquality
    have castEquality :=
      (ZMod.natCast_eq_natCast_iff _ _ modulus).2 reducedEquality
    simpa only [Nat.cast_mul, Nat.cast_add] using castEquality
  have targetCast :
      (r2InputValue representative : ZMod modulus) =
        representative * ifmaRadix * crossBasisInverse := by
    unfold r2InputValue
    have reduced :
        (representative * ifmaRadix * crossBasisInverse % Generated.M2) %
            modulus =
          (representative * ifmaRadix * crossBasisInverse) % modulus :=
      Nat.mod_mod_of_dvd _ (modulus2_dvd_product column)
    have castEquality :=
      (ZMod.natCast_eq_natCast_iff _ _ modulus).2 reduced
    simp only [montgomeryRadix, ifmaRadix] at castEquality ⊢
    simpa only [Nat.cast_mul] using castEquality
  have outputTimesRadix :
      (m2CombinedResidues left right m1Residues column : ZMod modulus) *
          ifmaRadix =
        (r2InputValue representative : ZMod modulus) * ifmaRadix := by
    rw [reducedContract, combinedCast, targetCast]
    have crossSmallTimes := congrArg
      (fun value : ZMod modulus => value * ifmaRadix) crossSmall
    calc
      _ = (((leftValue : ZMod modulus) * rightValue +
            Generated.target_modulus * r1RawLiftFromResidues m1Residues) *
          ifmaRadix * crossBasisInverse * crossBasisInverse) * ifmaRadix := by
        ring
      _ = ((representative : ZMod modulus) * ifmaRadix *
          crossBasisInverse) * ifmaRadix := crossSmallTimes
  have outputCast :=
    (ifma_radix_is_unit_modulus2 column).mul_right_cancel outputTimesRadix
  exact (ZMod.natCast_eq_natCast_iff _ _ modulus).1 outputCast

def m1ExpandedResidues
    (m2Residues : Fin 8 → Nat) (column : Fin 8) : Nat :=
  let wide := generatedR2ShiftedWideValue m2Residues column
  selectedReduceMasked (wide / ifmaRadix) (wide % ifmaRadix)
    (modulus1At column) (reducer1FactorAt column)

theorem m1_expanded_residues_contract
    (m2Residues : Fin 8 → Nat)
    (m2Fit : ∀ index, m2Residues index < 2 * modulus2At index) :
    ∀ column,
      (m1ExpandedResidues m2Residues column :
          ZMod (modulus1At column)) * ifmaRadix =
        generatedR2ShiftedWideValue m2Residues column ∧
      m1ExpandedResidues m2Residues column < 2 * modulus1At column := by
  intro column
  let wide := generatedR2ShiftedWideValue m2Residues column
  have wideBound : wide ≤ 32 * (modulus1At column) ^ 2 := by
    have result := generated_r2_shifted_wide_lt_32_square
      m2Residues (by simpa [modulus2At] using m2Fit) column
    simpa [wide, modulus1At] using result.le
  have split := Nat.mod_add_div wide ifmaRadix
  have splitReordered :
      (wide / ifmaRadix) * ifmaRadix + wide % ifmaRadix = wide := by
    simpa [Nat.add_comm, Nat.mul_comm] using split
  have reducer := selected_expansion_reducer_contract
    (wide / ifmaRadix) (wide % ifmaRadix)
    (modulus1At column) (reducer1FactorAt column)
    (modulus1_member column) (reducer1_factor_relation column)
    (by
      rw [splitReordered]
      simpa [Nat.pow_two, Nat.mul_assoc] using wideBound)
  have splitCast := congrArg
    (fun value : Nat => (value : ZMod (modulus1At column))) splitReordered
  push_cast at splitCast
  unfold m1ExpandedResidues
  dsimp only
  constructor
  · rw [reducer.1]
    simpa [wide] using splitCast
  · exact reducer.2

theorem m1_expanded_residues_represent
    (representative : Nat) (m2Residues : Fin 8 → Nat)
    (representativeBound : representative < 18 * Generated.target_modulus)
    (m2Fit : ∀ index, m2Residues index < 2 * modulus2At index)
    (m2Represent : ∀ index,
      Nat.ModEq (modulus2At index) (m2Residues index)
        (r2InputValue representative)) :
    ∀ column, Nat.ModEq (modulus1At column)
      (m1ExpandedResidues m2Residues column)
      (representative * ifmaRadix * Generated.qr_premult) := by
  intro column
  have reducerContract :=
    (m1_expanded_residues_contract m2Residues m2Fit column).1
  rw [modulus1At] at reducerContract
  have wideCast := generated_r2_shifted_wide_represents_lane
    m2Residues column
  have laneCast := generated_r2_standard_lane_represents_input
    representative m2Residues representativeBound
    (by simpa [modulus2At] using m2Fit)
    (by simpa [modulus2At] using m2Represent) column
  have wideCast' :
      (generatedR2ShiftedWideValue m2Residues column :
        ZMod (listAt8 Generated.moduli1 column.val)) =
        generatedR2LaneValue m2Residues column := by
    exact wideCast
  have laneCast' :
      generatedR2LaneValue m2Residues column =
        (Generated.r2_postmult :
          ZMod (listAt8 Generated.moduli1 column.val)) * representative := by
    exact laneCast
  have postmultRelation :
      (Generated.r2_postmult :
          ZMod (listAt8 Generated.moduli1 column.val)) * representative =
        representative * ifmaRadix * ifmaRadix * Generated.qr_premult := by
    rw [r2_postmult_relation]
    simp only [montgomeryRadix, ifmaRadix]
    push_cast
    ring
  have targetExpansion :
      ((representative * ifmaRadix * Generated.qr_premult : Nat) :
          ZMod (listAt8 Generated.moduli1 column.val)) =
        (representative : ZMod (listAt8 Generated.moduli1 column.val)) * ifmaRadix *
          Generated.qr_premult := by
    push_cast
    rfl
  have outputTimesRadix :
      (m1ExpandedResidues m2Residues column :
          ZMod (listAt8 Generated.moduli1 column.val)) * ifmaRadix =
        ((representative : ZMod (listAt8 Generated.moduli1 column.val)) *
          ifmaRadix * Generated.qr_premult) * ifmaRadix := by
    calc
      _ = (generatedR2ShiftedWideValue m2Residues column :
          ZMod (listAt8 Generated.moduli1 column.val)) :=
        reducerContract
      _ = generatedR2LaneValue m2Residues column := wideCast'
      _ = (Generated.r2_postmult :
          ZMod (listAt8 Generated.moduli1 column.val)) *
          representative := laneCast'
      _ = (representative :
          ZMod (listAt8 Generated.moduli1 column.val)) * ifmaRadix * ifmaRadix *
          Generated.qr_premult := postmultRelation
      _ = ((representative : ZMod (listAt8 Generated.moduli1 column.val)) *
          ifmaRadix * Generated.qr_premult) * ifmaRadix := by ring
  have radixUnit := ifma_radix_is_unit_modulus1 column
  rw [modulus1At] at radixUnit
  have outputCast :=
    radixUnit.mul_right_cancel outputTimesRadix
  have outputNatCast := outputCast.trans targetExpansion.symm
  have result := (ZMod.natCast_eq_natCast_iff _ _
    (listAt8 Generated.moduli1 column.val)).1 outputNatCast
  simpa [modulus1At] using result

theorem r2_input_represents_m2_encoding
    (representative : Nat) (index : Fin 8) :
    Nat.ModEq (modulus2At index) (r2InputValue representative)
      (representative * ifmaRadix * crossBasisInverse) := by
  unfold r2InputValue
  change
    (representative * montgomeryRadix * crossBasisInverse % Generated.M2) %
        modulus2At index =
      (representative * ifmaRadix * crossBasisInverse) % modulus2At index
  rw [Nat.mod_mod_of_dvd _ (modulus2_dvd_product index)]
  rfl

/--
The selected scalar-shaped MatrixNoK multiplication pipeline.  Inputs may use
either standard representative in each lane (`< 2m_i`).  The result returns
both bases below `2m_i`, both bases encode the same representative, and that
representative is the field product modulo p.
-/
theorem selected_base_field_pipeline
    (leftValue rightValue : Nat)
    (leftM1 rightM1 leftM2 rightM2 : Fin 8 → Nat)
    (leftValueBound : leftValue < 40 * Generated.target_modulus)
    (rightValueBound : rightValue < 40 * Generated.target_modulus)
    (leftM1Fit : ∀ index, leftM1 index < 2 * modulus1At index)
    (rightM1Fit : ∀ index, rightM1 index < 2 * modulus1At index)
    (leftM2Fit : ∀ index, leftM2 index < 2 * modulus2At index)
    (rightM2Fit : ∀ index, rightM2 index < 2 * modulus2At index)
    (leftM1Represent : ∀ index,
      Nat.ModEq (modulus1At index) (leftM1 index)
        (leftValue * ifmaRadix * Generated.qr_premult))
    (rightM1Represent : ∀ index,
      Nat.ModEq (modulus1At index) (rightM1 index)
        (rightValue * ifmaRadix * Generated.qr_premult))
    (leftM2Represent : ∀ index,
      Nat.ModEq (modulus2At index) (leftM2 index)
        (leftValue * ifmaRadix * crossBasisInverse))
    (rightM2Represent : ∀ index,
      Nat.ModEq (modulus2At index) (rightM2 index)
        (rightValue * ifmaRadix * crossBasisInverse)) :
    let intermediateM1 := m1ProductResidues leftM1 rightM1
    let representative := r1ProductFromResidues
      leftValue rightValue intermediateM1
    let resultM2 := m2CombinedResidues leftM2 rightM2 intermediateM1
    let resultM1 := m1ExpandedResidues resultM2
    representative < 18 * Generated.target_modulus ∧
      (representative * Generated.M1) % Generated.target_modulus =
        leftValue * rightValue % Generated.target_modulus ∧
      (∀ index, resultM1 index < 2 * modulus1At index ∧
        Nat.ModEq (modulus1At index) (resultM1 index)
          (representative * ifmaRadix * Generated.qr_premult)) ∧
      (∀ index, resultM2 index < 2 * modulus2At index ∧
        Nat.ModEq (modulus2At index) (resultM2 index)
          (representative * ifmaRadix * crossBasisInverse)) := by
  dsimp only
  let intermediateM1 := m1ProductResidues leftM1 rightM1
  let representative := r1ProductFromResidues
    leftValue rightValue intermediateM1
  let resultM2 := m2CombinedResidues leftM2 rightM2 intermediateM1
  let resultM1 := m1ExpandedResidues resultM2
  have intermediateFit := m1_product_residues_fit
    leftM1 rightM1 leftM1Fit rightM1Fit
  have intermediateRepresent := m1_product_residues_represent
    leftValue rightValue leftM1 rightM1 leftM1Represent rightM1Represent
    leftM1Fit rightM1Fit
  have representativeFacts := r1_product_from_standard_residues
    leftValue rightValue intermediateM1 intermediateRepresent
  have representativeBound := r1_product_from_standard_residues_lt_18p
    leftValue rightValue intermediateM1 leftValueBound rightValueBound
    intermediateFit
  have resultM2Contract := m2_combined_residues_contract
    leftM2 rightM2 intermediateM1 leftM2Fit rightM2Fit intermediateFit
  have resultM2Represent := m2_combined_residues_represent
    leftValue rightValue leftM2 rightM2 intermediateM1
    leftM2Represent rightM2Represent leftM2Fit rightM2Fit
    intermediateFit intermediateRepresent
  have resultM2Fit : ∀ index, resultM2 index < 2 * modulus2At index :=
    fun index => (resultM2Contract index).2
  have resultM1Fit := fun index =>
    (m1_expanded_residues_contract resultM2 resultM2Fit index).2
  have resultM1Represent := m1_expanded_residues_represent
    representative resultM2 representativeBound resultM2Fit resultM2Represent
  refine ⟨representativeBound, representativeFacts.2, ?_, ?_⟩
  · intro index
    exact ⟨resultM1Fit index, resultM1Represent index⟩
  · intro index
    refine ⟨resultM2Fit index, ?_⟩
    exact (resultM2Represent index).trans
      (r2_input_represents_m2_encoding representative index)

end CatBLST.AVX.MatrixNoK
