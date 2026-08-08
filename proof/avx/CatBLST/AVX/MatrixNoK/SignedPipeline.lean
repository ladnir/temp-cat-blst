/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Pipeline
import Mathlib.Tactic

/-!
# Signed MatrixNoK representatives

The physical RNS lanes are always unsigned residues.  Their common integer
lift need not be nonnegative: the fixed-point CRT correction used by canonical
import may select `u - p` instead of `u`.  This file models that signed lift
without changing the C++ representation or arithmetic.
-/

namespace CatBLST.AVX.MatrixNoK

/-- The least nonnegative natural representative of an integer modulo `m`. -/
def intModNat (value : Int) (modulus : Nat) : Nat :=
  (value % (modulus : Int)).toNat

theorem int_mod_nat_cast (value : Int) (modulus : Nat)
    (modulusPositive : 0 < modulus) :
    (intModNat value modulus : Int) = value % (modulus : Int) := by
  unfold intModNat
  rw [Int.toNat_of_nonneg]
  exact Int.emod_nonneg _ (by exact_mod_cast modulusPositive.ne')

theorem int_mod_nat_lt (value : Int) (modulus : Nat)
    (modulusPositive : 0 < modulus) :
    intModNat value modulus < modulus := by
  have modulusNonzero : (modulus : Int) ≠ 0 := by
    exact_mod_cast modulusPositive.ne'
  have bound := Int.emod_lt value modulusNonzero
  have cast := int_mod_nat_cast value modulus modulusPositive
  have castBound : (intModNat value modulus : Int) < (modulus : Int) := by
    rw [cast]
    simpa using bound
  exact_mod_cast castBound

theorem int_mod_nat_neg_of_pos_lt
    (value modulus : Nat) (valuePositive : 0 < value)
    (valueBound : value < modulus) :
    intModNat (-(value : Int)) modulus = modulus - value := by
  have modulusPositive : 0 < modulus := valuePositive.trans valueBound
  have castResult := int_mod_nat_cast (-(value : Int)) modulus modulusPositive
  have differenceBoundNat : modulus - value < modulus := by omega
  have differenceNonnegative :
      (0 : Int) ≤ ((modulus - value : Nat) : Int) := by positivity
  have differenceBound :
      ((modulus - value : Nat) : Int) < (modulus : Int) := by
    exact_mod_cast differenceBoundNat
  have negativeRewrite :
      -(value : Int) =
        ((modulus - value : Nat) : Int) - (modulus : Int) := by
    rw [Nat.cast_sub valueBound.le]
    ring
  have remainder :
      (-(value : Int)) % (modulus : Int) =
        ((modulus - value : Nat) : Int) := by
    rw [negativeRewrite, Int.sub_emod]
    simp only [Int.emod_self, sub_zero, Int.emod_emod]
    exact Int.emod_eq_of_lt differenceNonnegative differenceBound
  have castEquality :
      (intModNat (-(value : Int)) modulus : Int) =
        ((modulus - value : Nat) : Int) := castResult.trans remainder
  exact_mod_cast castEquality

theorem int_mod_nat_zmod (value : Int) (modulus : Nat)
    (modulusPositive : 0 < modulus) :
    (intModNat value modulus : ZMod modulus) = value := by
  calc
    (intModNat value modulus : ZMod modulus) =
        ((intModNat value modulus : Int) : ZMod modulus) := by norm_cast
    _ = (value : ZMod modulus) := by
      rw [ZMod.intCast_eq_intCast_iff']
      rw [int_mod_nat_cast value modulus modulusPositive]
      rw [Int.emod_emod]

theorem int_mod_nat_zmod_of_dvd
    (value : Int) (smallModulus largeModulus : Nat)
    (largePositive : 0 < largeModulus)
    (divides : smallModulus ∣ largeModulus) :
    (intModNat value largeModulus : ZMod smallModulus) = value := by
  have largeCongruence :
      (intModNat value largeModulus : Int) ≡ value
        [ZMOD (largeModulus : Int)] := by
    unfold Int.ModEq
    rw [int_mod_nat_cast value largeModulus largePositive]
    rw [Int.emod_emod]
  have castDivides : (smallModulus : Int) ∣ (largeModulus : Int) := by
    exact_mod_cast divides
  have smallCongruence :
      (intModNat value largeModulus : Int) ≡ value
        [ZMOD (smallModulus : Int)] :=
    Int.modEq_of_dvd (castDivides.trans largeCongruence.dvd)
  calc
    (intModNat value largeModulus : ZMod smallModulus) =
        ((intModNat value largeModulus : Int) : ZMod smallModulus) := by
      norm_cast
    _ = (value : ZMod smallModulus) := by
      exact (ZMod.intCast_eq_intCast_iff _ _ _).2 smallCongruence

/-- M1 residue represented by multiplying two signed QR lifts. -/
def signedR1ReducedProduct (left right : Int) : Nat :=
  intModNat
    (left * right * ifmaRadix * Generated.qr_premult *
      Generated.qr_premult)
    Generated.M1

theorem signed_r1_reduced_product_lane (left right : Int) (index : Fin 8) :
    (signedR1ReducedProduct left right : ZMod (modulus1At index)) =
      left * right * ifmaRadix * Generated.qr_premult *
        Generated.qr_premult := by
  unfold signedR1ReducedProduct
  have result := int_mod_nat_zmod_of_dvd
    (left * right * (ifmaRadix : Int) * Generated.qr_premult *
      Generated.qr_premult)
    (modulus1At index) Generated.M1
    (by norm_num [Generated.M1]) (modulus1_dvd_product index)
  push_cast at result
  exact result

theorem m1_product_residues_represent_signed
    (leftValue rightValue : Int) (left right : Fin 8 → Nat)
    (leftRepresent : ∀ index,
      (left index : ZMod (modulus1At index)) =
        leftValue * ifmaRadix * Generated.qr_premult)
    (rightRepresent : ∀ index,
      (right index : ZMod (modulus1At index)) =
        rightValue * ifmaRadix * Generated.qr_premult)
    (leftFit : ∀ index, left index < 2 * modulus1At index)
    (rightFit : ∀ index, right index < 2 * modulus1At index) :
    ∀ index, Nat.ModEq (modulus1At index)
      (m1ProductResidues left right index)
      (signedR1ReducedProduct leftValue rightValue) := by
  intro index
  let modulus := modulus1At index
  have reducedContract := (product_reduced_lane_contract
    (left index) (right index) modulus (reducer1FactorAt index)
    (modulus1_member index) (reducer1_factor_relation index)
    (leftFit index) (rightFit index)).1
  have targetCast := signed_r1_reduced_product_lane
    leftValue rightValue index
  have targetTimesRadix :
      (signedR1ReducedProduct leftValue rightValue : ZMod modulus) *
          ifmaRadix =
        leftValue * rightValue * ifmaRadix * ifmaRadix *
          Generated.qr_premult * Generated.qr_premult := by
    rw [targetCast]
    ring
  have equalTimesRadix :
      (m1ProductResidues left right index : ZMod modulus) * ifmaRadix =
        (signedR1ReducedProduct leftValue rightValue : ZMod modulus) *
          ifmaRadix := by
    unfold m1ProductResidues
    rw [reducedContract, leftRepresent index, rightRepresent index,
      targetTimesRadix]
    ring
  have equalCast :=
    (ifma_radix_is_unit_modulus1 index).mul_right_cancel equalTimesRadix
  exact (ZMod.natCast_eq_natCast_iff _ _ modulus).1 equalCast

def r1SignedNumerator
    (left right : Int) (residues : Fin 8 → Nat) : Int :=
  left * right + Generated.target_modulus * r1RawLiftFromResidues residues

def r1SignedProductFromResidues
    (left right : Int) (residues : Fin 8 → Nat) : Int :=
  r1SignedNumerator left right residues / Generated.M1

set_option maxRecDepth 10000 in
theorem r1_signed_numerator_divisible
    (left right : Int) (residues : Fin 8 → Nat)
    (residuesRepresent : ∀ index,
      Nat.ModEq (modulus1At index) (residues index)
        (signedR1ReducedProduct left right)) :
    (Generated.M1 : Int) ∣ r1SignedNumerator left right residues := by
  have rawRelation := r1_raw_lift_from_residues_relation
    residues (signedR1ReducedProduct left right)
      (by simpa [modulus1At] using residuesRepresent)
  have rawCast :
      (r1RawLiftFromResidues residues : ZMod Generated.M1) =
        Generated.r1_premult * signedR1ReducedProduct left right := by
    have castRelation := (ZMod.natCast_eq_natCast_iff
      (r1RawLiftFromResidues residues)
      (Generated.r1_premult * signedR1ReducedProduct left right)
      Generated.M1).2 rawRelation
    simpa only [Nat.cast_mul] using castRelation
  have reducedCast := int_mod_nat_zmod
    (left * right * (ifmaRadix : Int) * Generated.qr_premult *
      Generated.qr_premult)
    Generated.M1 (by norm_num [Generated.M1])
  change (signedR1ReducedProduct left right : ZMod Generated.M1) = _
    at reducedCast
  push_cast at reducedCast
  have inverseCast :
      ((Generated.qr_premult * Generated.qr_premult *
          Generated.r1_premult *
          (Generated.M1 - Generated.target_modulus) * ifmaRadix : Nat) :
          ZMod Generated.M1) = 1 := by
    apply (ZMod.natCast_eq_natCast_iff _ _ Generated.M1).2
    change
      (Generated.qr_premult * Generated.qr_premult *
          Generated.r1_premult *
          (Generated.M1 - Generated.target_modulus) * ifmaRadix) %
        Generated.M1 = 1 % Generated.M1
    have oneMod : 1 % Generated.M1 = 1 :=
      Nat.mod_eq_of_lt (by norm_num [Generated.M1])
    rw [oneMod]
    change
      (Generated.qr_premult * Generated.qr_premult *
          Generated.r1_premult *
          (Generated.M1 - Generated.target_modulus) * montgomeryRadix) %
        Generated.M1 = 1
    exact r1_premult_relation
  have inverseCastExpanded := inverseCast
  push_cast at inverseCastExpanded
  have targetBelowM1 : Generated.target_modulus ≤ Generated.M1 := by
    norm_num [Generated.target_modulus, Generated.M1]
  have subtractionCast :
      ((Generated.M1 - Generated.target_modulus : Nat) :
          ZMod Generated.M1) = -Generated.target_modulus := by
    rw [Nat.cast_sub targetBelowM1]
    simp
  have cancellation :
      (Generated.target_modulus : ZMod Generated.M1) *
          Generated.r1_premult * ifmaRadix * Generated.qr_premult *
            Generated.qr_premult = -1 := by
    calc
      _ = -((Generated.qr_premult : ZMod Generated.M1) *
          Generated.qr_premult * Generated.r1_premult *
          ((Generated.M1 - Generated.target_modulus : Nat) :
            ZMod Generated.M1) * ifmaRadix) := by
        rw [subtractionCast]
        ring
      _ = -1 := by rw [inverseCastExpanded]
  have numeratorCast :
      (r1SignedNumerator left right residues : ZMod Generated.M1) = 0 := by
    unfold r1SignedNumerator
    push_cast
    rw [rawCast, reducedCast]
    calc
      _ = (left : ZMod Generated.M1) * right *
          (1 + Generated.target_modulus * Generated.r1_premult *
            ifmaRadix * Generated.qr_premult * Generated.qr_premult) := by
        ring
      _ = 0 := by rw [cancellation]; ring
  exact (ZMod.intCast_zmod_eq_zero_iff_dvd _ _).1 numeratorCast

theorem r1_signed_product_reconstruction
    (left right : Int) (residues : Fin 8 → Nat)
    (residuesRepresent : ∀ index,
      Nat.ModEq (modulus1At index) (residues index)
        (signedR1ReducedProduct left right)) :
    r1SignedProductFromResidues left right residues * Generated.M1 =
      r1SignedNumerator left right residues := by
  unfold r1SignedProductFromResidues
  exact Int.ediv_mul_cancel
    (r1_signed_numerator_divisible left right residues residuesRepresent)

theorem r1_signed_product_is_field_multiplication
    (left right : Int) (residues : Fin 8 → Nat)
    (residuesRepresent : ∀ index,
      Nat.ModEq (modulus1At index) (residues index)
        (signedR1ReducedProduct left right)) :
    (r1SignedProductFromResidues left right residues :
        ZMod Generated.target_modulus) * Generated.M1 = left * right := by
  have reconstruction := r1_signed_product_reconstruction
    left right residues residuesRepresent
  have castEquality := congrArg
    (fun value : Int => (value : ZMod Generated.target_modulus)) reconstruction
  unfold r1SignedNumerator at castEquality
  push_cast at castEquality
  simpa using castEquality

theorem signed_product_abs_bound
    (left right bound : Int)
    (boundPositive : 0 < bound)
    (leftBound : |left| < bound)
    (rightBound : |right| < bound) :
    |left * right| < bound * bound := by
  rw [abs_mul]
  have first : |left| * |right| ≤ |left| * bound :=
    mul_le_mul_of_nonneg_left rightBound.le (abs_nonneg left)
  have second : |left| * bound < bound * bound :=
    mul_lt_mul_of_pos_right leftBound boundPositive
  exact first.trans_lt second

theorem r1_signed_product_from_standard_residues_range
    (left right : Int) (residues : Fin 8 → Nat)
    (leftBound : |left| < 40 * Generated.target_modulus)
    (rightBound : |right| < 40 * Generated.target_modulus)
    (residuesFit : ∀ index,
      residues index < 2 * listAt8 Generated.moduli1 index.val)
    (residuesRepresent : ∀ index,
      Nat.ModEq (modulus1At index) (residues index)
        (signedR1ReducedProduct left right)) :
    -(Generated.target_modulus : Int) <
        r1SignedProductFromResidues left right residues ∧
      r1SignedProductFromResidues left right residues <
        18 * Generated.target_modulus := by
  have targetPositive : (0 : Int) < Generated.target_modulus := by
    norm_num [Generated.target_modulus]
  have modulusPositive : (0 : Int) < Generated.M1 := by
    norm_num [Generated.M1]
  have dynamicRange :
      (1600 : Int) * Generated.target_modulus < Generated.M1 := by
    norm_num [Generated.target_modulus, Generated.M1]
  have productAbs := signed_product_abs_bound left right
    (40 * (Generated.target_modulus : Int)) (by positivity)
    (by simpa using leftBound) (by simpa using rightBound)
  have productBounds := (abs_lt.mp productAbs)
  have liftBoundNat := r1_standard_raw_lift_lt_17_product residues
    (by simpa [modulus1At] using residuesFit)
  have liftBound :
      (r1RawLiftFromResidues residues : Int) < 17 * Generated.M1 := by
    exact_mod_cast liftBoundNat
  have liftNonnegative :
      (0 : Int) ≤ r1RawLiftFromResidues residues := by
    exact_mod_cast Nat.zero_le (r1RawLiftFromResidues residues)
  have divisible := r1_signed_numerator_divisible
    left right residues residuesRepresent
  have numeratorLower :
      (Generated.M1 : Int) * (-Generated.target_modulus) <
        r1SignedNumerator left right residues := by
    unfold r1SignedNumerator
    nlinarith
  have numeratorUpper :
      r1SignedNumerator left right residues <
        (Generated.M1 : Int) * (18 * Generated.target_modulus) := by
    unfold r1SignedNumerator
    nlinarith
  constructor
  · unfold r1SignedProductFromResidues
    exact (Int.lt_div_iff_of_dvd_of_pos modulusPositive divisible).2
      numeratorLower
  · unfold r1SignedProductFromResidues
    exact (Int.div_lt_iff_of_dvd_of_pos modulusPositive divisible).2
      numeratorUpper

/-- The canonical M2 residue of a signed QR representative. -/
def signedR2InputValue (representative : Int) : Nat :=
  intModNat
    (representative * ifmaRadix * crossBasisInverse)
    Generated.M2

theorem signed_r2_input_lane (representative : Int) (index : Fin 8) :
    (signedR2InputValue representative : ZMod (modulus2At index)) =
      representative * ifmaRadix * crossBasisInverse := by
  unfold signedR2InputValue
  have result := int_mod_nat_zmod_of_dvd
    (representative * (ifmaRadix : Int) * crossBasisInverse)
    (modulus2At index) Generated.M2
    (by norm_num [Generated.M2]) (modulus2_dvd_product index)
  push_cast at result
  exact result

theorem signed_r2_input_scale_relation (representative : Int) :
    (Generated.r2_premult * signedR2InputValue representative :
        ZMod Generated.M2) = representative := by
  have inputCast := int_mod_nat_zmod
    (representative * (ifmaRadix : Int) * crossBasisInverse)
    Generated.M2 (by norm_num [Generated.M2])
  change (signedR2InputValue representative : ZMod Generated.M2) = _
    at inputCast
  push_cast at inputCast
  have premultCast :
      (Generated.r2_premult : ZMod Generated.M2) * ifmaRadix =
        Generated.M1 := by
    have castRelation :
        ((Generated.r2_premult * ifmaRadix : Nat) : ZMod Generated.M2) =
          (Generated.M1 : ZMod Generated.M2) := by
      apply (ZMod.natCast_eq_natCast_iff _ _ Generated.M2).2
      change
        (Generated.r2_premult * montgomeryRadix) % Generated.M2 =
          Generated.M1 % Generated.M2
      exact r2_premult_relation
    simpa only [Nat.cast_mul] using castRelation
  have crossCast :
      (Generated.M1 : ZMod Generated.M2) * crossBasisInverse = 1 := by
    have castRelation :
        ((Generated.M1 * crossBasisInverse : Nat) : ZMod Generated.M2) =
          (1 : ZMod Generated.M2) := by
      apply (ZMod.natCast_eq_natCast_iff _ _ Generated.M2).2
      change
        (Generated.M1 * crossBasisInverse) % Generated.M2 =
          1 % Generated.M2
      have oneMod : 1 % Generated.M2 = 1 :=
        Nat.mod_eq_of_lt (by norm_num [Generated.M2])
      rw [oneMod]
      exact cross_basis_inverse_relation
    simpa only [Nat.cast_mul] using castRelation
  rw [inputCast]
  calc
    _ = (representative : ZMod Generated.M2) *
        (Generated.r2_premult * ifmaRadix) * crossBasisInverse := by ring
    _ = (representative : ZMod Generated.M2) * Generated.M1 *
        crossBasisInverse := by rw [premultCast]
    _ = representative := by rw [mul_assoc, crossCast, mul_one]

theorem r2_raw_lift_from_residues_remainder_signed
    (representative : Int) (residues : Fin 8 → Nat)
    (residuesRepresent : ∀ index,
      Nat.ModEq (modulus2At index) (residues index)
        (signedR2InputValue representative)) :
    r2RawLiftFromResidues residues % Generated.M2 =
      intModNat representative Generated.M2 := by
  have reconstructed := exact_crt2_residues_mod_product residues
    (signedR2InputValue representative)
    (by simpa [modulus2At] using residuesRepresent)
  have scaled :=
    (Nat.ModEq.refl (n := Generated.M2) Generated.r2_premult).mul
      reconstructed
  have rawRelation :=
    (r2_raw_lift_from_residues_scaled residues).trans scaled
  have rawCast :
      (r2RawLiftFromResidues residues : ZMod Generated.M2) =
        Generated.r2_premult * signedR2InputValue representative := by
    have castRelation := (ZMod.natCast_eq_natCast_iff
      (r2RawLiftFromResidues residues)
      (Generated.r2_premult * signedR2InputValue representative)
      Generated.M2).2 rawRelation
    simpa only [Nat.cast_mul] using castRelation
  have targetCast := signed_r2_input_scale_relation representative
  have canonicalCast := int_mod_nat_zmod representative Generated.M2
    (by norm_num [Generated.M2])
  have equalityCast :
      (r2RawLiftFromResidues residues : ZMod Generated.M2) =
        intModNat representative Generated.M2 := by
    rw [rawCast, targetCast, canonicalCast]
  have congruence :=
    (ZMod.natCast_eq_natCast_iff _ _ Generated.M2).1 equalityCast
  unfold Nat.ModEq at congruence
  rw [Nat.mod_eq_of_lt (int_mod_nat_lt representative Generated.M2
    (by norm_num [Generated.M2]))] at congruence
  exact congruence

theorem signed_r2_input_of_nonnegative
    (representative : Int) (representativeNonnegative : 0 ≤ representative) :
    signedR2InputValue representative =
      r2InputValue representative.toNat := by
  have leftCast := int_mod_nat_cast
    (representative * (ifmaRadix : Int) * crossBasisInverse)
    Generated.M2 (by norm_num [Generated.M2])
  change (signedR2InputValue representative : Int) = _ at leftCast
  have rightCast :
      (r2InputValue representative.toNat : Int) =
        representative * ifmaRadix * crossBasisInverse % Generated.M2 := by
    unfold r2InputValue
    push_cast
    rw [Int.toNat_of_nonneg representativeNonnegative]
    rfl
  exact_mod_cast leftCast.trans rightCast.symm

theorem r2_standard_corrected_representative_is_exact_nonnegative
    (representative : Int) (residues : Fin 8 → Nat)
    (representativeNonnegative : 0 ≤ representative)
    (representativeBound : representative < 18 * Generated.target_modulus)
    (residuesFit : ∀ index,
      residues index < 2 * listAt8 Generated.moduli2 index.val)
    (residuesRepresent : ∀ index,
      Nat.ModEq (modulus2At index) (residues index)
        (signedR2InputValue representative)) :
    correctedCrtRepresentative residues r2ShiftedFactor Generated.M2 =
      representative := by
  have boundNat :
      representative.toNat < 18 * Generated.target_modulus := by
    rw [← Int.ofNat_lt]
    simpa [Int.toNat_of_nonneg representativeNonnegative] using
      representativeBound
  have representedNat : ∀ index,
      Nat.ModEq (listAt8 Generated.moduli2 index.val)
        (residues index) (r2InputValue representative.toNat) := by
    intro index
    simpa [modulus2At,
      signed_r2_input_of_nonnegative representative representativeNonnegative]
      using residuesRepresent index
  have exact := r2_standard_corrected_representative_is_exact
    representative.toNat residues boundNat residuesFit representedNat
  simpa [Int.toNat_of_nonneg representativeNonnegative] using exact

/-- The exact condition needed for a negative R2 lift: the aggregate upward
rounding error must cover the magnitude of the negative representative. -/
def r2NegativeMargin (magnitude : Nat) (residues : Fin 8 → Nat) : Prop :=
  magnitude * uint64Modulus ≤
    ∑ index, residues index *
      fixedPointCeilingError (r2ShiftedFactor index) Generated.M2

def canonicalSignedR2Residues
    (representative : Int) (index : Fin 8) : Nat :=
  signedR2InputValue representative % modulus2At index

theorem canonical_signed_r2_residue_le
    (representative : Int) (residues : Fin 8 → Nat)
    (residuesRepresent : ∀ index,
      Nat.ModEq (modulus2At index) (residues index)
        (signedR2InputValue representative)) :
    ∀ index, canonicalSignedR2Residues representative index ≤
      residues index := by
  intro index
  unfold canonicalSignedR2Residues
  have remainderEquality := residuesRepresent index
  unfold Nat.ModEq at remainderEquality
  rw [← remainderEquality]
  exact Nat.mod_le _ _

theorem r2_negative_margin_of_canonical
    (magnitude : Nat) (residues : Fin 8 → Nat)
    (residuesRepresent : ∀ index,
      Nat.ModEq (modulus2At index) (residues index)
        (signedR2InputValue (-(magnitude : Int))))
    (canonicalMargin : r2NegativeMargin magnitude
      (canonicalSignedR2Residues (-(magnitude : Int)))) :
    r2NegativeMargin magnitude residues := by
  unfold r2NegativeMargin at *
  exact canonicalMargin.trans (Finset.sum_le_sum fun index _ =>
    Nat.mul_le_mul_right
      (fixedPointCeilingError (r2ShiftedFactor index) Generated.M2)
      (canonical_signed_r2_residue_le
        (-(magnitude : Int)) residues residuesRepresent index))

theorem r2_standard_estimated_quotient_is_one_high_negative
    (magnitude : Nat) (residues : Fin 8 → Nat)
    (magnitudePositive : 0 < magnitude)
    (magnitudeBound : magnitude < Generated.M2)
    (residuesFit : ∀ index,
      residues index < 2 * listAt8 Generated.moduli2 index.val)
    (residuesRepresent : ∀ index,
      Nat.ModEq (modulus2At index) (residues index)
        (signedR2InputValue (-(magnitude : Int))))
    (margin : r2NegativeMargin magnitude residues) :
    estimatedCrtQuotient residues r2ShiftedFactor Generated.M2 =
      exactCrtQuotient residues r2ShiftedFactor Generated.M2 + 1 := by
  have remainder := r2_raw_lift_from_residues_remainder_signed
    (-(magnitude : Int)) residues residuesRepresent
  unfold r2RawLiftFromResidues at remainder
  rw [int_mod_nat_neg_of_pos_lt magnitude Generated.M2 magnitudePositive
    magnitudeBound] at remainder
  have scalarFits := r2_standard_residue_lt_ifma residues residuesFit
  apply estimated_quotient_one_high_of_margin
    residues r2ShiftedFactor Generated.M2
    (by norm_num [Generated.M2]) scalarFits
  unfold r2NegativeMargin at margin
  rw [remainder]
  have magnitudeLe : magnitude ≤ Generated.M2 := magnitudeBound.le
  simpa [Nat.sub_sub_self magnitudeLe] using margin

theorem r2_standard_corrected_representative_is_exact_negative
    (magnitude : Nat) (residues : Fin 8 → Nat)
    (magnitudePositive : 0 < magnitude)
    (magnitudeBound : magnitude < Generated.M2)
    (residuesFit : ∀ index,
      residues index < 2 * listAt8 Generated.moduli2 index.val)
    (residuesRepresent : ∀ index,
      Nat.ModEq (modulus2At index) (residues index)
        (signedR2InputValue (-(magnitude : Int))))
    (margin : r2NegativeMargin magnitude residues) :
    correctedCrtRepresentative residues r2ShiftedFactor Generated.M2 =
      -(magnitude : Int) := by
  have remainder := r2_raw_lift_from_residues_remainder_signed
    (-(magnitude : Int)) residues residuesRepresent
  unfold r2RawLiftFromResidues at remainder
  rw [int_mod_nat_neg_of_pos_lt magnitude Generated.M2 magnitudePositive
    magnitudeBound] at remainder
  have quotient := r2_standard_estimated_quotient_is_one_high_negative
    magnitude residues magnitudePositive magnitudeBound residuesFit
      residuesRepresent margin
  have split := Nat.mod_add_div
    (exactCrtNumerator residues r2ShiftedFactor) Generated.M2
  unfold correctedCrtRepresentative
  rw [quotient]
  unfold exactCrtQuotient
  have reconstructionNat :
      exactCrtNumerator residues r2ShiftedFactor =
        (Generated.M2 - magnitude) +
          (exactCrtNumerator residues r2ShiftedFactor / Generated.M2) *
            Generated.M2 := by
    calc
      exactCrtNumerator residues r2ShiftedFactor =
          exactCrtNumerator residues r2ShiftedFactor % Generated.M2 +
            Generated.M2 *
              (exactCrtNumerator residues r2ShiftedFactor / Generated.M2) :=
        split.symm
      _ = (Generated.M2 - magnitude) + Generated.M2 *
          (exactCrtNumerator residues r2ShiftedFactor / Generated.M2) := by
        rw [remainder]
      _ = (Generated.M2 - magnitude) +
          (exactCrtNumerator residues r2ShiftedFactor / Generated.M2) *
            Generated.M2 := by rw [Nat.mul_comm]
  have reconstructionInt :
      (exactCrtNumerator residues r2ShiftedFactor : Int) =
        ((Generated.M2 - magnitude : Nat) : Int) +
          ((exactCrtNumerator residues r2ShiftedFactor / Generated.M2 : Nat) :
            Int) * (Generated.M2 : Int) := by
    exact_mod_cast reconstructionNat
  rw [reconstructionInt, Nat.cast_sub magnitudeBound.le]
  push_cast
  ring

end CatBLST.AVX.MatrixNoK
