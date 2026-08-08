/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Fp2
import CatBLST.AVX.MatrixNoK.SignedPipeline
import Mathlib.Tactic

/-!
# Reachable signed boundary in selected Fp2 reduction

The selected Fp2 real coordinate is not guaranteed to choose a nonnegative
correction-free R1 representative.  Splitting the eight M1 moduli into two
four-factor field elements gives a valid multiplication whose real expression
is exactly `-M1`; every reduced M1 residue is zero and the chosen representative
is `-1`.  This is not a C++ failure: the existing R2 fixed-point schedule is
proved below to preserve this reachable negative representative.
-/

namespace CatBLST.AVX.MatrixNoK

def fp2M1LeftFactor : Nat :=
  listAt8 Generated.moduli1 0 * listAt8 Generated.moduli1 1 *
    listAt8 Generated.moduli1 2 * listAt8 Generated.moduli1 3

def fp2M1RightFactor : Nat :=
  listAt8 Generated.moduli1 4 * listAt8 Generated.moduli1 5 *
    listAt8 Generated.moduli1 6 * listAt8 Generated.moduli1 7

theorem fp2_m1_split_factors_are_field_values :
    fp2M1LeftFactor < Generated.target_modulus ∧
      fp2M1RightFactor < Generated.target_modulus := by
  norm_num [fp2M1LeftFactor, fp2M1RightFactor, listAt8,
    Generated.moduli1, Generated.target_modulus]

theorem fp2_m1_split_factors_multiply_to_m1 :
    fp2M1LeftFactor * fp2M1RightFactor = Generated.M1 := by
  norm_num [fp2M1LeftFactor, fp2M1RightFactor, listAt8,
    Generated.moduli1, Generated.M1]

def fp2NegativeUnitM1Residues (_index : Fin 8) : Nat := 0

theorem fp2_negative_unit_reduced_value_is_zero :
    signedR1ReducedValue (-(Generated.M1 : Int)) = 0 := by
  have castValue := int_mod_nat_cast
    (-(Generated.M1 : Int) * ifmaRadix * Generated.qr_premult *
      Generated.qr_premult)
    Generated.M1 (by norm_num [Generated.M1])
  change (signedR1ReducedValue (-(Generated.M1 : Int)) : Int) = _
    at castValue
  have remainderZero :
      (-(Generated.M1 : Int) * ifmaRadix * Generated.qr_premult *
        Generated.qr_premult) % (Generated.M1 : Int) = 0 := by
    rw [show -(Generated.M1 : Int) * ifmaRadix *
        Generated.qr_premult * Generated.qr_premult =
      (Generated.M1 : Int) *
        (-ifmaRadix * Generated.qr_premult * Generated.qr_premult) by ring]
    simp
  have castZero :
      (signedR1ReducedValue (-(Generated.M1 : Int)) : Int) = 0 :=
    castValue.trans remainderZero
  exact_mod_cast castZero

theorem fp2_negative_unit_residues_represent (index : Fin 8) :
    Nat.ModEq (modulus1At index) (fp2NegativeUnitM1Residues index)
      (signedR1ReducedValue (-(Generated.M1 : Int))) := by
  rw [fp2_negative_unit_reduced_value_is_zero]
  exact Nat.ModEq.refl 0

theorem fp2_real_expression_reaches_negative_m1 :
    (0 : Int) * 0 - fp2M1LeftFactor * fp2M1RightFactor =
      -(Generated.M1 : Int) := by
  push_cast
  rw [zero_sub]
  have productCast :
      (fp2M1LeftFactor : Int) * fp2M1RightFactor = Generated.M1 := by
    exact_mod_cast fp2_m1_split_factors_multiply_to_m1
  rw [productCast]

theorem fp2_correction_free_r1_can_select_negative_one :
    r1SignedValueFromResidues (-(Generated.M1 : Int))
      fp2NegativeUnitM1Residues = -1 := by
  change r1SignedValueFromResidues (-(Generated.M1 : Int))
    (fun _index => 0) = -1
  have reconstruction := r1_signed_value_reconstruction
    (-(Generated.M1 : Int)) fp2NegativeUnitM1Residues
    fp2_negative_unit_residues_represent
  unfold r1SignedValueNumerator fp2NegativeUnitM1Residues
    r1RawLiftFromResidues exactCrtNumerator at reconstruction
  simp at reconstruction
  apply mul_right_cancel₀ (show (Generated.M1 : Int) ≠ 0 by
    norm_num [Generated.M1])
  simpa using reconstruction

def fp2NegativeOneR2Residues (index : Fin 8) : Nat :=
  listAt8 [
    671982632611421,
    6880907853788,
    821190241501970,
    990310359749765,
    12021489756635,
    903164310377547,
    445358240441908,
    251549293611906
  ] index.val

theorem fp2_negative_one_r2_input_value :
    signedR2InputValue (-1) =
      1505393867229161485538073664013819807292133703609172883291565475159772376069635928860656827332187946571555257358576663718 := by
  unfold signedR2InputValue
  have castValue := int_mod_nat_cast
    (-(1 : Int) * ifmaRadix * crossBasisInverse)
    Generated.M2 (by norm_num [Generated.M2])
  norm_num [ifmaRadix, montgomeryRadix, crossBasisInverse,
    Generated.convert_to_postmult, Generated.M2] at castValue
  exact_mod_cast castValue

theorem fp2_negative_one_r2_residues_are_canonical (index : Fin 8) :
    fp2NegativeOneR2Residues index = canonicalSignedR2Residues (-1) index := by
  rw [canonicalSignedR2Residues, fp2_negative_one_r2_input_value]
  fin_cases index <;>
    norm_num [fp2NegativeOneR2Residues, Generated.moduli2,
      modulus2At, listAt8]

theorem fp2_negative_one_r2_is_preserved :
    correctedCrtRepresentative fp2NegativeOneR2Residues
        r2ShiftedFactor Generated.M2 = -1 := by
  norm_num [fp2NegativeOneR2Residues, correctedCrtRepresentative,
    exactCrtNumerator, estimatedCrtQuotient, approximateCrtNumerator,
    r2ShiftedFactor, ceilDiv, uint64Modulus,
    Generated.M2, Generated.moduli2, Generated.r2_premult, crtFactors2,
    crtFactors, crtInverses2, Generated.r2_shifted_quotient_estimations,
    listAt8, Fin.sum_univ_succ]

end CatBLST.AVX.MatrixNoK
