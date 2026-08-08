/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Constants
import Mathlib.Tactic

/-!
# MatrixNoK fixed-table regeneration

This file specifies the CRT table generator independently of the checked-in
tables.  The two inverse lists are proof witnesses, not trusted results: Lean
checks their defining congruences and the complete CRT basis residue pattern
before using them to regenerate every matrix, correction, and quotient entry.
-/

namespace CatBLST.AVX.MatrixNoK

set_option exponentiation.threshold 512
open scoped BigOperators

def crtInverses1 : List Nat := [
  936022637907613,
  501741240388833,
  802292453280408,
  783531998470659,
  611457125683783,
  384708461081466,
  325114214131677,
  158731496425870
]

def crtInverses2 : List Nat := [
  362325748323670,
  1062357115556509,
  1016115151632927,
  527391238345870,
  807765137545781,
  713177739680012,
  1092935200549882,
  47432202578403
]

def crtFactors (moduli inverses : List Nat) : List Nat :=
  List.zipWith (fun modulus inverse =>
    let rest := moduli.prod / modulus
    rest * inverse % moduli.prod) moduli inverses

def inverseWitnessesValid (moduli inverses : List Nat) : Prop :=
  moduli.length = inverses.length ∧
  List.Forall₂ (fun modulus inverse =>
    let rest := moduli.prod / modulus
    rest % modulus * inverse % modulus = 1) moduli inverses

theorem crt_inverse_witnesses1 :
    inverseWitnessesValid Generated.moduli1 crtInverses1 := by
  norm_num [inverseWitnessesValid, crtInverses1, Generated.moduli1]

theorem crt_inverse_witnesses2 :
    inverseWitnessesValid Generated.moduli2 crtInverses2 := by
  norm_num [inverseWitnessesValid, crtInverses2, Generated.moduli2]

def indices8 : List Nat := [0, 1, 2, 3, 4, 5, 6, 7]

def crtResidueMatrix8 (moduli factors : List Nat) : List (List Nat) :=
  indices8.map fun row =>
    indices8.map fun column =>
      factors.getD row 0 % moduli.getD column 1

def identityResidueMatrix8 : List (List Nat) :=
  indices8.map fun row =>
    indices8.map fun column =>
      if row = column then 1 else 0

def crtResiduePattern (moduli factors : List Nat) : Prop :=
  moduli.length = 8 ∧ factors.length = 8 ∧
  crtResidueMatrix8 moduli factors = identityResidueMatrix8

theorem crt_factors1_residue_pattern :
    crtResiduePattern Generated.moduli1
      (crtFactors Generated.moduli1 crtInverses1) := by
  norm_num [crtResiduePattern, crtResidueMatrix8, identityResidueMatrix8,
    indices8, crtFactors, crtInverses1, Generated.moduli1]

theorem crt_factors2_residue_pattern :
    crtResiduePattern Generated.moduli2
      (crtFactors Generated.moduli2 crtInverses2) := by
  norm_num [crtResiduePattern, crtResidueMatrix8, identityResidueMatrix8,
    indices8, crtFactors, crtInverses2, Generated.moduli2]

def radixDigits (value bits count : Nat) : List Nat :=
  (List.range count).map fun digit =>
    value / 2 ^ (bits * digit) % 2 ^ bits

def encodeOutput
    (value : Nat) (moduli : List Nat) (bits digits : Nat) : List Nat :=
  moduli.flatMap fun modulus => radixDigits (value % modulus) bits digits

def shiftedFactors
    (moduli inverses : List Nat)
    (premult inputDigits inputBits : Nat) : List Nat :=
  (List.range moduli.length).flatMap fun inputIndex =>
    (List.range inputDigits).map fun digit =>
      premult * 2 ^ (inputBits * digit) *
        (crtFactors moduli inverses).getD inputIndex 0 % moduli.prod

def regenerateMatrix
    (moduliIn inverses moduliOut : List Nat)
    (premult postmult inputDigits inputBits outputDigits outputBits : Nat) :
    List (List Nat) :=
  (shiftedFactors moduliIn inverses premult inputDigits inputBits).map fun shifted =>
    encodeOutput (shifted * postmult % moduliOut.prod)
      moduliOut outputBits outputDigits

def negMod (value modulus : Nat) : Nat :=
  (modulus - value % modulus) % modulus

def regenerateCorrection
    (moduliIn moduliOut : List Nat)
    (postmult outputDigits outputBits : Nat) : List Nat :=
  encodeOutput
    (negMod moduliIn.prod moduliOut.prod * postmult % moduliOut.prod)
    moduliOut outputBits outputDigits

def ceilDiv (numerator denominator : Nat) : Nat :=
  numerator / denominator + if numerator % denominator = 0 then 0 else 1

/-- A ceiling quotient encloses the exact numerator within one denominator. -/
theorem ceilDiv_error_bounds (numerator denominator : Nat)
    (denominatorPositive : 0 < denominator) :
    numerator ≤ ceilDiv numerator denominator * denominator ∧
      ceilDiv numerator denominator * denominator < numerator + denominator := by
  have division := Nat.mod_add_div numerator denominator
  have remainderBound := Nat.mod_lt numerator denominatorPositive
  unfold ceilDiv
  by_cases exactDivision : numerator % denominator = 0
  · simp [exactDivision] at division ⊢
    have quotientProduct : numerator / denominator * denominator = numerator := by
      rw [Nat.mul_comm]
      exact division
    rw [quotientProduct]
    omega
  · simp [exactDivision]
    have remainderPositive : 0 < numerator % denominator := Nat.pos_of_ne_zero exactDivision
    constructor <;> nlinarith

def regenerateQuotientEstimates
    (moduli inverses : List Nat)
    (premult inputDigits inputBits : Nat) : List Nat :=
  (shiftedFactors moduli inverses premult inputDigits inputBits).map fun shifted =>
    ceilDiv (shifted * 2 ^ 64) moduli.prod

def regenerateFloorQuotientEstimates
    (moduli inverses : List Nat)
    (premult inputDigits inputBits : Nat) : List Nat :=
  (shiftedFactors moduli inverses premult inputDigits inputBits).map fun shifted =>
    (shifted * 2 ^ 64) / moduli.prod

theorem r1_matrix_regenerates :
    regenerateMatrix Generated.moduli1 crtInverses1 Generated.moduli2
      Generated.r1_premult Generated.r1_postmult 1 52 1 52 =
        Generated.r1_rns_mat := by
  decide

theorem r2_matrix_regenerates :
    regenerateMatrix Generated.moduli2 crtInverses2 Generated.moduli1
      Generated.r2_premult Generated.r2_postmult 1 52 1 52 =
        Generated.r2_rns_mat := by
  decide

theorem convert_to_matrix_regenerates :
    regenerateMatrix [Generated.target_modulus] [1] Generated.moduli2
      Generated.convert_to_premult Generated.convert_to_postmult 8 50 1 52 =
        Generated.convert_to_rns_mat := by
  decide

theorem convert_from_matrix_regenerates :
    regenerateMatrix Generated.moduli2 crtInverses2 [Generated.target_modulus]
      Generated.convert_from_premult Generated.convert_from_postmult 1 52 8 52 =
        Generated.convert_from_rns_mat := by
  decide

theorem correction_vectors_regenerate :
    regenerateCorrection Generated.moduli1 Generated.moduli2
        Generated.r1_postmult 1 52 = Generated.r1_correction ∧
    regenerateCorrection Generated.moduli2 Generated.moduli1
        Generated.r2_postmult 1 52 = Generated.r2_correction ∧
    regenerateCorrection [Generated.target_modulus] Generated.moduli2
        Generated.convert_to_postmult 1 52 = Generated.convert_to_correction ∧
    regenerateCorrection Generated.moduli2 [Generated.target_modulus]
        Generated.convert_from_postmult 8 52 = Generated.convert_from_correction := by
  decide

theorem quotient_estimates_regenerate :
    regenerateQuotientEstimates Generated.moduli1 crtInverses1
        Generated.r1_premult 1 52 = Generated.r1_shifted_quotient_estimations ∧
    regenerateQuotientEstimates Generated.moduli2 crtInverses2
        Generated.r2_premult 1 52 = Generated.r2_shifted_quotient_estimations ∧
    regenerateFloorQuotientEstimates [Generated.target_modulus] [1]
        Generated.convert_to_premult 8 50 =
          Generated.convert_to_shifted_quotient_estimations ∧
    regenerateQuotientEstimates Generated.moduli2 crtInverses2
        Generated.convert_from_premult 1 52 =
          Generated.convert_from_shifted_quotient_estimations := by
  decide

def allShiftedFactorsInRange
    (moduli inverses : List Nat)
    (premult inputDigits inputBits : Nat) : Prop :=
  ∀ value ∈ shiftedFactors moduli inverses premult inputDigits inputBits,
    value < moduli.prod

theorem shiftedFactors_in_range
    (moduli inverses : List Nat)
    (premult inputDigits inputBits : Nat)
    (productPositive : 0 < moduli.prod) :
    allShiftedFactorsInRange moduli inverses premult inputDigits inputBits := by
  intro value membership
  simp only [shiftedFactors, List.mem_flatMap, List.mem_range,
    List.mem_map] at membership
  obtain ⟨inputIndex, _, digit, _, rfl⟩ := membership
  exact Nat.mod_lt _ productPositive

theorem shifted_factors_are_reduced :
    allShiftedFactorsInRange Generated.moduli1 crtInverses1
        Generated.r1_premult 1 52 ∧
    allShiftedFactorsInRange Generated.moduli2 crtInverses2
        Generated.r2_premult 1 52 ∧
    allShiftedFactorsInRange [Generated.target_modulus] [1]
        Generated.convert_to_premult 8 50 ∧
    allShiftedFactorsInRange Generated.moduli2 crtInverses2
        Generated.convert_from_premult 1 52 := by
  constructor
  · apply shiftedFactors_in_range
    norm_num [Generated.moduli1]
  constructor
  · apply shiftedFactors_in_range
    norm_num [Generated.moduli2]
  constructor
  · apply shiftedFactors_in_range
    norm_num [Generated.target_modulus]
  · apply shiftedFactors_in_range
    norm_num [Generated.moduli2]

/-- Eight 52-bit inputs times eight 64-bit estimates cannot wrap `uint128_t`. -/
theorem quotient_numerator_fits_u128 :
    8 * (2 ^ 52 - 1) * (2 ^ 64 - 1) < 2 ^ 128 := by
  norm_num

theorem quotient_numerator_fits_u128_for_inputs
    (scalars estimates : Fin 8 → Nat)
    (scalarsFit : ∀ index, scalars index < 2 ^ 52)
    (estimatesFit : ∀ index, estimates index < 2 ^ 64) :
    (∑ index, scalars index * estimates index) < 2 ^ 128 := by
  calc
    (∑ index, scalars index * estimates index) ≤
        ∑ _index : Fin 8, (2 ^ 52 - 1) * (2 ^ 64 - 1) := by
      apply Finset.sum_le_sum
      intro index _
      have scalarFits := scalarsFit index
      have estimateFits := estimatesFit index
      apply Nat.mul_le_mul
      · omega
      · omega
    _ = 8 * (2 ^ 52 - 1) * (2 ^ 64 - 1) := by norm_num
    _ < 2 ^ 128 := quotient_numerator_fits_u128

end CatBLST.AVX.MatrixNoK
