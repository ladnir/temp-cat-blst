/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.SignedPipeline
import CatBLST.AVX.MatrixNoK.CanonicalImport
import Mathlib.Tactic

/-!
# Boundary of the signed R2 contract

The fixed-point R2 expansion is not correct for every abstract negative lift
in `(-p, 0)`.  This exact counterexample prevents the reachable-pipeline proof
from silently assuming that wider statement.  It does not establish a C++
defect: the remaining obligation is to characterize the representatives that
the preceding R1 schedule can actually emit.
-/

namespace CatBLST.AVX.MatrixNoK

set_option exponentiation.threshold 512

def negativeR2CounterexampleMagnitude : Nat :=
  2652105521236084289755357371233918909542389279325034051769198635113044629510026440333065428518657785524845196455279

def negativeR2CounterexampleResidues (index : Fin 8) : Nat :=
  listAt8 [
    85131525,
    3976650,
    29473586786712,
    25140345,
    54691584,
    18979761,
    18067814,
    12280267
  ] index.val

theorem negative_r2_counterexample_is_in_signed_field_range :
    0 < negativeR2CounterexampleMagnitude ∧
      negativeR2CounterexampleMagnitude < Generated.target_modulus := by
  norm_num [negativeR2CounterexampleMagnitude, Generated.target_modulus]

theorem negative_r2_counterexample_input_value :
    signedR2InputValue (-(negativeR2CounterexampleMagnitude : Int)) =
      1167277358685670452076721592973614357942727014784180263116639308050535784177362747036649306092735403932395542566642358238 := by
  unfold signedR2InputValue
  have castValue := int_mod_nat_cast
    (-(negativeR2CounterexampleMagnitude : Int) * ifmaRadix *
      crossBasisInverse) Generated.M2 (by norm_num [Generated.M2])
  norm_num [negativeR2CounterexampleMagnitude, ifmaRadix, montgomeryRadix,
    crossBasisInverse, Generated.convert_to_postmult, Generated.M2] at castValue
  exact_mod_cast castValue

theorem negative_r2_counterexample_residues_are_canonical
    (index : Fin 8) :
    negativeR2CounterexampleResidues index =
      canonicalSignedR2Residues
        (-(negativeR2CounterexampleMagnitude : Int)) index := by
  rw [canonicalSignedR2Residues,
    negative_r2_counterexample_input_value]
  fin_cases index <;>
    norm_num [negativeR2CounterexampleResidues, Generated.moduli2,
      modulus2At, listAt8]

theorem negative_r2_counterexample_selects_lower_quotient :
    estimatedCrtQuotient negativeR2CounterexampleResidues r2ShiftedFactor
        Generated.M2 =
      exactCrtQuotient negativeR2CounterexampleResidues r2ShiftedFactor
        Generated.M2 := by
  norm_num [estimatedCrtQuotient, approximateCrtNumerator,
    exactCrtQuotient, exactCrtNumerator, negativeR2CounterexampleResidues,
    r2ShiftedFactor, ceilDiv, uint64Modulus, Generated.M2,
    Generated.moduli2, Generated.r2_premult, crtFactors2, crtFactors,
    crtInverses2, Generated.r2_shifted_quotient_estimations, listAt8,
    Fin.sum_univ_succ]

theorem negative_r2_counterexample_is_not_preserved :
    correctedCrtRepresentative negativeR2CounterexampleResidues
        r2ShiftedFactor Generated.M2 ≠
      -(negativeR2CounterexampleMagnitude : Int) := by
  norm_num [correctedCrtRepresentative, exactCrtNumerator,
    estimatedCrtQuotient, approximateCrtNumerator,
    negativeR2CounterexampleResidues, negativeR2CounterexampleMagnitude,
    r2ShiftedFactor, ceilDiv, uint64Modulus, Generated.M2,
    Generated.moduli2, Generated.r2_premult, crtFactors2, crtFactors,
    crtInverses2, Generated.r2_shifted_quotient_estimations, listAt8,
    Fin.sum_univ_succ]

/-! ## Reachable canonical-import counterexample -/

def reachableR2CounterexampleMagnitude : Nat :=
  235613689341516701328257298793030861480281319249132620885720573137412674019990281389671108000420691382528538236

def reachableR2CounterexampleCanonicalValue : Nat :=
  1267253245981839142264863935160608057326608172767180188451362181929717619810821363803276830181705476116338914753674

def reachableR2CounterexampleDigits (index : Fin 8) : Nat :=
  listAt8 [
    500616628141194,
    1029601075259562,
    686812189624822,
    769198589660543,
    581522751939018,
    26173626664305,
    89497429892085,
    552541535
  ] index.val

def reachableR2CounterexampleResidues (index : Fin 8) : Nat :=
  listAt8 [
    390364834,
    14226308,
    136333126987290,
    19524759,
    13127718,
    733971630,
    560530109,
    244433910
  ] index.val

theorem reachable_r2_counterexample_digits_fit (index : Fin 8) :
    reachableR2CounterexampleDigits index < 2 ^ 50 := by
  fin_cases index <;>
    norm_num [reachableR2CounterexampleDigits, listAt8]

theorem reachable_r2_counterexample_digits_encode_canonical_value :
    canonicalDigitValue reachableR2CounterexampleDigits =
      reachableR2CounterexampleCanonicalValue := by
  norm_num [canonicalDigitValue, reachableR2CounterexampleDigits, listAt8,
    reachableR2CounterexampleCanonicalValue, Fin.sum_univ_succ]

theorem reachable_r2_counterexample_canonical_value_is_valid :
    reachableR2CounterexampleCanonicalValue < Generated.target_modulus := by
  norm_num [reachableR2CounterexampleCanonicalValue,
    Generated.target_modulus]

theorem reachable_r2_counterexample_import_is_negative :
    convertToCorrectedRepresentative reachableR2CounterexampleDigits =
      -(reachableR2CounterexampleMagnitude : Int) := by
  norm_num [convertToCorrectedRepresentative, correctedCrtRepresentative,
    exactCrtNumerator, estimatedCrtQuotient, approximateCrtNumerator,
    reachableR2CounterexampleDigits, reachableR2CounterexampleMagnitude,
    convertToShiftedFactor, ceilDiv, uint64Modulus,
    Generated.convert_to_premult, Generated.target_modulus, listAt8,
    Fin.sum_univ_succ]

theorem reachable_r2_counterexample_import_represents_canonical_value :
    (-(reachableR2CounterexampleMagnitude : Int) :
        ZMod Generated.target_modulus) =
      Generated.convert_to_premult *
        reachableR2CounterexampleCanonicalValue := by
  calc
    _ = (convertToCorrectedRepresentative reachableR2CounterexampleDigits :
          ZMod Generated.target_modulus) := by
      have castEquality := congrArg
        (fun value : Int => (value : ZMod Generated.target_modulus))
        reachable_r2_counterexample_import_is_negative
      simpa using castEquality.symm
    _ = Generated.convert_to_premult *
          canonicalDigitValue reachableR2CounterexampleDigits :=
      convert_to_corrected_represents_digits reachableR2CounterexampleDigits
    _ = _ := by
      rw [reachable_r2_counterexample_digits_encode_canonical_value]

theorem reachable_r2_counterexample_input_value :
    signedR2InputValue (-(reachableR2CounterexampleMagnitude : Int)) =
      616327079703076581242832078304941328016735964488459999848989771852762966190805087676499037199305175125193321929304050194 := by
  unfold signedR2InputValue
  have castValue := int_mod_nat_cast
    (-(reachableR2CounterexampleMagnitude : Int) * ifmaRadix *
      crossBasisInverse) Generated.M2 (by norm_num [Generated.M2])
  norm_num [reachableR2CounterexampleMagnitude, ifmaRadix, montgomeryRadix,
    crossBasisInverse, Generated.convert_to_postmult, Generated.M2] at castValue
  exact_mod_cast castValue

theorem reachable_r2_counterexample_residues_are_canonical
    (index : Fin 8) :
    reachableR2CounterexampleResidues index =
      canonicalSignedR2Residues
        (-(reachableR2CounterexampleMagnitude : Int)) index := by
  rw [canonicalSignedR2Residues, reachable_r2_counterexample_input_value]
  fin_cases index <;>
    norm_num [reachableR2CounterexampleResidues, Generated.moduli2,
      modulus2At, listAt8]

/-- This is the exact Montgomery boundary used by
`from_canonical_50bit_digits`: reducing each ConvertTo accumulator by M2 emits
the stated canonical residue. -/
theorem reachable_r2_counterexample_convert_to_reduces_to_residue
    (index : Fin 8) :
    (reachableR2CounterexampleResidues index :
        ZMod (modulus2At index)) * ifmaRadix =
      generatedConvertToLaneValue reachableR2CounterexampleDigits index := by
  rw [reachable_r2_counterexample_residues_are_canonical]
  unfold canonicalSignedR2Residues
  rw [ZMod.natCast_mod]
  rw [signed_r2_input_lane]
  rw [generated_convert_to_lane_represents_signed_import]
  rw [reachable_r2_counterexample_import_is_negative]
  have postmult := conversion_multiplier_relations.2.1
  rw [postmult]
  simp only [ifmaRadix, montgomeryRadix]
  push_cast
  simp [mul_comm, mul_left_comm]
  rfl

theorem reachable_r2_counterexample_selects_lower_quotient :
    estimatedCrtQuotient reachableR2CounterexampleResidues r2ShiftedFactor
        Generated.M2 =
      exactCrtQuotient reachableR2CounterexampleResidues r2ShiftedFactor
        Generated.M2 := by
  norm_num [estimatedCrtQuotient, approximateCrtNumerator,
    exactCrtQuotient, exactCrtNumerator, reachableR2CounterexampleResidues,
    r2ShiftedFactor, ceilDiv, uint64Modulus, Generated.M2,
    Generated.moduli2, Generated.r2_premult, crtFactors2, crtFactors,
    crtInverses2, Generated.r2_shifted_quotient_estimations, listAt8,
    Fin.sum_univ_succ]

theorem reachable_r2_counterexample_corrects_to_positive_complement :
    correctedCrtRepresentative reachableR2CounterexampleResidues
        r2ShiftedFactor Generated.M2 =
      (Generated.M2 : Int) - reachableR2CounterexampleMagnitude := by
  norm_num [correctedCrtRepresentative, exactCrtNumerator,
    estimatedCrtQuotient, approximateCrtNumerator,
    reachableR2CounterexampleResidues, reachableR2CounterexampleMagnitude,
    r2ShiftedFactor, ceilDiv, uint64Modulus, Generated.M2,
    Generated.moduli2, Generated.r2_premult, crtFactors2, crtFactors,
    crtInverses2, Generated.r2_shifted_quotient_estimations, listAt8,
    Fin.sum_univ_succ]

theorem reachable_r2_counterexample_m1_lane_zero_is_wrong :
    generatedR2LaneValue reachableR2CounterexampleResidues (0 : Fin 8) ≠
      (Generated.r2_postmult :
          ZMod (listAt8 Generated.moduli1 (0 : Fin 8).val)) *
        (-(reachableR2CounterexampleMagnitude : Int)) := by
  rw [generated_r2_lane_eq_fixed]
  rw [fixed_matrix_lane_represents_corrected_crt _ _ _ _ _
    (by norm_num [Generated.moduli1, listAt8])]
  rw [reachable_r2_counterexample_corrects_to_positive_complement]
  intro equality
  push_cast at equality
  have impossible :
      (Generated.r2_postmult :
          ZMod (listAt8 Generated.moduli1 (0 : Fin 8).val)) *
        Generated.M2 = 0 := by
    linear_combination equality
  have impossibleNat :
      ((Generated.r2_postmult * Generated.M2 : Nat) :
          ZMod (listAt8 Generated.moduli1 (0 : Fin 8).val)) = 0 := by
    simpa only [Nat.cast_mul] using impossible
  have divides := (ZMod.natCast_eq_zero_iff
    (Generated.r2_postmult * Generated.M2)
    (listAt8 Generated.moduli1 (0 : Fin 8).val)).1 impossibleNat
  norm_num [Generated.M2, Generated.r2_postmult, Generated.moduli1,
    listAt8] at divides

end CatBLST.AVX.MatrixNoK
