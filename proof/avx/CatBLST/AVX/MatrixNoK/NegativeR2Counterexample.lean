/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.SignedPipeline
import Mathlib.Tactic

/-!
# Boundary of the signed R2 contract

The fixed-point R2 expansion is not correct for every abstract negative lift
in `(-p, 0)`. This exact counterexample records that narrow contract boundary.
Canonical ConvertTo no longer reaches it: its repaired downward quotient
coefficients select a nonnegative representative in `[0, 2p)`.
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

end CatBLST.AVX.MatrixNoK
