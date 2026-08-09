/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Fp12FrobeniusConstants
import Mathlib.Tactic

/-!
# Physical bounds for the selected Fp12 Frobenius schedules

The selected maps form at most two products in any reduced coordinate.  An
input standard lane is below `2m`; every generated Frobenius-constant lane is
below `m`.  Consequently each product is below `2m²` and each paired sum is
below `4m²`, far inside the literal C++ `ready<MAX_ADD=800>` boundary.

The coordinates negated without multiplication use the literal
`standard_negate` path: `2m - value`, addition of the generated `40p` standard
zero, and `reduce_auto<2>`'s one conditional modulus subtraction.
-/

namespace CatBLST.AVX.MatrixNoK

def fp12FrobeniusMaximumProducts : Nat := 2
def fp12FrobeniusPhysicalBound : Nat := 4
def fp12FrobeniusRepresentativeBound : Nat := 80
def fp12FrobeniusStandardZeroRepresentative : Nat :=
  40 * Generated.target_modulus

theorem fp12_frobenius_bound_ledger :
    fp12FrobeniusPhysicalBound =
        fp12FrobeniusMaximumProducts * standardBound ∧
      fp12FrobeniusPhysicalBound ≤ fp12PhysicalCompletedBound ∧
      fp12FrobeniusPhysicalBound ≤ fp12ReadyBound ∧
      fp12FrobeniusRepresentativeBound < fp12MaximumCompletedRnsBound := by
  norm_num [fp12FrobeniusPhysicalBound, fp12FrobeniusMaximumProducts,
    fp12FrobeniusRepresentativeBound, standardBound,
    fp12PhysicalCompletedBound, fp12ReadyBound,
    fp12MaximumCompletedRnsBound]

theorem fp12_frobenius_single_product_bound
    (input constant modulus : Nat)
    (inputFit : input < 2 * modulus)
    (constantFit : constant < modulus) :
    input * constant < 2 * modulus * modulus := by
  nlinarith

theorem fp12_frobenius_two_product_bound
    (input1 constant1 input2 constant2 modulus : Nat)
    (input1Fit : input1 < 2 * modulus)
    (constant1Fit : constant1 < modulus)
    (input2Fit : input2 < 2 * modulus)
    (constant2Fit : constant2 < modulus) :
    input1 * constant1 + input2 * constant2 <
      fp12FrobeniusPhysicalBound * modulus * modulus := by
  have first := fp12_frobenius_single_product_bound
    input1 constant1 modulus input1Fit constant1Fit
  have second := fp12_frobenius_single_product_bound
    input2 constant2 modulus input2Fit constant2Fit
  norm_num [fp12FrobeniusPhysicalBound]
  nlinarith

theorem fp12_frobenius_single_product_encoding
    (input constant inputValue constantValue scale modulus : Nat)
    (inputRepresents :
      (input : ZMod modulus) = inputValue * scale)
    (constantRepresents :
      (constant : ZMod modulus) = constantValue * scale) :
    ((input * constant : Nat) : ZMod modulus) =
      (inputValue * constantValue : ZMod modulus) * scale * scale := by
  push_cast
  rw [inputRepresents, constantRepresents]
  ring

theorem fp12_frobenius_two_product_encoding
    (input1 constant1 input2 constant2
      inputValue1 constantValue1 inputValue2 constantValue2 scale modulus : Nat)
    (input1Represents :
      (input1 : ZMod modulus) = inputValue1 * scale)
    (constant1Represents :
      (constant1 : ZMod modulus) = constantValue1 * scale)
    (input2Represents :
      (input2 : ZMod modulus) = inputValue2 * scale)
    (constant2Represents :
      (constant2 : ZMod modulus) = constantValue2 * scale) :
    ((input1 * constant1 + input2 * constant2 : Nat) : ZMod modulus) =
      (inputValue1 * constantValue1 + inputValue2 * constantValue2 :
        ZMod modulus) * scale * scale := by
  push_cast
  rw [input1Represents, constant1Represents,
    input2Represents, constant2Represents]
  ring

def generatedStandardNegateOffsetM1 (lane : Fin 8) : Nat :=
  listAt8 Generated.modulus_in_rns_m1 lane.val

def generatedStandardNegateOffsetM2 (lane : Fin 8) : Nat :=
  listAt8 Generated.modulus_in_rns_m2 lane.val

theorem generated_standard_negate_offset_contract :
    ∀ lane : Fin 8,
      generatedStandardNegateOffsetM1 lane < modulus1At lane ∧
        (generatedStandardNegateOffsetM1 lane : ZMod (modulus1At lane)) =
          (fp12FrobeniusStandardZeroRepresentative * ifmaRadix *
            Generated.qr_premult : Nat) ∧
        generatedStandardNegateOffsetM2 lane < modulus2At lane ∧
        (generatedStandardNegateOffsetM2 lane : ZMod (modulus2At lane)) =
          (fp12FrobeniusStandardZeroRepresentative * ifmaRadix *
            crossBasisInverse : Nat) := by
  native_decide

/-- Literal scalar model of `min_sub_reduce<3, 2>`. -/
def cppMinSubToTwo (value modulus : Nat) : Nat :=
  if value < modulus then value else value - modulus

theorem cpp_min_sub_to_two_contract
    (value modulus : Nat)
    (_modulusPositive : 0 < modulus)
    (valueFit : value < 3 * modulus) :
    cppMinSubToTwo value modulus < 2 * modulus ∧
      (cppMinSubToTwo value modulus : ZMod modulus) = value := by
  by_cases below : value < modulus
  · simp [cppMinSubToTwo, below]
    omega
  · have modulusLe : modulus ≤ value := Nat.le_of_not_gt below
    have differenceFit : value - modulus < 2 * modulus := by omega
    constructor
    · simpa [cppMinSubToTwo, below] using differenceFit
    · simp only [cppMinSubToTwo, below, if_false]
      rw [Nat.cast_sub modulusLe]
      simp

/-- Literal scalar model of the standard-RNS-bound `standard_negate` lane. -/
def cppStandardNegateLane (value offset modulus : Nat) : Nat :=
  cppMinSubToTwo (2 * modulus - value + offset) modulus

theorem cpp_standard_negate_lane_contract
    (value offset modulus representative scale zeroRepresentative : Nat)
    (modulusPositive : 0 < modulus)
    (valueFit : value < 2 * modulus)
    (offsetFit : offset < modulus)
    (representativeLe : representative ≤ zeroRepresentative)
    (valueRepresents :
      (value : ZMod modulus) = representative * scale)
    (offsetRepresents :
      (offset : ZMod modulus) = zeroRepresentative * scale) :
    cppStandardNegateLane value offset modulus < 2 * modulus ∧
      (cppStandardNegateLane value offset modulus : ZMod modulus) =
        ((zeroRepresentative - representative : Nat) : ZMod modulus) * scale := by
  have valueLe : value ≤ 2 * modulus := by omega
  have preFit : 2 * modulus - value + offset < 3 * modulus := by omega
  have reduced := cpp_min_sub_to_two_contract
    (2 * modulus - value + offset) modulus modulusPositive preFit
  constructor
  · exact reduced.1
  · change
      (cppMinSubToTwo (2 * modulus - value + offset) modulus : ZMod modulus) = _
    rw [reduced.2]
    push_cast
    rw [Nat.cast_sub valueLe, valueRepresents, offsetRepresents,
      Nat.cast_sub representativeLe]
    simp
    ring

theorem fp12_frobenius_selected_reduction_is_exact
    (completed : Int)
    (rawM1 m1High m1Low rawM2 : Fin 8 → Nat)
    (completedNonnegative : 0 ≤ completed)
    (completedBound : completed <
      fp12FrobeniusRepresentativeBound * Generated.target_modulus *
        Generated.target_modulus)
    (rawM1Represents : ∀ index,
      (rawM1 index : ZMod (modulus1At index)) =
        completed * ifmaRadix * ifmaRadix *
          Generated.qr_premult * Generated.qr_premult)
    (rawM1Fit : ∀ index, rawM1 index ≤
      fp12FrobeniusPhysicalBound * modulus1At index * modulus1At index)
    (m1WideSplit : ∀ index,
      m1High index * ifmaRadix + m1Low index = rawM1 index)
    (rawM2Represents : ∀ index,
      (rawM2 index : ZMod (modulus2At index)) =
        completed * ifmaRadix * ifmaRadix *
          crossBasisInverse * crossBasisInverse)
    (rawM2Fit : ∀ index, rawM2 index <
      fp12FrobeniusPhysicalBound * modulus2At index * modulus2At index) :
    let m1Residues := cppFp12M1Residue m1High m1Low
    correctedCrtRepresentative
        (cppFp12M2Residue rawM2 m1Residues)
        r2ShiftedFactor Generated.M2 =
      r1SignedValueFromResidues completed m1Residues := by
  apply fp12_selected_redundant_reduction_is_exact
  · exact completedNonnegative
  · have ledger := fp12_frobenius_bound_ledger.2.2.2
    have modulusPositive : 0 < Generated.target_modulus := by
      norm_num [Generated.target_modulus]
    nlinarith
  · exact rawM1Represents
  · intro index
    have fit := rawM1Fit index
    have ledger := fp12_frobenius_bound_ledger.2.2.1
    have scaled := Nat.mul_le_mul_right
      (modulus1At index * modulus1At index) ledger
    have scaled' :
        fp12FrobeniusPhysicalBound * modulus1At index * modulus1At index ≤
          fp12ReadyBound * modulus1At index * modulus1At index := by
      simpa [Nat.mul_assoc] using scaled
    exact fit.trans scaled'
  · exact m1WideSplit
  · exact rawM2Represents
  · intro index
    have fit := rawM2Fit index
    have ledger := fp12_frobenius_bound_ledger.2.1
    have scaled := Nat.mul_le_mul_right
      (modulus2At index * modulus2At index) ledger
    have scaled' :
        fp12FrobeniusPhysicalBound * modulus2At index * modulus2At index ≤
          fp12PhysicalCompletedBound * modulus2At index * modulus2At index := by
      simpa [Nat.mul_assoc] using scaled
    exact fit.trans_le scaled'

end CatBLST.AVX.MatrixNoK
