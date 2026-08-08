/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Conversion
import Mathlib.Data.ZMod.Basic
import Mathlib.Tactic

/-!
# RNS conversion semantics

This file states the algebra performed by one output lane of a fixed RNS
matrix.  Table layout and literal equality are proved in `Regeneration`; the
theorems here show why the regenerated entry and correction formulas represent
the corrected CRT integer modulo each output modulus.
-/

namespace CatBLST.AVX.MatrixNoK

open scoped BigOperators

theorem negMod_cast (value modulus : Nat) (modulusPositive : 0 < modulus) :
    (negMod value modulus : ZMod modulus) = -(value : ZMod modulus) := by
  unfold negMod
  rw [ZMod.natCast_mod]
  have remainderAtMost : value % modulus ≤ modulus :=
    (Nat.mod_lt value modulusPositive).le
  rw [Nat.cast_sub remainderAtMost]
  simp [ZMod.natCast_mod]

def fixedMatrixEntry
    (shiftedFactor postmult outputModulus : Nat) : Nat :=
  shiftedFactor * postmult % outputModulus

def fixedCorrection
    (inputModulus postmult outputModulus : Nat) : Nat :=
  negMod inputModulus outputModulus * postmult % outputModulus

def fixedMatrixLaneValue
    (scalars shiftedFactors : Fin 8 → Nat)
    (inputModulus postmult outputModulus : Nat) : ZMod outputModulus :=
  (∑ index, (scalars index : ZMod outputModulus) *
      fixedMatrixEntry (shiftedFactors index) postmult outputModulus) +
    (estimatedCrtQuotient scalars shiftedFactors inputModulus :
      ZMod outputModulus) *
      fixedCorrection inputModulus postmult outputModulus

def floorFixedMatrixLaneValue
    (scalars shiftedFactors : Fin 8 → Nat)
    (inputModulus postmult outputModulus : Nat) : ZMod outputModulus :=
  (∑ index, (scalars index : ZMod outputModulus) *
      fixedMatrixEntry (shiftedFactors index) postmult outputModulus) +
    (floorEstimatedCrtQuotient scalars shiftedFactors inputModulus :
      ZMod outputModulus) *
      fixedCorrection inputModulus postmult outputModulus

def rawFixedMatrixLaneValue
    (scalars shiftedFactors : Fin 8 → Nat)
    (postmult outputModulus : Nat) : ZMod outputModulus :=
  ∑ index, (scalars index : ZMod outputModulus) *
    fixedMatrixEntry (shiftedFactors index) postmult outputModulus

theorem raw_fixed_matrix_lane_represents_crt
    (scalars shiftedFactors : Fin 8 → Nat)
    (postmult outputModulus : Nat) :
    rawFixedMatrixLaneValue scalars shiftedFactors postmult outputModulus =
      (postmult : ZMod outputModulus) *
        exactCrtNumerator scalars shiftedFactors := by
  unfold rawFixedMatrixLaneValue fixedMatrixEntry exactCrtNumerator
  simp only [ZMod.natCast_mod]
  push_cast
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro index _
  ring

theorem fixed_matrix_lane_represents_corrected_crt
    (scalars shiftedFactors : Fin 8 → Nat)
    (inputModulus postmult outputModulus : Nat)
    (outputModulusPositive : 0 < outputModulus) :
    fixedMatrixLaneValue scalars shiftedFactors
        inputModulus postmult outputModulus =
      (postmult : ZMod outputModulus) *
        (correctedCrtRepresentative
          scalars shiftedFactors inputModulus : Int) := by
  unfold fixedMatrixLaneValue fixedMatrixEntry fixedCorrection
  simp only [ZMod.natCast_mod]
  push_cast
  rw [negMod_cast inputModulus outputModulus outputModulusPositive]
  unfold correctedCrtRepresentative exactCrtNumerator
  push_cast
  have matrixSum :
      (∑ index, (scalars index : ZMod outputModulus) *
        ((shiftedFactors index : ZMod outputModulus) * postmult)) =
      (∑ index, (scalars index : ZMod outputModulus) * shiftedFactors index) *
        postmult := by
    rw [Finset.sum_mul]
    apply Finset.sum_congr rfl
    intro index _
    ring
  rw [matrixSum]
  ring

theorem floor_fixed_matrix_lane_represents_corrected_crt
    (scalars shiftedFactors : Fin 8 → Nat)
    (inputModulus postmult outputModulus : Nat)
    (outputModulusPositive : 0 < outputModulus) :
    floorFixedMatrixLaneValue scalars shiftedFactors
        inputModulus postmult outputModulus =
      (postmult : ZMod outputModulus) *
        floorCorrectedCrtRepresentative
          scalars shiftedFactors inputModulus := by
  unfold floorFixedMatrixLaneValue fixedMatrixEntry fixedCorrection
  simp only [ZMod.natCast_mod]
  push_cast
  rw [negMod_cast inputModulus outputModulus outputModulusPositive]
  unfold floorCorrectedCrtRepresentative exactCrtNumerator
  have correctionBelow :
      floorEstimatedCrtQuotient scalars shiftedFactors inputModulus *
          inputModulus ≤
        ∑ index, scalars index * shiftedFactors index := by
    by_cases inputPositive : 0 < inputModulus
    · have quotientBound := floor_estimated_quotient_le_exact
        scalars shiftedFactors inputModulus inputPositive
      have productBound := Nat.mul_le_mul_right inputModulus quotientBound
      have exactBound := Nat.div_mul_le_self
        (∑ index, scalars index * shiftedFactors index) inputModulus
      exact productBound.trans (by
        simpa [exactCrtQuotient, exactCrtNumerator, Nat.mul_comm] using exactBound)
    · simp at inputPositive
      simp [inputPositive, floorEstimatedCrtQuotient,
        approximateCrtNumeratorFloor]
  rw [Nat.cast_sub correctionBelow]
  push_cast
  have matrixSum :
      (∑ index, (scalars index : ZMod outputModulus) *
        ((shiftedFactors index : ZMod outputModulus) * postmult)) =
      (∑ index, (scalars index : ZMod outputModulus) * shiftedFactors index) *
        postmult := by
    rw [Finset.sum_mul]
    apply Finset.sum_congr rfl
    intro index _
    ring
  rw [matrixSum]
  ring

end CatBLST.AVX.MatrixNoK
