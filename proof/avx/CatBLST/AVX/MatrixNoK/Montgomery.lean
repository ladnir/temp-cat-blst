/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.KernelBounds
import Mathlib.Tactic

/-!
# MatrixNoK one-lane Montgomery reduction

This is the mathematical specification of `FixedMontgomeryReduce::reduce_small`.
`reduce_wide` first folds the stored low word's bits above radix `2^52` into
the high word; this file therefore specifies the remaining step using the low
radix digit.
-/

namespace CatBLST.AVX.MatrixNoK

def montgomeryQuotient (low factor : Nat) : Nat :=
  (low % ifmaRadix * factor) % ifmaRadix

def montgomeryReductionCarry (low modulus factor : Nat) : Nat :=
  montgomeryQuotient low factor * modulus / ifmaRadix

def montgomeryReduceStep
    (high low modulus factor : Nat) : Nat :=
  high + modulus - montgomeryReductionCarry low modulus factor

theorem montgomery_quotient_lt
    (low factor : Nat) :
    montgomeryQuotient low factor < ifmaRadix := by
  exact Nat.mod_lt _ (by norm_num [ifmaRadix])

theorem montgomery_quotient_cancels_low
    (low modulus factor : Nat)
    (inverseRelation : modulus * factor % ifmaRadix = 1) :
    montgomeryQuotient low factor * modulus % ifmaRadix =
      low % ifmaRadix := by
  unfold montgomeryQuotient
  calc
    ((low % ifmaRadix * factor) % ifmaRadix * modulus) % ifmaRadix =
        (low % ifmaRadix * factor * modulus) % ifmaRadix := by
      exact Nat.mod_mul_mod _ _ _
    _ = (low % ifmaRadix * (modulus * factor)) % ifmaRadix := by
      congr 1
      ring
    _ = (low % ifmaRadix * (modulus * factor % ifmaRadix)) %
          ifmaRadix := (Nat.mul_mod_mod _ _ _).symm
    _ = low % ifmaRadix := by
      rw [inverseRelation]
      simp

theorem montgomery_carry_lt_modulus
    (low modulus factor : Nat)
    (modulusPositive : 0 < modulus) :
    montgomeryReductionCarry low modulus factor < modulus := by
  unfold montgomeryReductionCarry
  apply (Nat.div_lt_iff_lt_mul (by norm_num [ifmaRadix])).2
  have quotientBound := montgomery_quotient_lt low factor
  simpa [Nat.mul_comm] using
    (Nat.mul_lt_mul_right modulusPositive).2 quotientBound

theorem montgomery_step_bound
    (high low modulus factor highBound : Nat)
    (modulusPositive : 0 < modulus)
    (highFits : high < highBound * modulus) :
    montgomeryReduceStep high low modulus factor <
      (highBound + 1) * modulus := by
  have carryBound := montgomery_carry_lt_modulus
    low modulus factor modulusPositive
  unfold montgomeryReduceStep
  calc
    high + modulus - montgomeryReductionCarry low modulus factor ≤
        high + modulus := Nat.sub_le _ _
    _ < highBound * modulus + modulus := Nat.add_lt_add_right highFits modulus
    _ = (highBound + 1) * modulus := by ring

theorem montgomery_step_reconstruction
    (high low modulus factor : Nat)
    (modulusPositive : 0 < modulus)
    (inverseRelation : modulus * factor % ifmaRadix = 1) :
    montgomeryReduceStep high low modulus factor * ifmaRadix +
        montgomeryQuotient low factor * modulus =
      high * ifmaRadix + modulus * ifmaRadix + low % ifmaRadix := by
  have cancellation := montgomery_quotient_cancels_low
    low modulus factor inverseRelation
  have productSplit := Nat.mod_add_div
    (montgomeryQuotient low factor * modulus) ifmaRadix
  have carryBound := montgomery_carry_lt_modulus
    low modulus factor modulusPositive
  have splitWithCarry :
      low % ifmaRadix + ifmaRadix *
          montgomeryReductionCarry low modulus factor =
        montgomeryQuotient low factor * modulus := by
    simpa [montgomeryReductionCarry, cancellation] using productSplit
  rw [Nat.mul_comm ifmaRadix
    (montgomeryReductionCarry low modulus factor)] at splitWithCarry
  have carryAtMost :
      montgomeryReductionCarry low modulus factor ≤ high + modulus := by
    omega
  have carryProductAtMost :
      montgomeryReductionCarry low modulus factor * ifmaRadix ≤
        (high + modulus) * ifmaRadix :=
    Nat.mul_le_mul_right ifmaRadix carryAtMost
  unfold montgomeryReduceStep
  calc
    (high + modulus - montgomeryReductionCarry low modulus factor) * ifmaRadix +
        montgomeryQuotient low factor * modulus =
      ((high + modulus) * ifmaRadix -
          montgomeryReductionCarry low modulus factor * ifmaRadix) +
        montgomeryQuotient low factor * modulus := by
      rw [Nat.sub_mul]
    _ = (((high + modulus) * ifmaRadix -
          montgomeryReductionCarry low modulus factor * ifmaRadix) +
        montgomeryReductionCarry low modulus factor * ifmaRadix) +
          low % ifmaRadix := by
      rw [← splitWithCarry]
      ac_rfl
    _ = (high + modulus) * ifmaRadix + low % ifmaRadix := by
      rw [Nat.sub_add_cancel carryProductAtMost]
    _ = high * ifmaRadix + modulus * ifmaRadix + low % ifmaRadix := by
      ring

theorem montgomery_step_congruent
    (high low modulus factor : Nat)
    (modulusPositive : 0 < modulus)
    (inverseRelation : modulus * factor % ifmaRadix = 1) :
    (montgomeryReduceStep high low modulus factor * ifmaRadix) % modulus =
      (high * ifmaRadix + low % ifmaRadix) % modulus := by
  have reconstruction := montgomery_step_reconstruction
    high low modulus factor modulusPositive inverseRelation
  have reduced := congrArg (fun value => value % modulus) reconstruction
  simpa [Nat.add_mod, Nat.mul_mod] using reduced

/-- Folding the stored low word into the high word preserves the wide value. -/
theorem fold_low_word_reconstruction (storedHigh storedLow : Nat) :
    (storedHigh + storedLow / ifmaRadix) * ifmaRadix +
        storedLow % ifmaRadix =
      storedHigh * ifmaRadix + storedLow := by
  have split := Nat.mod_add_div storedLow ifmaRadix
  nlinarith

theorem reduce_wide_step_congruent
    (storedHigh storedLow modulus factor : Nat)
    (modulusPositive : 0 < modulus)
    (inverseRelation : modulus * factor % ifmaRadix = 1) :
    (montgomeryReduceStep
        (storedHigh + storedLow / ifmaRadix)
        storedLow modulus factor * ifmaRadix) % modulus =
      (storedHigh * ifmaRadix + storedLow) % modulus := by
  rw [montgomery_step_congruent _ _ _ _ modulusPositive inverseRelation]
  rw [fold_low_word_reconstruction]

end CatBLST.AVX.MatrixNoK
