/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Montgomery
import Mathlib.Tactic

/-!
# Fixed-point CRT quotient semantics

The fixed matrices estimate the CRT quotient with 64 fractional bits.  This
file proves the aggregate estimate used by the selected eight-lane kernels is
either the exact floor quotient or one greater.  Consequently the generic
corrected CRT representative lies in `[-M, M)`; later RNS proofs choose the
appropriate signed or nonnegative encoding for each consumer.
-/

namespace CatBLST.AVX.MatrixNoK

open scoped BigOperators

def exactCrtNumerator
    (scalars shiftedFactors : Fin 8 → Nat) : Nat :=
  ∑ index, scalars index * shiftedFactors index

def approximateCrtNumerator
    (scalars shiftedFactors : Fin 8 → Nat) (modulus : Nat) : Nat :=
  ∑ index, scalars index *
    ceilDiv (shiftedFactors index * uint64Modulus) modulus

/-- The nonnegative error introduced by rounding one fixed-point CRT factor
upward to 64 fractional bits. -/
def fixedPointCeilingError (shiftedFactor modulus : Nat) : Nat :=
  ceilDiv (shiftedFactor * uint64Modulus) modulus * modulus -
    shiftedFactor * uint64Modulus

theorem fixed_point_ceiling_error_decomposition
    (shiftedFactor modulus : Nat) (modulusPositive : 0 < modulus) :
    shiftedFactor * uint64Modulus +
        fixedPointCeilingError shiftedFactor modulus =
      ceilDiv (shiftedFactor * uint64Modulus) modulus * modulus := by
  unfold fixedPointCeilingError
  exact Nat.add_sub_of_le
    (ceilDiv_error_bounds
      (shiftedFactor * uint64Modulus) modulus modulusPositive).1

/-- The aggregate fixed-point numerator is the exact CRT numerator plus the
sum of the per-lane upward-rounding errors.  Keeping this as an equality is
useful for proving when the selected quotient is specifically one high. -/
theorem aggregate_ceiling_error_identity
    (scalars shiftedFactors : Fin 8 → Nat) (modulus : Nat)
    (modulusPositive : 0 < modulus) :
    exactCrtNumerator scalars shiftedFactors * uint64Modulus +
        ∑ index, scalars index *
          fixedPointCeilingError (shiftedFactors index) modulus =
      approximateCrtNumerator scalars shiftedFactors modulus * modulus := by
  unfold exactCrtNumerator approximateCrtNumerator
  simp only [Finset.sum_mul, ← Finset.sum_add_distrib]
  apply Finset.sum_congr rfl
  intro index _
  calc
    scalars index * shiftedFactors index * uint64Modulus +
        scalars index * fixedPointCeilingError (shiftedFactors index) modulus =
      scalars index *
        (shiftedFactors index * uint64Modulus +
          fixedPointCeilingError (shiftedFactors index) modulus) := by ring
    _ = scalars index *
        (ceilDiv (shiftedFactors index * uint64Modulus) modulus * modulus) := by
      rw [fixed_point_ceiling_error_decomposition
        (shiftedFactors index) modulus modulusPositive]
    _ = scalars index *
        ceilDiv (shiftedFactors index * uint64Modulus) modulus * modulus := by
      ring

theorem aggregate_quotient_enclosure
    (scalars shiftedFactors : Fin 8 → Nat) (modulus : Nat)
    (modulusPositive : 0 < modulus) :
    exactCrtNumerator scalars shiftedFactors * uint64Modulus ≤
        approximateCrtNumerator scalars shiftedFactors modulus * modulus ∧
    approximateCrtNumerator scalars shiftedFactors modulus * modulus ≤
        exactCrtNumerator scalars shiftedFactors * uint64Modulus +
          (modulus - 1) * ∑ index, scalars index := by
  have perIndexLower : ∀ index : Fin 8,
      scalars index * (shiftedFactors index * uint64Modulus) ≤
        scalars index *
          (ceilDiv (shiftedFactors index * uint64Modulus) modulus * modulus) := by
    intro index
    exact Nat.mul_le_mul_left _
      (ceilDiv_error_bounds
        (shiftedFactors index * uint64Modulus) modulus modulusPositive).1
  have perIndexUpper : ∀ index : Fin 8,
      scalars index *
          (ceilDiv (shiftedFactors index * uint64Modulus) modulus * modulus) ≤
        scalars index * (shiftedFactors index * uint64Modulus) +
          scalars index * (modulus - 1) := by
    intro index
    have strictUpper := (ceilDiv_error_bounds
      (shiftedFactors index * uint64Modulus) modulus modulusPositive).2
    have nonStrictUpper :
        ceilDiv (shiftedFactors index * uint64Modulus) modulus * modulus ≤
          shiftedFactors index * uint64Modulus + (modulus - 1) := by
      omega
    calc
      scalars index *
          (ceilDiv (shiftedFactors index * uint64Modulus) modulus * modulus) ≤
        scalars index *
          (shiftedFactors index * uint64Modulus + (modulus - 1)) :=
            Nat.mul_le_mul_left _ nonStrictUpper
      _ = scalars index * (shiftedFactors index * uint64Modulus) +
          scalars index * (modulus - 1) := by ring
  constructor
  · simp only [exactCrtNumerator, approximateCrtNumerator, Finset.sum_mul]
    apply Finset.sum_le_sum
    intro index _
    simpa [Nat.mul_assoc] using perIndexLower index
  · simp only [exactCrtNumerator, approximateCrtNumerator,
      Finset.sum_mul, Finset.mul_sum]
    rw [← Finset.sum_add_distrib]
    apply Finset.sum_le_sum
    intro index _
    simpa [Nat.mul_assoc, Nat.mul_comm, Nat.mul_left_comm] using
      perIndexUpper index

theorem scalar_sum_below_fixed_point
    (scalars : Fin 8 → Nat)
    (scalarsFit : ∀ index, scalars index < ifmaRadix) :
    (∑ index, scalars index) < uint64Modulus := by
  have sumBound := fin8_sum_lt scalarsFit
  have envelope : 8 * ifmaRadix < uint64Modulus := by
    norm_num [ifmaRadix, uint64Modulus]
  exact sumBound.trans envelope

def exactCrtQuotient
    (scalars shiftedFactors : Fin 8 → Nat) (modulus : Nat) : Nat :=
  exactCrtNumerator scalars shiftedFactors / modulus

def estimatedCrtQuotient
    (scalars shiftedFactors : Fin 8 → Nat) (modulus : Nat) : Nat :=
  approximateCrtNumerator scalars shiftedFactors modulus / uint64Modulus

theorem estimated_quotient_bracket
    (scalars shiftedFactors : Fin 8 → Nat) (modulus : Nat)
    (modulusPositive : 0 < modulus)
    (scalarsFit : ∀ index, scalars index < ifmaRadix) :
    exactCrtQuotient scalars shiftedFactors modulus ≤
        estimatedCrtQuotient scalars shiftedFactors modulus ∧
    estimatedCrtQuotient scalars shiftedFactors modulus ≤
        exactCrtQuotient scalars shiftedFactors modulus + 1 := by
  have enclosure := aggregate_quotient_enclosure
    scalars shiftedFactors modulus modulusPositive
  have scalarSumBound := scalar_sum_below_fixed_point scalars scalarsFit
  have exactSplit := Nat.mod_add_div
    (exactCrtNumerator scalars shiftedFactors) modulus
  have exactRemainder := Nat.mod_lt
    (exactCrtNumerator scalars shiftedFactors) modulusPositive
  have approximateSplit := Nat.mod_add_div
    (approximateCrtNumerator scalars shiftedFactors modulus) uint64Modulus
  have wordPositive : 0 < uint64Modulus := by norm_num [uint64Modulus]
  have approximateRemainder := Nat.mod_lt
    (approximateCrtNumerator scalars shiftedFactors modulus) wordPositive
  unfold exactCrtQuotient estimatedCrtQuotient
  constructor
  · by_contra contradiction
    have quotientOrder :
        approximateCrtNumerator scalars shiftedFactors modulus / uint64Modulus <
          exactCrtNumerator scalars shiftedFactors / modulus := by omega
    have quotientStep :
        approximateCrtNumerator scalars shiftedFactors modulus / uint64Modulus + 1 ≤
          exactCrtNumerator scalars shiftedFactors / modulus := by omega
    have approximateUpper :
        approximateCrtNumerator scalars shiftedFactors modulus <
          uint64Modulus *
            (approximateCrtNumerator scalars shiftedFactors modulus /
              uint64Modulus + 1) := by
      calc
        approximateCrtNumerator scalars shiftedFactors modulus =
            approximateCrtNumerator scalars shiftedFactors modulus %
                uint64Modulus +
              uint64Modulus *
                (approximateCrtNumerator scalars shiftedFactors modulus /
                  uint64Modulus) := approximateSplit.symm
        _ < uint64Modulus + uint64Modulus *
              (approximateCrtNumerator scalars shiftedFactors modulus /
                uint64Modulus) :=
              Nat.add_lt_add_right approximateRemainder _
        _ = uint64Modulus *
              (approximateCrtNumerator scalars shiftedFactors modulus /
                uint64Modulus + 1) := by ring
    have exactLower :
        modulus * (exactCrtNumerator scalars shiftedFactors / modulus) ≤
          exactCrtNumerator scalars shiftedFactors := by
      calc
        modulus * (exactCrtNumerator scalars shiftedFactors / modulus) ≤
            exactCrtNumerator scalars shiftedFactors % modulus +
              modulus * (exactCrtNumerator scalars shiftedFactors / modulus) :=
                Nat.le_add_left _ _
        _ = exactCrtNumerator scalars shiftedFactors := exactSplit
    have quotientScaled := Nat.mul_le_mul_left uint64Modulus quotientStep
    have approximateScaled :=
      (Nat.mul_lt_mul_right modulusPositive).2 approximateUpper
    have exactScaled := Nat.mul_le_mul_right uint64Modulus exactLower
    have impossible :
        approximateCrtNumerator scalars shiftedFactors modulus * modulus <
          exactCrtNumerator scalars shiftedFactors * uint64Modulus := by
      calc
        approximateCrtNumerator scalars shiftedFactors modulus * modulus <
            (uint64Modulus *
              (approximateCrtNumerator scalars shiftedFactors modulus /
                uint64Modulus + 1)) * modulus := approximateScaled
        _ ≤ (uint64Modulus *
              (exactCrtNumerator scalars shiftedFactors / modulus)) *
                modulus := Nat.mul_le_mul_right modulus quotientScaled
        _ = (modulus *
              (exactCrtNumerator scalars shiftedFactors / modulus)) *
                uint64Modulus := by ring
        _ ≤ exactCrtNumerator scalars shiftedFactors * uint64Modulus :=
              exactScaled
    exact (not_lt_of_ge enclosure.1) impossible
  · by_contra contradiction
    have quotientOrder :
        exactCrtNumerator scalars shiftedFactors / modulus + 1 <
          approximateCrtNumerator scalars shiftedFactors modulus / uint64Modulus := by
      omega
    have quotientStep :
        exactCrtNumerator scalars shiftedFactors / modulus + 2 ≤
          approximateCrtNumerator scalars shiftedFactors modulus / uint64Modulus := by
      omega
    have approximateLower :
        uint64Modulus *
            (approximateCrtNumerator scalars shiftedFactors modulus /
              uint64Modulus) ≤
          approximateCrtNumerator scalars shiftedFactors modulus := by
      calc
        uint64Modulus *
            (approximateCrtNumerator scalars shiftedFactors modulus /
              uint64Modulus) ≤
          approximateCrtNumerator scalars shiftedFactors modulus %
              uint64Modulus +
            uint64Modulus *
              (approximateCrtNumerator scalars shiftedFactors modulus /
                uint64Modulus) := Nat.le_add_left _ _
        _ = approximateCrtNumerator scalars shiftedFactors modulus :=
              approximateSplit
    have exactUpper :
        exactCrtNumerator scalars shiftedFactors <
          modulus * (exactCrtNumerator scalars shiftedFactors / modulus + 1) := by
      calc
        exactCrtNumerator scalars shiftedFactors =
            exactCrtNumerator scalars shiftedFactors % modulus +
              modulus * (exactCrtNumerator scalars shiftedFactors / modulus) :=
                exactSplit.symm
        _ < modulus +
              modulus * (exactCrtNumerator scalars shiftedFactors / modulus) :=
                Nat.add_lt_add_right exactRemainder _
        _ = modulus *
              (exactCrtNumerator scalars shiftedFactors / modulus + 1) := by ring
    have exactScaled := (Nat.mul_lt_mul_right wordPositive).2 exactUpper
    have modulusPred : modulus - 1 + 1 = modulus := by omega
    have scalarPred :
        (∑ index, scalars index) ≤ uint64Modulus - 1 := by omega
    have errorNonStrict := Nat.mul_le_mul
      (Nat.le_refl (modulus - 1)) scalarPred
    have errorEnvelope :
        (modulus - 1) * (uint64Modulus - 1) <
          modulus * uint64Modulus := by
      have wordPred : uint64Modulus - 1 + 1 = uint64Modulus := by omega
      nlinarith
    have errorBound :
        (modulus - 1) * ∑ index, scalars index <
          modulus * uint64Modulus :=
      errorNonStrict.trans_lt errorEnvelope
    have quotientScaled := Nat.mul_le_mul_left uint64Modulus quotientStep
    have lowerApproximate := Nat.mul_le_mul_right modulus approximateLower
    have lowerQuotient := Nat.mul_le_mul_right modulus quotientScaled
    have correctedUpper :
        exactCrtNumerator scalars shiftedFactors * uint64Modulus +
            (modulus - 1) * ∑ index, scalars index <
          (uint64Modulus *
              (exactCrtNumerator scalars shiftedFactors / modulus + 2)) *
            modulus := by
      calc
        exactCrtNumerator scalars shiftedFactors * uint64Modulus +
            (modulus - 1) * ∑ index, scalars index <
          (modulus *
              (exactCrtNumerator scalars shiftedFactors / modulus + 1)) *
                uint64Modulus + modulus * uint64Modulus :=
            Nat.add_lt_add exactScaled errorBound
        _ = (uint64Modulus *
              (exactCrtNumerator scalars shiftedFactors / modulus + 2)) *
                modulus := by ring
    have impossible :
        exactCrtNumerator scalars shiftedFactors * uint64Modulus +
            (modulus - 1) * ∑ index, scalars index <
          approximateCrtNumerator scalars shiftedFactors modulus * modulus := by
      exact correctedUpper.trans_le (lowerQuotient.trans lowerApproximate)
    exact (not_lt_of_ge enclosure.2) impossible

theorem estimated_quotient_exact_of_margin
    (scalars shiftedFactors : Fin 8 → Nat) (modulus : Nat)
    (modulusPositive : 0 < modulus)
    (scalarsFit : ∀ index, scalars index < ifmaRadix)
    (margin :
      (exactCrtNumerator scalars shiftedFactors % modulus) * uint64Modulus +
          (modulus - 1) * ∑ index, scalars index <
        modulus * uint64Modulus) :
    estimatedCrtQuotient scalars shiftedFactors modulus =
      exactCrtQuotient scalars shiftedFactors modulus := by
  have bracket := estimated_quotient_bracket
    scalars shiftedFactors modulus modulusPositive scalarsFit
  have enclosure := aggregate_quotient_enclosure
    scalars shiftedFactors modulus modulusPositive
  have exactSplit := Nat.mod_add_div
    (exactCrtNumerator scalars shiftedFactors) modulus
  have exactDecomposition :
      exactCrtNumerator scalars shiftedFactors =
        exactCrtQuotient scalars shiftedFactors modulus * modulus +
          exactCrtNumerator scalars shiftedFactors % modulus := by
    unfold exactCrtQuotient
    nlinarith
  have numeratorUpper :
      approximateCrtNumerator scalars shiftedFactors modulus * modulus <
        (exactCrtQuotient scalars shiftedFactors modulus + 1) *
          uint64Modulus * modulus := by
    calc
      approximateCrtNumerator scalars shiftedFactors modulus * modulus ≤
          exactCrtNumerator scalars shiftedFactors * uint64Modulus +
            (modulus - 1) * ∑ index, scalars index := enclosure.2
      _ = (exactCrtQuotient scalars shiftedFactors modulus * modulus +
            exactCrtNumerator scalars shiftedFactors % modulus) *
              uint64Modulus +
            (modulus - 1) * ∑ index, scalars index := by
          nth_rewrite 1 [exactDecomposition]
          rfl
      _ < (exactCrtQuotient scalars shiftedFactors modulus * modulus) *
              uint64Modulus + modulus * uint64Modulus := by
          nlinarith
      _ = (exactCrtQuotient scalars shiftedFactors modulus + 1) *
            uint64Modulus * modulus := by ring
  have approximateUpper :
      approximateCrtNumerator scalars shiftedFactors modulus <
        (exactCrtQuotient scalars shiftedFactors modulus + 1) *
          uint64Modulus := by
    exact (Nat.mul_lt_mul_right modulusPositive).mp (by
      simpa [Nat.mul_assoc, Nat.mul_comm, Nat.mul_left_comm] using numeratorUpper)
  have quotientUpper :
      estimatedCrtQuotient scalars shiftedFactors modulus <
        exactCrtQuotient scalars shiftedFactors modulus + 1 := by
    unfold estimatedCrtQuotient
    exact (Nat.div_lt_iff_lt_mul (by norm_num [uint64Modulus])).2
      (by simpa [Nat.mul_assoc, Nat.mul_comm] using approximateUpper)
  omega

/-- If the aggregate upward-rounding error covers the distance from the exact
CRT remainder to the next multiple of `modulus`, the fixed-point quotient is
the upper member of the generic two-value bracket. -/
theorem estimated_quotient_one_high_of_margin
    (scalars shiftedFactors : Fin 8 → Nat) (modulus : Nat)
    (modulusPositive : 0 < modulus)
    (scalarsFit : ∀ index, scalars index < ifmaRadix)
    (margin :
      (modulus - exactCrtNumerator scalars shiftedFactors % modulus) *
          uint64Modulus ≤
        ∑ index, scalars index *
          fixedPointCeilingError (shiftedFactors index) modulus) :
    estimatedCrtQuotient scalars shiftedFactors modulus =
      exactCrtQuotient scalars shiftedFactors modulus + 1 := by
  have bracket := estimated_quotient_bracket
    scalars shiftedFactors modulus modulusPositive scalarsFit
  have exactSplit := Nat.mod_add_div
    (exactCrtNumerator scalars shiftedFactors) modulus
  have exactRemainder := Nat.mod_lt
    (exactCrtNumerator scalars shiftedFactors) modulusPositive
  have exactDecomposition :
      exactCrtNumerator scalars shiftedFactors =
        exactCrtQuotient scalars shiftedFactors modulus * modulus +
          exactCrtNumerator scalars shiftedFactors % modulus := by
    unfold exactCrtQuotient
    nlinarith
  have errorIdentity := aggregate_ceiling_error_identity
    scalars shiftedFactors modulus modulusPositive
  have remainderSplit :
      exactCrtNumerator scalars shiftedFactors % modulus +
          (modulus - exactCrtNumerator scalars shiftedFactors % modulus) =
        modulus := by
    exact Nat.add_sub_of_le exactRemainder.le
  have scaledLower :
      (exactCrtQuotient scalars shiftedFactors modulus + 1) *
          uint64Modulus * modulus ≤
        approximateCrtNumerator scalars shiftedFactors modulus * modulus := by
    rw [← errorIdentity]
    rw [exactDecomposition]
    nlinarith
  have approximateLower :
      (exactCrtQuotient scalars shiftedFactors modulus + 1) *
          uint64Modulus ≤
        approximateCrtNumerator scalars shiftedFactors modulus := by
    nlinarith
  have approximateSplit := Nat.mod_add_div
    (approximateCrtNumerator scalars shiftedFactors modulus) uint64Modulus
  have wordPositive : 0 < uint64Modulus := by norm_num [uint64Modulus]
  have approximateRemainder := Nat.mod_lt
    (approximateCrtNumerator scalars shiftedFactors modulus) wordPositive
  have quotientLower :
      exactCrtQuotient scalars shiftedFactors modulus + 1 ≤
        estimatedCrtQuotient scalars shiftedFactors modulus := by
    unfold estimatedCrtQuotient
    by_contra contradiction
    have quotientOrder :
        approximateCrtNumerator scalars shiftedFactors modulus /
            uint64Modulus <
          exactCrtQuotient scalars shiftedFactors modulus + 1 := by omega
    have approximateUpper :
        approximateCrtNumerator scalars shiftedFactors modulus <
          (exactCrtQuotient scalars shiftedFactors modulus + 1) *
            uint64Modulus := by
      calc
        approximateCrtNumerator scalars shiftedFactors modulus =
            approximateCrtNumerator scalars shiftedFactors modulus %
                uint64Modulus +
              uint64Modulus *
                (approximateCrtNumerator scalars shiftedFactors modulus /
                  uint64Modulus) := approximateSplit.symm
        _ < uint64Modulus + uint64Modulus *
              (approximateCrtNumerator scalars shiftedFactors modulus /
                uint64Modulus) :=
              Nat.add_lt_add_right approximateRemainder _
        _ ≤ uint64Modulus + uint64Modulus *
              exactCrtQuotient scalars shiftedFactors modulus := by
              nlinarith
        _ = (exactCrtQuotient scalars shiftedFactors modulus + 1) *
              uint64Modulus := by ring
    exact (not_lt_of_ge approximateLower) approximateUpper
  omega

def correctedCrtRepresentative
    (scalars shiftedFactors : Fin 8 → Nat) (modulus : Nat) : Int :=
  (exactCrtNumerator scalars shiftedFactors : Int) -
    (estimatedCrtQuotient scalars shiftedFactors modulus : Int) * modulus

theorem corrected_representative_range
    (scalars shiftedFactors : Fin 8 → Nat) (modulus : Nat)
    (modulusPositive : 0 < modulus)
    (scalarsFit : ∀ index, scalars index < ifmaRadix) :
    -(modulus : Int) ≤
        correctedCrtRepresentative scalars shiftedFactors modulus ∧
      correctedCrtRepresentative scalars shiftedFactors modulus < modulus := by
  have bracket := estimated_quotient_bracket
    scalars shiftedFactors modulus modulusPositive scalarsFit
  have quotientCases :
      estimatedCrtQuotient scalars shiftedFactors modulus =
          exactCrtQuotient scalars shiftedFactors modulus ∨
      estimatedCrtQuotient scalars shiftedFactors modulus =
          exactCrtQuotient scalars shiftedFactors modulus + 1 := by
    omega
  have exactSplit := Nat.mod_add_div
    (exactCrtNumerator scalars shiftedFactors) modulus
  have exactRemainder := Nat.mod_lt
    (exactCrtNumerator scalars shiftedFactors) modulusPositive
  have exactSplitInt :
      (exactCrtNumerator scalars shiftedFactors : Int) =
        (exactCrtNumerator scalars shiftedFactors % modulus : Int) +
          (modulus : Int) *
            (exactCrtQuotient scalars shiftedFactors modulus : Int) := by
    unfold exactCrtQuotient
    exact_mod_cast exactSplit.symm
  rcases quotientCases with exactQuotient | oneHigh
  · unfold correctedCrtRepresentative
    rw [exactQuotient, exactSplitInt]
    constructor <;> ring_nf <;> omega
  · unfold correctedCrtRepresentative
    rw [oneHigh, exactSplitInt]
    constructor <;> push_cast <;> ring_nf <;> omega

end CatBLST.AVX.MatrixNoK
