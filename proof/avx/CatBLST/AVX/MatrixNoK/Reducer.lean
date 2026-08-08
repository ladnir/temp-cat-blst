/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Montgomery
import Mathlib.Tactic

/-!
# MatrixNoK selected reducer bounds

This file models the exact `reduce_wide` shape used by the fixed MatrixNoK
reducer.  A stored pair represents `high * 2^52 + low`; folding the excess low
bits into `high` is exact.  The selected BLS12-381 moduli lie just below
`2^50`, so a wide bound of `A * m²` becomes a high-word bound of
`ceil(A / 4) * m` before the one-lane Montgomery step.
-/

namespace CatBLST.AVX.MatrixNoK

def foldedWideHigh (high low : Nat) : Nat := high + low / ifmaRadix

def selectedMontgomeryStep
    (high low modulus factor : Nat) : Nat :=
  montgomeryReduceStep (foldedWideHigh high low) low modulus factor

def selectedReduceMasked
    (high low modulus factor : Nat) : Nat :=
  maskReduceValue (selectedMontgomeryStep high low modulus factor) modulus

theorem folded_high_scaled_le_wide (high low : Nat) :
    foldedWideHigh high low * ifmaRadix ≤ high * ifmaRadix + low := by
  unfold foldedWideHigh
  have division := Nat.div_mul_le_self low ifmaRadix
  nlinarith

theorem folded_high_lt_quarter_bound
    (high low modulus wideBound quarterBound : Nat)
    (modulusPositive : 0 < modulus)
    (modulusWindow : 4 * modulus < ifmaRadix)
    (quarterPositive : 0 < quarterBound)
    (wideFits : high * ifmaRadix + low ≤ wideBound * modulus * modulus)
    (wideBoundFits : wideBound ≤ 4 * quarterBound) :
    foldedWideHigh high low < quarterBound * modulus := by
  have foldedScaled := folded_high_scaled_le_wide high low
  by_contra contradiction
  have foldedLower : quarterBound * modulus ≤ foldedWideHigh high low := by
    omega
  have lowerScaled := Nat.mul_le_mul_right ifmaRadix foldedLower
  have strictEnvelope :
      (4 * quarterBound) * modulus * modulus <
        (quarterBound * modulus) * ifmaRadix := by
    have scaled := (Nat.mul_lt_mul_left
      (Nat.mul_pos quarterPositive modulusPositive)).2 modulusWindow
    simpa [Nat.mul_assoc, Nat.mul_comm, Nat.mul_left_comm] using scaled
  have wideEnvelope :
      wideBound * modulus * modulus ≤
        (4 * quarterBound) * modulus * modulus := by
    simpa [Nat.mul_assoc] using
      (Nat.mul_le_mul_right (modulus * modulus) wideBoundFits)
  have wideStrict :
      high * ifmaRadix + low <
        (quarterBound * modulus) * ifmaRadix := by
    calc
      high * ifmaRadix + low ≤ wideBound * modulus * modulus := wideFits
      _ ≤ (4 * quarterBound) * modulus * modulus := wideEnvelope
      _ < (quarterBound * modulus) * ifmaRadix := strictEnvelope
  have foldedStrict :
      foldedWideHigh high low * ifmaRadix <
        (quarterBound * modulus) * ifmaRadix :=
    foldedScaled.trans_lt wideStrict
  have impossible :
      (quarterBound * modulus) * ifmaRadix <
        (quarterBound * modulus) * ifmaRadix :=
    lowerScaled.trans_lt foldedStrict
  exact (Nat.lt_irrefl _) impossible

theorem selected_modulus_quarter_window
    (modulus : Nat) (modulusMember : modulus ∈ allModuli) :
    0 < modulus ∧ 4 * modulus < ifmaRadix := by
  have window := all_moduli_element_window modulus modulusMember
  constructor
  · have radixPositive : 0 < elementRadix - 113 := by
      norm_num [elementRadix]
    omega
  · have radixRelation : 4 * elementRadix = ifmaRadix := by
      norm_num [elementRadix, ifmaRadix]
    nlinarith

theorem selected_step_congruent_wide
    (high low modulus factor : Nat)
    (modulusPositive : 0 < modulus)
    (inverseRelation : modulus * factor % ifmaRadix = 1) :
    (selectedMontgomeryStep high low modulus factor * ifmaRadix) % modulus =
      (high * ifmaRadix + low) % modulus := by
  unfold selectedMontgomeryStep foldedWideHigh
  exact reduce_wide_step_congruent high low modulus factor
    modulusPositive inverseRelation

theorem selected_step_bound
    (high low modulus factor wideBound quarterBound : Nat)
    (modulusPositive : 0 < modulus)
    (modulusWindow : 4 * modulus < ifmaRadix)
    (quarterPositive : 0 < quarterBound)
    (wideFits : high * ifmaRadix + low ≤ wideBound * modulus * modulus)
    (wideBoundFits : wideBound ≤ 4 * quarterBound) :
    selectedMontgomeryStep high low modulus factor <
      (quarterBound + 1) * modulus := by
  unfold selectedMontgomeryStep
  apply montgomery_step_bound _ _ modulus factor quarterBound
    modulusPositive
  exact folded_high_lt_quarter_bound high low modulus wideBound quarterBound
    modulusPositive modulusWindow quarterPositive wideFits wideBoundFits

theorem selected_masked_reducer_contract
    (high low modulus factor wideBound quarterBound : Nat)
    (modulusMember : modulus ∈ allModuli)
    (inverseRelation : modulus * factor % ifmaRadix = 1)
    (quarterPositive : 0 < quarterBound)
    (wideFits : high * ifmaRadix + low ≤ wideBound * modulus * modulus)
    (wideBoundFits : wideBound ≤ 4 * quarterBound)
    (outputFitsWord : (quarterBound + 1) * elementRadix < uint64Modulus) :
    (selectedReduceMasked high low modulus factor : ZMod modulus) *
        ifmaRadix = high * ifmaRadix + low ∧
      selectedReduceMasked high low modulus factor < 2 * modulus := by
  have modulusFacts := selected_modulus_quarter_window modulus modulusMember
  have stepBound := selected_step_bound high low modulus factor
    wideBound quarterBound modulusFacts.1 modulusFacts.2
    quarterPositive wideFits wideBoundFits
  have stepFitsWord :
      selectedMontgomeryStep high low modulus factor < uint64Modulus := by
    have modulusBelow := (all_moduli_element_window modulus modulusMember).2
    have scaledBelow :
        (quarterBound + 1) * modulus <
          (quarterBound + 1) * elementRadix := by
      exact (Nat.mul_lt_mul_left (by omega)).2 modulusBelow
    exact stepBound.trans (scaledBelow.trans outputFitsWord)
  have masked := concrete_mask_reduce_contract stepFitsWord modulusMember
  constructor
  · have maskedCast :
        (selectedReduceMasked high low modulus factor : ZMod modulus) =
          (selectedMontgomeryStep high low modulus factor : ZMod modulus) := by
      rw [ZMod.natCast_eq_natCast_iff]
      exact masked.1
    rw [maskedCast]
    have stepCongruent := selected_step_congruent_wide high low modulus factor
      modulusFacts.1 inverseRelation
    have stepCast :
        (selectedMontgomeryStep high low modulus factor * ifmaRadix :
            ZMod modulus) =
          (high * ifmaRadix + low : ZMod modulus) := by
      have castEquality :=
        (ZMod.natCast_eq_natCast_iff
          (selectedMontgomeryStep high low modulus factor * ifmaRadix)
          (high * ifmaRadix + low) modulus).2 stepCongruent
      simpa only [Nat.cast_mul, Nat.cast_add] using castEquality
    exact stepCast
  · exact masked.2

theorem selected_product_step_bound
    (high low modulus factor : Nat)
    (modulusMember : modulus ∈ allModuli)
    (wideFits : high * ifmaRadix + low ≤ 4 * modulus * modulus) :
    selectedMontgomeryStep high low modulus factor < 2 * modulus := by
  have facts := selected_modulus_quarter_window modulus modulusMember
  exact selected_step_bound high low modulus factor 4 1 facts.1 facts.2
    (by norm_num) wideFits (by norm_num)

theorem selected_combined_reducer_contract
    (high low modulus factor : Nat)
    (modulusMember : modulus ∈ allModuli)
    (inverseRelation : modulus * factor % ifmaRadix = 1)
    (wideFits : high * ifmaRadix + low ≤ 20 * modulus * modulus) :
    (selectedReduceMasked high low modulus factor : ZMod modulus) *
        ifmaRadix = high * ifmaRadix + low ∧
      selectedReduceMasked high low modulus factor < 2 * modulus := by
  exact selected_masked_reducer_contract high low modulus factor 20 5
    modulusMember inverseRelation (by norm_num) wideFits (by norm_num)
    (by norm_num [elementRadix, uint64Modulus])

theorem selected_expansion_reducer_contract
    (high low modulus factor : Nat)
    (modulusMember : modulus ∈ allModuli)
    (inverseRelation : modulus * factor % ifmaRadix = 1)
    (wideFits : high * ifmaRadix + low ≤ 32 * modulus * modulus) :
    (selectedReduceMasked high low modulus factor : ZMod modulus) *
        ifmaRadix = high * ifmaRadix + low ∧
      selectedReduceMasked high low modulus factor < 2 * modulus := by
  exact selected_masked_reducer_contract high low modulus factor 32 8
    modulusMember inverseRelation (by norm_num) wideFits (by norm_num)
    (by norm_num [elementRadix, uint64Modulus])

/--
The exact reducer selected by the generic delayed path in
`FpRingCandidate::batch_reduce_expand`.  C++ assigns the accumulated M2 value
the conservative bound `MAX_ADD + delayed_r1_type_bound = 800 + 32 = 832`.
`reduce_wide` therefore returns a value below `209m`; the following
`reduce_auto<2>` takes its mask-reduction branch and produces the value modeled
by `selectedReduceMasked`.
-/
theorem selected_delayed_reducer_contract
    (high low modulus factor : Nat)
    (modulusMember : modulus ∈ allModuli)
    (inverseRelation : modulus * factor % ifmaRadix = 1)
    (wideFits : high * ifmaRadix + low ≤
      delayedCombinedTypeBound * modulus * modulus) :
    (selectedReduceMasked high low modulus factor : ZMod modulus) *
        ifmaRadix = high * ifmaRadix + low ∧
      selectedMontgomeryStep high low modulus factor <
        reducerOutputBound * modulus ∧
      selectedReduceMasked high low modulus factor < 2 * modulus := by
  have modulusFacts := selected_modulus_quarter_window modulus modulusMember
  have stepBound := selected_step_bound high low modulus factor
    delayedCombinedTypeBound reducerReboundedHighBound
    modulusFacts.1 modulusFacts.2
    (by norm_num [reducerReboundedHighBound, delayedCombinedTypeBound,
      delayedInputBound, delayedR1TypeBound, r1ExactAccumulationBound,
      standardBound])
    wideFits
    (by norm_num [reducerReboundedHighBound, delayedCombinedTypeBound,
      delayedInputBound, delayedR1TypeBound, r1ExactAccumulationBound,
      standardBound])
  have masked := selected_masked_reducer_contract high low modulus factor
    delayedCombinedTypeBound reducerReboundedHighBound modulusMember
    inverseRelation
    (by norm_num [reducerReboundedHighBound, delayedCombinedTypeBound,
      delayedInputBound, delayedR1TypeBound, r1ExactAccumulationBound,
      standardBound])
    wideFits
    (by norm_num [reducerReboundedHighBound, delayedCombinedTypeBound,
      delayedInputBound, delayedR1TypeBound, r1ExactAccumulationBound,
      standardBound])
    (by norm_num [reducerReboundedHighBound, delayedCombinedTypeBound,
      delayedInputBound, delayedR1TypeBound, r1ExactAccumulationBound,
      standardBound, elementRadix, uint64Modulus])
  refine ⟨masked.1, ?_, masked.2⟩
  simpa [reducerOutputBound] using stepBound

end CatBLST.AVX.MatrixNoK
