/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Regeneration
import Mathlib.Tactic

/-!
# MatrixNoK selected-kernel bounds

This file records the numeric envelopes used by the selected BLS12-381
MatrixNoK kernels.  It deliberately distinguishes the exact multiplication
bound (`4 + 16 = 20`) from the generic delayed-path C++ type allowance
(`800 + 32 = 832`).  The latter is conservative and is not silently narrowed
by the proof.
-/

namespace CatBLST.AVX.MatrixNoK

open scoped BigOperators

def elementRadix : Nat := 2 ^ 50
def ifmaRadix : Nat := 2 ^ 52
def uint64Modulus : Nat := 2 ^ 64
def uint128Modulus : Nat := 2 ^ 128

def standardBound : Nat := 2
def multiplicationProductBound : Nat := standardBound * standardBound
def r1ExactAccumulationBound : Nat := 8 * standardBound
def delayedR1TypeBound : Nat := 2 * r1ExactAccumulationBound
def multiplicationCombinedBound : Nat :=
  multiplicationProductBound + r1ExactAccumulationBound
def delayedInputBound : Nat := 800
def delayedCombinedTypeBound : Nat := delayedInputBound + delayedR1TypeBound
def reducerReboundedHighBound : Nat :=
  (delayedCombinedTypeBound + 3) / 4
def reducerOutputBound : Nat := reducerReboundedHighBound + 1

/-!
## Completed RNS range

VROOM's BLS12-381 schedule tracks represented multiples of `p` in
`[-1932, 2377]`.  The C++ completion step adds `1932 * p` whenever the lower
bound is negative, translating the complete reviewed interval to `[0, 4309]`.
The final theorem below restates VROOM's conservative two-basis capacity
check for the specialized CatBLST constants.
-/

def minTrackedRnsMultiple : Int := -1932
def maxTrackedRnsMultiple : Int := 2377
def negativeCompletionMultiple : Int := -minTrackedRnsMultiple
def maxCompletedRnsMultiple : Int :=
  maxTrackedRnsMultiple - minTrackedRnsMultiple
def baseRedundantMultiple : Nat := 40
def wideRedundantMultiple : Nat :=
  baseRedundantMultiple * baseRedundantMultiple
def rnsBits : Nat := 50 * 8 - 1

theorem completed_rns_range_values :
    negativeCompletionMultiple = 1932 ∧
    maxCompletedRnsMultiple = 4309 ∧
    wideRedundantMultiple = 1600 ∧
    rnsBits = 399 := by
  norm_num [negativeCompletionMultiple, minTrackedRnsMultiple,
    maxCompletedRnsMultiple, maxTrackedRnsMultiple,
    wideRedundantMultiple, baseRedundantMultiple, rnsBits]

/-- The entire completed range, including its wide redundant encoding, fits
strictly inside the exact product of the two generated RNS bases. -/
theorem completed_rns_range_fits_combined_bases :
    (4309 : Nat) * wideRedundantMultiple *
        Generated.target_modulus * Generated.target_modulus <
      Generated.M1 * Generated.M2 := by
  norm_num [wideRedundantMultiple, baseRedundantMultiple,
    Generated.target_modulus, Generated.M1, Generated.M2]

theorem selected_bound_values :
    multiplicationProductBound = 4 ∧
    r1ExactAccumulationBound = 16 ∧
    delayedR1TypeBound = 32 ∧
    multiplicationCombinedBound = 20 ∧
    delayedCombinedTypeBound = 832 ∧
    reducerReboundedHighBound = 208 ∧
    reducerOutputBound = 209 := by
  norm_num [multiplicationProductBound, r1ExactAccumulationBound,
    delayedR1TypeBound, multiplicationCombinedBound, delayedInputBound,
    delayedCombinedTypeBound, reducerReboundedHighBound, reducerOutputBound,
    standardBound]

theorem multiplication_combined_below_generated_base :
    multiplicationCombinedBound < 40 := by
  norm_num [multiplicationCombinedBound, multiplicationProductBound,
    r1ExactAccumulationBound, standardBound]

theorem all_moduli_element_window :
    ∀ modulus ∈ allModuli,
      elementRadix - 113 ≤ modulus ∧ modulus < elementRadix := by
  norm_num [allModuli, Generated.moduli1, Generated.moduli2, elementRadix]

theorem all_moduli_ifma_compatible :
    ∀ modulus ∈ allModuli, 2 * modulus < ifmaRadix := by
  norm_num [allModuli, Generated.moduli1, Generated.moduli2, ifmaRadix]

def maskReduceValue (value modulus : Nat) : Nat :=
  value % elementRadix +
    value / elementRadix * (elementRadix - modulus)

theorem mask_reduce_reconstruction {value modulus : Nat}
    (modulusAtMostRadix : modulus ≤ elementRadix) :
    maskReduceValue value modulus +
        (value / elementRadix) * modulus = value := by
  have split := Nat.mod_add_div value elementRadix
  unfold maskReduceValue
  have radixSplit : elementRadix - modulus + modulus = elementRadix :=
    Nat.sub_add_cancel modulusAtMostRadix
  calc
    value % elementRadix + value / elementRadix * (elementRadix - modulus) +
        value / elementRadix * modulus =
        value % elementRadix + value / elementRadix *
          ((elementRadix - modulus) + modulus) := by ring
    _ = value % elementRadix + value / elementRadix * elementRadix := by
      rw [radixSplit]
    _ = value := by simpa [Nat.mul_comm] using split

theorem mask_reduce_congruent {value modulus : Nat}
    (modulusAtMostRadix : modulus ≤ elementRadix) :
    maskReduceValue value modulus % modulus = value % modulus := by
  have reconstruction := mask_reduce_reconstruction (value := value)
    modulusAtMostRadix
  have reduced := congrArg (fun result => result % modulus) reconstruction
  simpa [Nat.add_mod] using reduced

theorem mask_reduce_below_two_moduli {value modulus : Nat}
    (valueFitsWord : value < uint64Modulus)
    (modulusWindow : elementRadix - 113 ≤ modulus)
    (modulusAtMostRadix : modulus ≤ elementRadix) :
    maskReduceValue value modulus < 2 * modulus := by
  have radixPositive : 0 < elementRadix := by norm_num [elementRadix]
  have remainderBound := Nat.mod_lt value radixPositive
  have quotientBound : value / elementRadix < 2 ^ 14 := by
    apply (Nat.div_lt_iff_lt_mul radixPositive).2
    norm_num [uint64Modulus, elementRadix] at valueFitsWord ⊢
    exact valueFitsWord
  have deltaBound : elementRadix - modulus ≤ 113 := by omega
  have quotientBoundNonStrict : value / elementRadix ≤ 2 ^ 14 - 1 := by omega
  have productBound := Nat.mul_le_mul quotientBoundNonStrict deltaBound
  have numericGap :
      elementRadix + (2 ^ 14 - 1) * 113 <
        2 * (elementRadix - 113) := by
    norm_num [elementRadix]
  unfold maskReduceValue
  nlinarith

theorem concrete_mask_reduce_contract {value modulus : Nat}
    (valueFitsWord : value < uint64Modulus)
    (modulusMember : modulus ∈ allModuli) :
    maskReduceValue value modulus % modulus = value % modulus ∧
      maskReduceValue value modulus < 2 * modulus := by
  have window := all_moduli_element_window modulus modulusMember
  constructor
  · exact mask_reduce_congruent (le_of_lt window.2)
  · exact mask_reduce_below_two_moduli valueFitsWord window.1 (le_of_lt window.2)

theorem fin8_sum_lt {values : Fin 8 → Nat} {bound : Nat}
    (eachFits : ∀ index, values index < bound) :
    (∑ index, values index) < 8 * bound := by
  calc
    (∑ index, values index) ≤ ∑ _index : Fin 8, (bound - 1) := by
      apply Finset.sum_le_sum
      intro index _
      have := eachFits index
      omega
    _ = 8 * (bound - 1) := by simp
    _ < 8 * bound := by
      have boundPositive : 0 < bound := by
        have := eachFits (0 : Fin 8)
        omega
      omega

/-- Worst selected delayed accumulator plus eight IFMA digits cannot wrap. -/
theorem selected_eight_term_accumulator_fits_word
    (initial : Nat) (contributions : Fin 8 → Nat)
    (initialFits : initial < delayedInputBound * elementRadix)
    (contributionsFit : ∀ index, contributions index < ifmaRadix) :
    initial + ∑ index, contributions index < uint64Modulus := by
  have contributionSum := fin8_sum_lt contributionsFit
  have numericEnvelope :
      delayedInputBound * elementRadix + 8 * ifmaRadix < uint64Modulus := by
    norm_num [delayedInputBound, elementRadix, ifmaRadix, uint64Modulus]
  omega

def quotientNumerator
    (scalars estimates : Fin 8 → Nat) : Nat :=
  ∑ index, scalars index * estimates index

theorem quotient_numerator_fits_119_bits
    (scalars estimates : Fin 8 → Nat)
    (scalarsFit : ∀ index, scalars index < ifmaRadix)
    (estimatesFit : ∀ index, estimates index < uint64Modulus) :
    quotientNumerator scalars estimates < 2 ^ 119 := by
  unfold quotientNumerator
  calc
    (∑ index, scalars index * estimates index) ≤
        ∑ _index : Fin 8, (ifmaRadix - 1) * (uint64Modulus - 1) := by
      apply Finset.sum_le_sum
      intro index _
      have scalarFits := scalarsFit index
      have estimateFits := estimatesFit index
      apply Nat.mul_le_mul <;> omega
    _ = 8 * (ifmaRadix - 1) * (uint64Modulus - 1) := by
      simp
      ring
    _ < 2 ^ 119 := by
      norm_num [ifmaRadix, uint64Modulus]

theorem quotient_words_fit
    (scalars estimates : Fin 8 → Nat)
    (scalarsFit : ∀ index, scalars index < ifmaRadix)
    (estimatesFit : ∀ index, estimates index < uint64Modulus) :
    quotientNumerator scalars estimates < uint128Modulus ∧
    quotientNumerator scalars estimates / uint64Modulus < 2 ^ 55 ∧
    (quotientNumerator scalars estimates / uint64Modulus) / ifmaRadix < 8 := by
  have numeratorBound := quotient_numerator_fits_119_bits
    scalars estimates scalarsFit estimatesFit
  have uint128Bound : quotientNumerator scalars estimates < uint128Modulus := by
    exact numeratorBound.trans (by norm_num [uint128Modulus])
  have quotientBound :
      quotientNumerator scalars estimates / uint64Modulus < 2 ^ 55 := by
    apply (Nat.div_lt_iff_lt_mul (by norm_num [uint64Modulus])).2
    norm_num [uint64Modulus]
    exact numeratorBound
  have highBound :
      (quotientNumerator scalars estimates / uint64Modulus) /
          ifmaRadix < 8 := by
    apply (Nat.div_lt_iff_lt_mul (by norm_num [ifmaRadix])).2
    norm_num [ifmaRadix]
    exact quotientBound
  exact ⟨uint128Bound, quotientBound, highBound⟩

def lowDigitProduct (left right : Nat) : Nat :=
  (left % ifmaRadix) * (right % ifmaRadix) % ifmaRadix

def highDigitProduct (left right : Nat) : Nat :=
  (left % ifmaRadix) * (right % ifmaRadix) / ifmaRadix

theorem digit_products_fit
    (left right : Nat) :
    lowDigitProduct left right < ifmaRadix ∧
      highDigitProduct left right < ifmaRadix := by
  have radixPositive : 0 < ifmaRadix := by norm_num [ifmaRadix]
  have leftDigit := Nat.mod_lt left radixPositive
  have rightDigit := Nat.mod_lt right radixPositive
  have productBound := Nat.mul_le_mul
    (show left % ifmaRadix ≤ ifmaRadix - 1 by omega)
    (show right % ifmaRadix ≤ ifmaRadix - 1 by omega)
  have predecessor : ifmaRadix - 1 + 1 = ifmaRadix := by omega
  have productStrict :
      (ifmaRadix - 1) * (ifmaRadix - 1) < ifmaRadix * ifmaRadix := by
    nlinarith
  constructor
  · exact Nat.mod_lt _ radixPositive
  · apply (Nat.div_lt_iff_lt_mul radixPositive).2
    exact productBound.trans_lt productStrict

/-- Recombine the two radix-`2^52` digits emitted by one IFMA product. -/
theorem digit_product_reconstruction (left right : Nat) :
    lowDigitProduct left right + ifmaRadix * highDigitProduct left right =
      (left % ifmaRadix) * (right % ifmaRadix) := by
  simpa [lowDigitProduct, highDigitProduct, Nat.mul_comm] using
    Nat.mod_add_div
      ((left % ifmaRadix) * (right % ifmaRadix)) ifmaRadix

/--
The non-shifted fixed-matrix correction is the larger of the two correction
schedules: its upper quotient word is reconstructed in the high accumulator.
This theorem therefore also covers the shifted schedule, whose final two IFMA
digits are each individually below `2^52`.
-/
theorem fixed_matrix_correction_fits_word
    (initialHigh initialLow correction quotient : Nat)
    (initialHighFits : initialHigh < 8 * ifmaRadix)
    (initialLowFits : initialLow < 8 * ifmaRadix)
    (correctionFits : correction < ifmaRadix)
    (quotientFits : quotient < 2 ^ 55) :
    initialLow + lowDigitProduct correction quotient < uint64Modulus ∧
    initialHigh + highDigitProduct correction quotient +
        correction * (quotient / ifmaRadix) < uint64Modulus := by
  have productsFit := digit_products_fit correction quotient
  have quotientHighFits : quotient / ifmaRadix < 8 := by
    apply (Nat.div_lt_iff_lt_mul (by norm_num [ifmaRadix])).2
    norm_num [ifmaRadix]
    exact quotientFits
  have upperContribution : correction * (quotient / ifmaRadix) <
      ifmaRadix * 8 := by
    nlinarith
  have lowEnvelope : 9 * ifmaRadix < uint64Modulus := by
    norm_num [ifmaRadix, uint64Modulus]
  have highEnvelope : 17 * ifmaRadix < uint64Modulus := by
    norm_num [ifmaRadix, uint64Modulus]
  constructor <;> nlinarith

/--
The shifted fixed-matrix schedule adds two independent IFMA products to each
accumulator: the low quotient digit times `correction`, then the upper quotient
digit times `correctionShift`.  This is the exact source shape used by
`FixedRNSMatrix::reduce_raw_shifted` and its fixed-batch peer.
-/
theorem shifted_matrix_correction_fits_word
    (initialHigh initialLow correction correctionShift quotient : Nat)
    (initialHighFits : initialHigh < 8 * ifmaRadix)
    (initialLowFits : initialLow < 8 * ifmaRadix) :
    initialLow + lowDigitProduct correction quotient +
          lowDigitProduct correctionShift (quotient / ifmaRadix) <
        uint64Modulus ∧
      initialHigh + highDigitProduct correction quotient +
          highDigitProduct correctionShift (quotient / ifmaRadix) <
        uint64Modulus := by
  have firstProducts := digit_products_fit correction quotient
  have shiftedProducts :=
    digit_products_fit correctionShift (quotient / ifmaRadix)
  have envelope : 10 * ifmaRadix < uint64Modulus := by
    norm_num [ifmaRadix, uint64Modulus]
  constructor <;> omega

theorem exact_multiplication_r1_accumulator_fits_word
    (initial : Nat) (contributions : Fin 8 → Nat)
    (initialFits : initial < ifmaRadix)
    (contributionsFit : ∀ index, contributions index < ifmaRadix) :
    initial + ∑ index, contributions index < uint64Modulus := by
  have contributionSum := fin8_sum_lt contributionsFit
  have envelope : 9 * ifmaRadix < uint64Modulus := by
    norm_num [ifmaRadix, uint64Modulus]
  omega

/-- The largest post-Montgomery bound retained by the selected reducer fits. -/
theorem reducer_output_fits_word :
    reducerOutputBound * elementRadix < uint64Modulus := by
  norm_num [reducerOutputBound, reducerReboundedHighBound,
    delayedCombinedTypeBound, delayedInputBound, delayedR1TypeBound,
    r1ExactAccumulationBound, standardBound, elementRadix, uint64Modulus]

end CatBLST.AVX.MatrixNoK
