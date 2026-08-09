/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.GtExponentiationPhysical
import Mathlib.Tactic

set_option maxHeartbeats 800000

/-!
# Balanced GT scalar split

This begins the remaining compiled-FFT P3 row.  It models the four rounded
divisions in `gt_split_scalar_4d`, the final `z^4 = z^2 - 1 (mod r)` carry
fold, and the odd-coordinate sign change needed for the un-conjugated
Frobenius basis.  Schedule compilation and bytecode execution are separate
obligations built on this split.
-/

namespace CatBLST.AVX.MatrixNoK

def gtCenteredRemainder (value : Nat) : Int :=
  let remainder := value % gtRadix
  if remainder > gtRadix / 2 then
    (remainder : Int) - (gtRadix : Int)
  else
    (remainder : Int)

def gtCenteredQuotient (value : Nat) : Nat :=
  value / gtRadix + if value % gtRadix > gtRadix / 2 then 1 else 0

theorem gt_centered_division_identity (value : Nat) :
    (value : Int) = gtCenteredRemainder value +
      (gtRadix : Int) * (gtCenteredQuotient value : Int) := by
  have division := Nat.mod_add_div value gtRadix
  have castDivision := congrArg (fun x : Nat => (x : Int)) division
  push_cast at castDivision
  by_cases carry : value % gtRadix > gtRadix / 2
  · simp only [gtCenteredRemainder, gtCenteredQuotient, carry, if_true]
    rw [← castDivision]
    push_cast
    ring
  · simp only [gtCenteredRemainder, gtCenteredQuotient, carry, if_false,
      add_zero]
    exact castDivision.symm

def gtCenterOffset : Nat := gtRadix / 2 - 1

theorem gt_centered_quotient_closed (value : Nat) :
    gtCenteredQuotient value = (value + gtCenterOffset) / gtRadix := by
  let remainder := value % gtRadix
  let quotient := value / gtRadix
  have division : remainder + gtRadix * quotient = value :=
    Nat.mod_add_div value gtRadix
  have rearranged :
      value + gtCenterOffset =
        (remainder + gtCenterOffset) + gtRadix * quotient := by
    omega
  simp only [gtCenteredQuotient]
  change quotient + (if remainder > gtRadix / 2 then 1 else 0) = _
  rw [rearranged, Nat.add_mul_div_left _ _ gt_radix_positive]
  have remainderBound : remainder < gtRadix :=
    Nat.mod_lt value gt_radix_positive
  have radixTwiceHalf : gtRadix = 2 * (gtRadix / 2) := by
    native_decide
  have centerOffsetRelation : gtCenterOffset + 1 = gtRadix / 2 := by
    native_decide
  by_cases carry : remainder > gtRadix / 2
  · simp only [carry, if_true]
    have low : 1 * gtRadix <= remainder + gtCenterOffset := by
      omega
    have high : remainder + gtCenterOffset < (1 + 1) * gtRadix := by
      omega
    rw [Nat.div_eq_of_lt_le low high]
    omega
  · simp only [carry, if_false, add_zero]
    have low : remainder + gtCenterOffset < gtRadix := by
      omega
    rw [Nat.div_eq_of_lt low]
    omega

theorem gt_centered_quotient_mono {left right : Nat} (h : left <= right) :
    gtCenteredQuotient left <= gtCenteredQuotient right := by
  rw [gt_centered_quotient_closed, gt_centered_quotient_closed]
  exact Nat.div_le_div_right (Nat.add_le_add_right h gtCenterOffset)

theorem gt_centered_remainder_bounds (value : Nat) :
    -((gtRadix / 2 : Nat) : Int) <= gtCenteredRemainder value ∧
      gtCenteredRemainder value <= ((gtRadix / 2 : Nat) : Int) := by
  let remainder := value % gtRadix
  have remainderBound : remainder < gtRadix :=
    Nat.mod_lt value gt_radix_positive
  have radixTwiceHalf : gtRadix = 2 * (gtRadix / 2) := by
    native_decide
  have remainderBoundInt : (remainder : Int) < (gtRadix : Int) := by
    exact_mod_cast remainderBound
  have radixTwiceHalfInt :
      (gtRadix : Int) = 2 * ((gtRadix / 2 : Nat) : Int) := by
    exact_mod_cast radixTwiceHalf
  by_cases carry : remainder > gtRadix / 2
  · have carryInt : ((gtRadix / 2 : Nat) : Int) < (remainder : Int) := by
      exact_mod_cast carry
    simp only [gtCenteredRemainder, remainder, carry, if_true]
    omega
  · have noCarry : remainder <= gtRadix / 2 := by omega
    have noCarryInt : (remainder : Int) <= ((gtRadix / 2 : Nat) : Int) := by
      exact_mod_cast noCarry
    have remainderNonnegative : (0 : Int) <= (remainder : Int) := by omega
    simp only [gtCenteredRemainder, remainder, carry, if_false]
    omega

def gtBalancedQ1 (scalar : Nat) : Nat := gtCenteredQuotient scalar
def gtBalancedQ2 (scalar : Nat) : Nat := gtCenteredQuotient (gtBalancedQ1 scalar)
def gtBalancedQ3 (scalar : Nat) : Nat := gtCenteredQuotient (gtBalancedQ2 scalar)
def gtBalancedTail (scalar : Nat) : Nat := gtCenteredQuotient (gtBalancedQ3 scalar)

theorem gt_balanced_tail_bound (scalar : Nat) (scalarBound : scalar < 2 ^ 256) :
    gtBalancedTail scalar <= 2 := by
  let maximum := 2 ^ 256 - 1
  have inputBound : scalar <= maximum := by
    omega
  have h1 := gt_centered_quotient_mono inputBound
  have h2 := gt_centered_quotient_mono h1
  have h3 := gt_centered_quotient_mono h2
  have h4 := gt_centered_quotient_mono h3
  have maximumTail : gtBalancedTail maximum = 2 := by
    native_decide
  change gtBalancedTail scalar <= 2
  rw [← maximumTail]
  simpa [gtBalancedTail, gtBalancedQ3, gtBalancedQ2, gtBalancedQ1] using h4

/-- Four signed digits before the source folds the remaining quotient. -/
def gtBalancedRawDigits (scalar : Nat) : Fin 4 -> Int :=
  ![gtCenteredRemainder scalar,
    gtCenteredRemainder (gtBalancedQ1 scalar),
    gtCenteredRemainder (gtBalancedQ2 scalar),
    gtCenteredRemainder (gtBalancedQ3 scalar)]

theorem gt_balanced_raw_split_recompose (scalar : Nat) :
    (scalar : Int) =
      gtBalancedRawDigits scalar 0 +
      (gtRadix : Int) * gtBalancedRawDigits scalar 1 +
      (gtRadix : Int) ^ 2 * gtBalancedRawDigits scalar 2 +
      (gtRadix : Int) ^ 3 * gtBalancedRawDigits scalar 3 +
      (gtRadix : Int) ^ 4 * (gtBalancedTail scalar : Int) := by
  have h0 := gt_centered_division_identity scalar
  have h1 := gt_centered_division_identity (gtBalancedQ1 scalar)
  have h2 := gt_centered_division_identity (gtBalancedQ2 scalar)
  have h3 := gt_centered_division_identity (gtBalancedQ3 scalar)
  change (scalar : Int) = gtCenteredRemainder scalar +
    (gtRadix : Int) * (gtBalancedQ1 scalar : Int) at h0
  change (gtBalancedQ1 scalar : Int) =
    gtCenteredRemainder (gtBalancedQ1 scalar) +
    (gtRadix : Int) * (gtBalancedQ2 scalar : Int) at h1
  change (gtBalancedQ2 scalar : Int) =
    gtCenteredRemainder (gtBalancedQ2 scalar) +
    (gtRadix : Int) * (gtBalancedQ3 scalar : Int) at h2
  change (gtBalancedQ3 scalar : Int) =
    gtCenteredRemainder (gtBalancedQ3 scalar) +
    (gtRadix : Int) * (gtBalancedTail scalar : Int) at h3
  simp only [gtBalancedRawDigits]
  change (scalar : Int) =
    gtCenteredRemainder scalar +
    (gtRadix : Int) * gtCenteredRemainder (gtBalancedQ1 scalar) +
    (gtRadix : Int) ^ 2 * gtCenteredRemainder (gtBalancedQ2 scalar) +
    (gtRadix : Int) ^ 3 * gtCenteredRemainder (gtBalancedQ3 scalar) +
    (gtRadix : Int) ^ 4 * (gtBalancedTail scalar : Int)
  rw [h0, h1, h2, h3]
  ring

/-- Exact final carry fold and odd-coordinate sign changes from the source. -/
def gtBalancedSourceDigits (scalar : Nat) : Fin 4 -> Int :=
  let raw := gtBalancedRawDigits scalar
  let tail := (gtBalancedTail scalar : Int)
  ![raw 0 - tail, -raw 1, raw 2 + tail, -raw 3]

theorem gt_balanced_source_digits_fit_int64
    (scalar : Nat) (scalarBound : scalar < 2 ^ 256) :
    forall index : Fin 4,
      -((2 ^ 63 : Nat) : Int) < gtBalancedSourceDigits scalar index ∧
      gtBalancedSourceDigits scalar index < ((2 ^ 63 : Nat) : Int) := by
  have raw0 := gt_centered_remainder_bounds scalar
  have raw1 := gt_centered_remainder_bounds (gtBalancedQ1 scalar)
  have raw2 := gt_centered_remainder_bounds (gtBalancedQ2 scalar)
  have raw3 := gt_centered_remainder_bounds (gtBalancedQ3 scalar)
  have tail := gt_balanced_tail_bound scalar scalarBound
  have tailInt : (gtBalancedTail scalar : Int) <= 2 := by
    exact_mod_cast tail
  have halfFits :
      ((gtRadix / 2 : Nat) : Int) + 2 < ((2 ^ 63 : Nat) : Int) := by
    native_decide
  intro index
  fin_cases index <;>
    simp [gtBalancedSourceDigits, gtBalancedRawDigits] <;>
    omega

theorem gt_scalar_order_radix_polynomial :
    blsScalarOrder = gtRadix ^ 4 - gtRadix ^ 2 + 1 := by
  native_decide

theorem gt_radix_fourth_eq_radix_squared_sub_one :
    (gtRadix : GtExponent) ^ 4 = (gtRadix : GtExponent) ^ 2 - 1 := by
  native_decide

def gtBalancedSplitExponent (scalar : Nat) : GtExponent :=
  Finset.univ.sum (fun index : Fin 4 =>
    (gtBalancedSourceDigits scalar index : GtExponent) *
      gtFrobeniusBasisExponent index)

/-- The balanced source coefficients reconstruct every 256-bit scalar mod r. -/
theorem gt_balanced_split_exponent_correct (scalar : Nat) :
    gtBalancedSplitExponent scalar = (scalar : GtExponent) := by
  have recomposed := gt_balanced_raw_split_recompose scalar
  have castRecomposed := congrArg
    (fun value : Int => (value : GtExponent)) recomposed
  push_cast at castRecomposed
  simp only [gtBalancedSplitExponent, Fin.sum_univ_four]
  simp [gtBalancedSourceDigits, gtFrobeniusBasisExponent]
  rw [gt_frobenius_eigenvalue_is_negative_z]
  rw [gt_radix_fourth_eq_radix_squared_sub_one] at castRecomposed
  calc
    _ =
        (gtBalancedRawDigits scalar 0 : GtExponent) +
        (gtRadix : GtExponent) * (gtBalancedRawDigits scalar 1 : GtExponent) +
        (gtRadix : GtExponent) ^ 2 * (gtBalancedRawDigits scalar 2 : GtExponent) +
        (gtRadix : GtExponent) ^ 3 * (gtBalancedRawDigits scalar 3 : GtExponent) +
        ((gtRadix : GtExponent) ^ 2 - 1) *
          (gtBalancedTail scalar : GtExponent) := by ring
    _ = (scalar : GtExponent) := castRecomposed.symm

end CatBLST.AVX.MatrixNoK
