/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.CurveScalarSplit
import Mathlib.Tactic

set_option maxHeartbeats 800000

/-!
# Fixed-window Booth executor used by G1 GLV and G2 GLS

This file models the selected `table_log_size = 4` production path in
`PrecomputedTable::build`, `booth_encode`, and `table_msm`.  Digits are stored
low-to-high, while `boothExecute` is Horner evaluation from the recursive tail
and therefore matches the source's high-to-low add-and-five-doubles loop.
-/

namespace CatBLST.AVX.MatrixNoK

def boothHalf (width : Nat) : Nat := 2 ^ width

def boothBase (width : Nat) : Nat := 2 ^ (width + 1)

theorem booth_base_eq_two_half (width : Nat) :
    boothBase width = 2 * boothHalf width := by
  simp [boothBase, boothHalf, pow_succ, Nat.mul_comm]

theorem booth_half_positive (width : Nat) : 0 < boothHalf width := by
  exact pow_pos (by norm_num) width

theorem booth_base_positive (width : Nat) : 0 < boothBase width := by
  exact pow_pos (by norm_num) (width + 1)

/-- Signed value represented by the source's `(booth_index, booth_sign)`. -/
def boothDigit (width value : Nat) : Int :=
  let raw := value % boothBase width
  if raw < boothHalf width then
    (raw : Int)
  else
    (raw : Int) - (boothBase width : Int)

/-- Carry/quotient state for the next source window. -/
def boothNext (width value : Nat) : Nat :=
  let raw := value % boothBase width
  value / boothBase width + if raw < boothHalf width then 0 else 1

theorem booth_step_recompose (width value : Nat) :
    (value : Int) = boothDigit width value +
      (boothBase width : Int) * boothNext width value := by
  have split := Nat.mod_add_div value (boothBase width)
  have castSplit := congrArg (fun x : Nat => (x : Int)) split
  push_cast at castSplit
  by_cases low : value % boothBase width < boothHalf width
  · simp [boothDigit, boothNext, low]
    linarith
  · simp [boothDigit, boothNext, low]
    linarith

theorem booth_next_eq_rounded_div (width value : Nat) :
    boothNext width value =
      (value + boothHalf width) / boothBase width := by
  let base := boothBase width
  let half := boothHalf width
  let raw := value % base
  let quotient := value / base
  have split : raw + base * quotient = value := Nat.mod_add_div value base
  have rawBound : raw < base := Nat.mod_lt _ (booth_base_positive width)
  have baseShape : base = 2 * half := booth_base_eq_two_half width
  by_cases low : raw < half
  · have sumBound : raw + half < base := by omega
    have decompose : value + half = (raw + half) + base * quotient := by
      omega
    rw [boothNext]
    simp only [raw, base, half] at low ⊢
    simp only [low, if_true, add_zero]
    rw [decompose, Nat.add_mul_div_left _ _ (booth_base_positive width)]
    change quotient = (raw + half) / base + quotient
    rw [Nat.div_eq_of_lt sumBound]
    simp
  · have halfLe : half ≤ raw := Nat.le_of_not_gt low
    have sumLower : base ≤ raw + half := by omega
    have sumUpper : raw + half < 2 * base := by omega
    have remainderBound : raw + half - base < base := by omega
    have firstDecompose : value + half = raw + half + base * quotient := by
      omega
    have decompose : value + half =
        (raw + half - base) + base * (quotient + 1) := by
      rw [firstDecompose, Nat.mul_add]
      omega
    rw [boothNext]
    simp only [raw, base, half] at low ⊢
    simp only [low, if_false]
    rw [decompose, Nat.add_mul_div_left _ _ (booth_base_positive width)]
    rw [Nat.div_eq_of_lt remainderBound]
    simp only [zero_add]
    rfl

/-- Low-to-high signed Booth digits for a fixed number of source windows. -/
def boothRecode (width : Nat) : Nat -> Nat -> List Int
  | 0, _ => []
  | windows + 1, value =>
      boothDigit width value ::
        boothRecode width windows (boothNext width value)

def boothAfter (width : Nat) : Nat -> Nat -> Nat
  | 0, value => value
  | windows + 1, value =>
      boothAfter width windows (boothNext width value)

/-- High-to-low Horner execution of the stored low-to-high digit list. -/
def boothExecute (width : Nat) : List Int -> Int
  | [] => 0
  | digit :: digits =>
      digit + (boothBase width : Int) * boothExecute width digits

theorem booth_recode_recompose (width windows value : Nat) :
    (value : Int) = boothExecute width (boothRecode width windows value) +
      (boothBase width : Int) ^ windows *
        (boothAfter width windows value : Int) := by
  induction windows generalizing value with
  | zero => simp [boothRecode, boothExecute, boothAfter]
  | succ windows ih =>
      simp only [boothRecode, boothExecute, boothAfter]
      have step := booth_step_recompose width value
      have tail := ih (boothNext width value)
      calc
        (value : Int) = boothDigit width value +
            (boothBase width : Int) * (boothNext width value : Int) := step
        _ = boothDigit width value + (boothBase width : Int) *
            (boothExecute width
                (boothRecode width windows (boothNext width value)) +
              (boothBase width : Int) ^ windows *
                (boothAfter width windows (boothNext width value) : Int)) := by
                  rw [← tail]
        _ = boothDigit width value + (boothBase width : Int) *
              boothExecute width
                (boothRecode width windows (boothNext width value)) +
            (boothBase width : Int) ^ (windows + 1) *
              (boothAfter width windows (boothNext width value) : Int) := by
                rw [pow_succ]
                ring

/-- Accumulated rounding bias after `windows` centered-radix steps. -/
def boothBias (width : Nat) : Nat -> Nat
  | 0 => 0
  | windows + 1 =>
      boothHalf width + boothBase width * boothBias width windows

theorem booth_bias_succ_alt (width windows : Nat) :
    boothBias width (windows + 1) =
      boothHalf width + boothBase width * boothBias width windows := by
  rfl

theorem booth_after_closed_form (width windows value : Nat) :
    boothAfter width windows value =
      (value + boothBias width windows) / boothBase width ^ windows := by
  induction windows generalizing value with
  | zero => simp [boothAfter, boothBias]
  | succ windows ih =>
      rw [boothAfter, ih, booth_next_eq_rounded_div]
      rw [booth_bias_succ_alt, pow_succ]
      calc
        ((value + boothHalf width) / boothBase width +
              boothBias width windows) / boothBase width ^ windows =
            ((value + boothHalf width +
                boothBase width * boothBias width windows) /
              boothBase width) / boothBase width ^ windows := by
                rw [Nat.add_mul_div_left _ _ (booth_base_positive width)]
        _ = (value + boothHalf width +
              boothBase width * boothBias width windows) /
              (boothBase width * boothBase width ^ windows) := by
                rw [Nat.div_div_eq_div_mul]
        _ = (value + (boothHalf width +
              boothBase width * boothBias width windows)) /
              (boothBase width ^ windows * boothBase width) := by
                congr 1
                · omega
                · ring

theorem booth_after_eq_zero_of_capacity
    (width windows value bound : Nat)
    (valueBound : value < bound)
    (capacity : bound - 1 + boothBias width windows <
      boothBase width ^ windows) :
    boothAfter width windows value = 0 := by
  rw [booth_after_closed_form]
  apply Nat.div_eq_of_lt
  omega

theorem g1_booth_capacity :
    curveRadixSquared - 1 + boothBias 4 26 < boothBase 4 ^ 26 := by
  native_decide

theorem g2_booth_capacity :
    curveRadix - 1 + boothBias 4 13 < boothBase 4 ^ 13 := by
  native_decide

theorem g1_booth_after_zero (value : Nat)
    (bound : value < curveRadixSquared) :
    boothAfter 4 26 value = 0 :=
  booth_after_eq_zero_of_capacity 4 26 value curveRadixSquared
    bound g1_booth_capacity

theorem g2_booth_after_zero (value : Nat)
    (bound : value < curveRadix) :
    boothAfter 4 13 value = 0 :=
  booth_after_eq_zero_of_capacity 4 13 value curveRadix
    bound g2_booth_capacity

theorem g1_booth_execute_scalar (value : Nat)
    (bound : value < curveRadixSquared) :
    boothExecute 4 (boothRecode 4 26 value) = (value : Int) := by
  have recomposed := booth_recode_recompose 4 26 value
  rw [g1_booth_after_zero value bound] at recomposed
  simpa using recomposed.symm

theorem g2_booth_execute_scalar (value : Nat)
    (bound : value < curveRadix) :
    boothExecute 4 (boothRecode 4 13 value) = (value : Int) := by
  have recomposed := booth_recode_recompose 4 13 value
  rw [g2_booth_after_zero value bound] at recomposed
  simpa using recomposed.symm

/-- The source table loop constructs row `i` as `(i + 1)P`. -/
def precomputedTableExponent (base : CurveExponent) (index : Nat) : CurveExponent :=
  (index + 1 : Nat) * base

theorem precomputed_table_odd_step (base : CurveExponent) (index : Nat)
    (odd : index % 2 = 1) :
    precomputedTableExponent base index =
      2 * precomputedTableExponent base (index / 2) := by
  simp only [precomputedTableExponent]
  have split := Nat.mod_add_div index 2
  rw [odd] at split
  have indexShape : index + 1 = 2 * (index / 2 + 1) := by omega
  rw [indexShape]
  push_cast
  ring

theorem precomputed_table_even_step (base : CurveExponent) (index : Nat)
    (positive : 0 < index) (even : index % 2 = 0) :
    precomputedTableExponent base index =
      precomputedTableExponent base (index - 1) + base := by
  simp only [precomputedTableExponent]
  have indexShape : index + 1 = (index - 1 + 1) + 1 := by omega
  rw [indexShape]
  push_cast
  ring

def boothComponentExponent
    (windows _bound : Nat) (basis : CurveExponent) (scalar : Nat) : CurveExponent :=
  (boothExecute 4 (boothRecode 4 windows scalar) : CurveExponent) * basis

theorem g1_booth_executor_correct (scalar : Nat)
    (canonical : scalar < blsScalarOrder) :
    Finset.univ.sum (fun index : Fin 2 =>
      boothComponentExponent 26 curveRadixSquared
        (g1GlvBasisExponent index) (g1GlvSplit scalar index)) =
      (scalar : CurveExponent) := by
  have bounds := g1_canonical_split_bounds scalar canonical
  simp only [boothComponentExponent, Fin.sum_univ_two]
  rw [g1_booth_execute_scalar _ (bounds 0),
    g1_booth_execute_scalar _ (bounds 1)]
  simpa [g1GlvSplitExponent, Fin.sum_univ_two] using
    (g1_glv_split_exponent_correct scalar)

theorem g2_booth_executor_correct (scalar : Nat)
    (canonical : scalar < blsScalarOrder) :
    Finset.univ.sum (fun index : Fin 4 =>
      boothComponentExponent 13 curveRadix
        (g2GlsBasisExponent index) (g2GlsSplit scalar index)) =
      (scalar : CurveExponent) := by
  have bounds := g2_canonical_split_bounds scalar canonical
  simp only [boothComponentExponent, Fin.sum_univ_four]
  rw [g2_booth_execute_scalar _ (bounds 0),
    g2_booth_execute_scalar _ (bounds 1),
    g2_booth_execute_scalar _ (bounds 2),
    g2_booth_execute_scalar _ (bounds 3)]
  simpa [g2GlsSplitExponent, Fin.sum_univ_four] using
    (g2_gls_split_exponent_correct scalar)

end CatBLST.AVX.MatrixNoK
