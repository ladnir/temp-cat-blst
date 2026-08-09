/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.CurveScalarSplit
import Mathlib.Tactic

set_option maxHeartbeats 800000

/-!
# Variable-time wNAF executors used by G1 GLV and G2 GLS

This file models `recode_wnaf_vartime` and the selected AVX executors.  The
production choices are width 5 with 129 digit slots for each 128-bit G1
component and width 4 with 65 slots for each 64-bit G2 component.

Digits are stored low-to-high.  `wnafExecute` evaluates that storage from the
recursive tail, matching the production high-to-low loop.  The separate
`wnafAfter` term makes exhaustion explicit: no theorem silently assumes that
the fixed source buffer was long enough.
-/

namespace CatBLST.AVX.MatrixNoK

def wnafBase (width : Nat) : Nat := 2 ^ width

def wnafHalf (width : Nat) : Nat := 2 ^ (width - 1)

theorem wnaf_base_positive (width : Nat) : 0 < wnafBase width := by
  exact pow_pos (by norm_num) width

theorem wnaf_base_eq_two_half (width : Nat) (widthPositive : 0 < width) :
    wnafBase width = 2 * wnafHalf width := by
  obtain ⟨width, rfl⟩ := Nat.exists_eq_succ_of_ne_zero (by omega : width ≠ 0)
  simp [wnafBase, wnafHalf, pow_succ, Nat.mul_comm]

def wnafDigit (width value : Nat) : Int :=
  if value % 2 = 0 then
    0
  else
    let raw := value % wnafBase width
    if raw < wnafHalf width then
      (raw : Int)
    else
      (raw : Int) - (wnafBase width : Int)

/-!
`wnafNext` is the arithmetic effect of the source's add/subtract followed by
the multi-limb right shift.  The quotient-shaped definition avoids importing
machine-limb details; `wnaf_step_recompose` proves the exact integer relation.
-/
def wnafNext (width value : Nat) : Nat :=
  if value % 2 = 0 then
    value / 2
  else
    let raw := value % wnafBase width
    let quotient := value / wnafBase width
    if raw < wnafHalf width then
      wnafHalf width * quotient
    else
      wnafHalf width * (quotient + 1)

theorem wnaf_step_recompose (width value : Nat) (widthPositive : 0 < width) :
    (value : Int) = wnafDigit width value + 2 * wnafNext width value := by
  have splitTwo := Nat.mod_add_div value 2
  have splitBase := Nat.mod_add_div value (wnafBase width)
  have baseShape := wnaf_base_eq_two_half width widthPositive
  by_cases even : value % 2 = 0
  · simp only [wnafDigit, wnafNext, even, if_true]
    have natural : value = 2 * (value / 2) := by omega
    have castNatural : (value : Int) = 2 * (value / 2 : Nat) := by
      exact_mod_cast natural
    simpa using castNatural
  · let raw := value % wnafBase width
    let quotient := value / wnafBase width
    have sourceSplit : raw + wnafBase width * quotient = value := splitBase
    have castSplit : (raw : Int) + (wnafBase width : Int) *
        (quotient : Int) = (value : Int) := by
      exact_mod_cast sourceSplit
    have castBase : (wnafBase width : Int) =
        2 * (wnafHalf width : Int) := by
      exact_mod_cast baseShape
    by_cases low : raw < wnafHalf width
    · simp only [wnafDigit, wnafNext, even, if_false]
      simp only [raw] at low ⊢
      simp only [low, if_true]
      change (value : Int) = (raw : Int) +
        2 * (wnafHalf width * quotient : Nat)
      push_cast
      rw [castBase] at castSplit
      linarith
    · simp only [wnafDigit, wnafNext, even, if_false]
      simp only [raw] at low ⊢
      simp only [low, if_false]
      change (value : Int) = (raw : Int) - (wnafBase width : Int) +
        2 * (wnafHalf width * (quotient + 1) : Nat)
      push_cast
      rw [castBase] at castSplit ⊢
      linarith

def wnafRecode (width : Nat) : Nat → Nat → List Int
  | 0, _ => []
  | slots + 1, value =>
      wnafDigit width value :: wnafRecode width slots (wnafNext width value)

def wnafAfter (width : Nat) : Nat → Nat → Nat
  | 0, value => value
  | slots + 1, value => wnafAfter width slots (wnafNext width value)

def wnafExecute : List Int → Int
  | [] => 0
  | digit :: digits => digit + 2 * wnafExecute digits

theorem wnaf_recode_recompose
    (width slots value : Nat) (widthPositive : 0 < width) :
    (value : Int) = wnafExecute (wnafRecode width slots value) +
      2 ^ slots * (wnafAfter width slots value : Int) := by
  induction slots generalizing value with
  | zero => simp [wnafRecode, wnafExecute, wnafAfter]
  | succ slots ih =>
      simp only [wnafRecode, wnafExecute, wnafAfter]
      have step := wnaf_step_recompose width value widthPositive
      have tail := ih (wnafNext width value)
      calc
        (value : Int) = wnafDigit width value +
            2 * (wnafNext width value : Int) := step
        _ = wnafDigit width value + 2 *
            (wnafExecute (wnafRecode width slots (wnafNext width value)) +
              2 ^ slots *
                (wnafAfter width slots (wnafNext width value) : Int)) := by
                  rw [← tail]
        _ = wnafDigit width value +
              2 * wnafExecute
                (wnafRecode width slots (wnafNext width value)) +
            2 ^ (slots + 1) *
              (wnafAfter width slots (wnafNext width value) : Int) := by
                rw [pow_succ]
                ring

theorem wnaf_execute_of_exhausted
    (width slots value : Nat) (widthPositive : 0 < width)
    (exhausted : wnafAfter width slots value = 0) :
    wnafExecute (wnafRecode width slots value) = (value : Int) := by
  have recomposed := wnaf_recode_recompose width slots value widthPositive
  rw [exhausted] at recomposed
  simpa using recomposed.symm

theorem wnaf_after_add (width first second value : Nat) :
    wnafAfter width (first + second) value =
      wnafAfter width second (wnafAfter width first value) := by
  induction first generalizing value with
  | zero => simp [wnafAfter]
  | succ first ih =>
      rw [Nat.succ_add]
      simp only [wnafAfter]
      exact ih (wnafNext width value)

/-!
Backward minimum-predecessor bound.  A positive residual can come from an even
predecessor `2m`, or from a centered odd predecessor.  The `+ 1` is essential:
the centered half-radix is even for the selected widths, while the source only
takes this branch for an odd scalar.
-/
def wnafRequiredStep (width residual : Nat) : Nat :=
  min (2 * residual)
    (wnafBase width *
        ((residual + wnafHalf width - 1) / wnafHalf width - 1) +
      wnafHalf width + 1)

def wnafRequired (width : Nat) : Nat → Nat
  | 0 => 1
  | slots + 1 => wnafRequiredStep width (wnafRequired width slots)

theorem wnaf_required_step_five (residual value : Nat)
    (nextBound : residual ≤ wnafNext 5 value) :
    wnafRequiredStep 5 residual ≤ value := by
  by_cases even : value % 2 = 0
  · have candidate : 2 * residual ≤ value := by
      simp only [wnafNext, even, if_true] at nextBound
      omega
    exact (Nat.min_le_left _ _).trans candidate
  · let raw := value % 32
    let quotient := value / 32
    have split : raw + 32 * quotient = value := Nat.mod_add_div value 32
    have rawBound : raw < 32 := Nat.mod_lt _ (by norm_num)
    have oddRemainder : raw % 2 = 1 := by
      have valueParity : value % 2 = 1 := by omega
      omega
    simp only [wnafNext, even, if_false] at nextBound
    norm_num [wnafBase, wnafHalf] at nextBound
    change residual ≤
      (if raw < 16 then 16 * quotient else 16 * (quotient + 1)) at nextBound
    by_cases low : raw < 16
    · have quotientBound : residual ≤ 16 * quotient := by
        simpa only [low, if_true] using nextBound
      by_cases residualZero : residual = 0
      · subst residual
        simp [wnafRequiredStep]
      have candidate :
          32 * ((residual + 15) / 16 - 1) + 17 ≤ value := by
        omega
      unfold wnafRequiredStep
      apply (Nat.min_le_right _ _).trans
      norm_num [wnafBase, wnafHalf]
      omega
    · have quotientBound : residual ≤ 16 * (quotient + 1) := by
        simpa only [low, if_false] using nextBound
      by_cases residualZero : residual = 0
      · subst residual
        simp [wnafRequiredStep]
      have candidate :
          32 * ((residual + 15) / 16 - 1) + 17 ≤ value := by
        omega
      unfold wnafRequiredStep
      apply (Nat.min_le_right _ _).trans
      norm_num [wnafBase, wnafHalf]
      omega

theorem wnaf_required_step_four (residual value : Nat)
    (nextBound : residual ≤ wnafNext 4 value) :
    wnafRequiredStep 4 residual ≤ value := by
  by_cases even : value % 2 = 0
  · have candidate : 2 * residual ≤ value := by
      simp only [wnafNext, even, if_true] at nextBound
      omega
    exact (Nat.min_le_left _ _).trans candidate
  · let raw := value % 16
    let quotient := value / 16
    have split : raw + 16 * quotient = value := Nat.mod_add_div value 16
    have rawBound : raw < 16 := Nat.mod_lt _ (by norm_num)
    have oddRemainder : raw % 2 = 1 := by
      have valueParity : value % 2 = 1 := by omega
      omega
    simp only [wnafNext, even, if_false] at nextBound
    norm_num [wnafBase, wnafHalf] at nextBound
    change residual ≤
      (if raw < 8 then 8 * quotient else 8 * (quotient + 1)) at nextBound
    by_cases low : raw < 8
    · have quotientBound : residual ≤ 8 * quotient := by
        simpa only [low, if_true] using nextBound
      by_cases residualZero : residual = 0
      · subst residual
        simp [wnafRequiredStep]
      have candidate :
          16 * ((residual + 7) / 8 - 1) + 9 ≤ value := by
        omega
      unfold wnafRequiredStep
      apply (Nat.min_le_right _ _).trans
      norm_num [wnafBase, wnafHalf]
      omega
    · have quotientBound : residual ≤ 8 * (quotient + 1) := by
        simpa only [low, if_false] using nextBound
      by_cases residualZero : residual = 0
      · subst residual
        simp [wnafRequiredStep]
      have candidate :
          16 * ((residual + 7) / 8 - 1) + 9 ≤ value := by
        omega
      unfold wnafRequiredStep
      apply (Nat.min_le_right _ _).trans
      norm_num [wnafBase, wnafHalf]
      omega

theorem wnaf_required_five_lower (slots value : Nat)
    (positive : 0 < wnafAfter 5 slots value) :
    wnafRequired 5 slots ≤ value := by
  induction slots generalizing value with
  | zero =>
      change 1 ≤ value
      change 0 < value at positive
      omega
  | succ slots ih =>
      simp only [wnafAfter] at positive
      have nextBound := ih (wnafNext 5 value) positive
      simpa only [wnafRequired] using
        wnaf_required_step_five (wnafRequired 5 slots) value nextBound

theorem wnaf_required_four_lower (slots value : Nat)
    (positive : 0 < wnafAfter 4 slots value) :
    wnafRequired 4 slots ≤ value := by
  induction slots generalizing value with
  | zero =>
      change 1 ≤ value
      change 0 < value at positive
      omega
  | succ slots ih =>
      simp only [wnafAfter] at positive
      have nextBound := ih (wnafNext 4 value) positive
      simpa only [wnafRequired] using
        wnaf_required_step_four (wnafRequired 4 slots) value nextBound

theorem g1_wnaf_required_exceeds_component :
    curveRadixSquared < wnafRequired 5 129 := by
  native_decide

theorem g2_wnaf_required_exceeds_component :
    curveRadix < wnafRequired 4 65 := by
  native_decide

theorem g1_wnaf_exhausted (value : Nat)
    (bound : value < curveRadixSquared) :
    wnafAfter 5 129 value = 0 := by
  by_contra nonzero
  have positive : 0 < wnafAfter 5 129 value := Nat.pos_of_ne_zero nonzero
  have required := wnaf_required_five_lower 129 value positive
  have exceeds := g1_wnaf_required_exceeds_component
  omega

theorem g2_wnaf_exhausted (value : Nat)
    (bound : value < curveRadix) :
    wnafAfter 4 65 value = 0 := by
  by_contra nonzero
  have positive : 0 < wnafAfter 4 65 value := Nat.pos_of_ne_zero nonzero
  have required := wnaf_required_four_lower 65 value positive
  have exceeds := g2_wnaf_required_exceeds_component
  omega

/-! The production odd table contains `P, 3P, 5P, ...`. -/
def wnafOddTableExponent
    (base : CurveExponent) (index : Nat) : CurveExponent :=
  (2 * index + 1 : Nat) * base

theorem wnaf_odd_table_zero (base : CurveExponent) :
    wnafOddTableExponent base 0 = base := by
  simp [wnafOddTableExponent]

theorem wnaf_odd_table_step (base : CurveExponent) (index : Nat) :
    wnafOddTableExponent base (index + 1) =
      wnafOddTableExponent base index + 2 * base := by
  simp only [wnafOddTableExponent]
  push_cast
  ring

def wnafComponentExponent
    (width slots : Nat) (basis : CurveExponent) (scalar : Nat) : CurveExponent :=
  (wnafExecute (wnafRecode width slots scalar) : CurveExponent) * basis

theorem g1_wnaf_executor_correct (scalar : Nat)
    (canonical : scalar < blsScalarOrder) :
    Finset.univ.sum (fun index : Fin 2 =>
      wnafComponentExponent 5 129
        (g1GlvBasisExponent index) (g1GlvSplit scalar index)) =
      (scalar : CurveExponent) := by
  have bounds := g1_canonical_split_bounds scalar canonical
  simp only [wnafComponentExponent, Fin.sum_univ_two]
  rw [wnaf_execute_of_exhausted 5 129 _ (by norm_num)
        (g1_wnaf_exhausted _ (bounds 0)),
    wnaf_execute_of_exhausted 5 129 _ (by norm_num)
        (g1_wnaf_exhausted _ (bounds 1))]
  simpa [g1GlvSplitExponent, Fin.sum_univ_two] using
    (g1_glv_split_exponent_correct scalar)

theorem g2_wnaf_executor_correct (scalar : Nat)
    (canonical : scalar < blsScalarOrder) :
    Finset.univ.sum (fun index : Fin 4 =>
      wnafComponentExponent 4 65
        (g2GlsBasisExponent index) (g2GlsSplit scalar index)) =
      (scalar : CurveExponent) := by
  have bounds := g2_canonical_split_bounds scalar canonical
  simp only [wnafComponentExponent, Fin.sum_univ_four]
  rw [wnaf_execute_of_exhausted 4 65 _ (by norm_num)
        (g2_wnaf_exhausted _ (bounds 0)),
    wnaf_execute_of_exhausted 4 65 _ (by norm_num)
        (g2_wnaf_exhausted _ (bounds 1)),
    wnaf_execute_of_exhausted 4 65 _ (by norm_num)
        (g2_wnaf_exhausted _ (bounds 2)),
    wnaf_execute_of_exhausted 4 65 _ (by norm_num)
        (g2_wnaf_exhausted _ (bounds 3))]
  simpa [g2GlsSplitExponent, Fin.sum_univ_four] using
    (g2_gls_split_exponent_correct scalar)

end CatBLST.AVX.MatrixNoK
