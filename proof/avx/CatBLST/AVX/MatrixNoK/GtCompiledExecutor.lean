/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.GtCompiledSchedule
import Mathlib.Tactic

set_option maxHeartbeats 800000

/-!
# Compiled GT GLV4 executor

The model below evaluates the exact high/low-nibble schedule.  It proves the
compiler's delayed square runs, square-only chunks, first-table initialization,
and trusted executor all refine the direct 64-bit joint scan.
-/

namespace CatBLST.AVX.MatrixNoK

def gtBalancedSignedBasis (scalar : Nat) (index : Fin 4) : GtExponent :=
  if gtBalancedSourceDigits scalar index < 0 then
    -gtFrobeniusBasisExponent index
  else
    gtFrobeniusBasisExponent index

theorem gt_magnitude_signed_basis
    (scalar : Nat) (index : Fin 4) :
    (gtBalancedMagnitudes scalar index : GtExponent) *
        gtBalancedSignedBasis scalar index =
      (gtBalancedSourceDigits scalar index : GtExponent) *
        gtFrobeniusBasisExponent index := by
  let digit := gtBalancedSourceDigits scalar index
  change (digit.natAbs : GtExponent) *
      (if digit < 0 then -gtFrobeniusBasisExponent index
       else gtFrobeniusBasisExponent index) =
    (digit : GtExponent) * gtFrobeniusBasisExponent index
  by_cases negative : digit < 0
  · have magnitudeInt : (digit.natAbs : Int) = -digit := by
      calc
        (digit.natAbs : Int) = ((-digit).natAbs : Int) := by
          rw [Int.natAbs_neg]
        _ = -digit := Int.natAbs_of_nonneg (by omega)
    have magnitudeGt : (digit.natAbs : GtExponent) = -(digit : GtExponent) := by
      calc
        (digit.natAbs : GtExponent) =
            (((digit.natAbs : Nat) : Int) : GtExponent) := by norm_num
        _ = ((-digit : Int) : GtExponent) := by rw [magnitudeInt]
        _ = -(digit : GtExponent) := by push_cast; rfl
    simp only [negative, if_true]
    rw [magnitudeGt]
    ring
  · have nonnegative : 0 <= digit := by omega
    have magnitudeInt : (digit.natAbs : Int) = digit :=
      Int.natAbs_of_nonneg nonnegative
    have magnitudeGt : (digit.natAbs : GtExponent) = (digit : GtExponent) := by
      calc
        (digit.natAbs : GtExponent) =
            (((digit.natAbs : Nat) : Int) : GtExponent) := by norm_num
        _ = (digit : GtExponent) := by rw [magnitudeInt]
    simp only [negative, if_false]
    rw [magnitudeGt]

theorem gt_balanced_magnitude_bound
    (scalar : Nat) (scalarBound : scalar < 2 ^ 256) (index : Fin 4) :
    gtBalancedMagnitudes scalar index < 2 ^ 63 := by
  let digit := gtBalancedSourceDigits scalar index
  have signedBounds := gt_balanced_source_digits_fit_int64 scalar scalarBound index
  by_cases negative : digit < 0
  · have magnitudeInt : (digit.natAbs : Int) = -digit := by
      calc
        (digit.natAbs : Int) = ((-digit).natAbs : Int) := by
          rw [Int.natAbs_neg]
        _ = -digit := Int.natAbs_of_nonneg (by omega)
    change digit.natAbs < 2 ^ 63
    exact_mod_cast (show (digit.natAbs : Int) < ((2 ^ 63 : Nat) : Int) by
      rw [magnitudeInt]
      omega)
  · have nonnegative : 0 <= digit := by omega
    have magnitudeInt : (digit.natAbs : Int) = digit :=
      Int.natAbs_of_nonneg nonnegative
    change digit.natAbs < 2 ^ 63
    exact_mod_cast (show (digit.natAbs : Int) < ((2 ^ 63 : Nat) : Int) by
      rw [magnitudeInt]
      exact signedBounds.2)

def gtBalancedSignMask (scalar : Nat) : Nat :=
  Finset.univ.sum (fun index : Fin 4 =>
    if gtBalancedSourceDigits scalar index < 0 then 2 ^ index.val else 0)

theorem gt_balanced_sign_mask_lt_sixteen (scalar : Nat) :
    gtBalancedSignMask scalar < 16 := by
  simp only [gtBalancedSignMask, Fin.sum_univ_four]
  by_cases h0 : gtBalancedSourceDigits scalar 0 < 0 <;>
    by_cases h1 : gtBalancedSourceDigits scalar 1 < 0 <;>
    by_cases h2 : gtBalancedSourceDigits scalar 2 < 0 <;>
    by_cases h3 : gtBalancedSourceDigits scalar 3 < 0 <;>
    simp [h0, h1, h2, h3]

theorem gt_balanced_sign_mask_bits (scalar : Nat) :
    gtBalancedSignMask scalar % 2 =
        (if gtBalancedSourceDigits scalar 0 < 0 then 1 else 0) ∧
    gtBalancedSignMask scalar / 2 % 2 =
        (if gtBalancedSourceDigits scalar 1 < 0 then 1 else 0) ∧
    gtBalancedSignMask scalar / 4 % 2 =
        (if gtBalancedSourceDigits scalar 2 < 0 then 1 else 0) ∧
    gtBalancedSignMask scalar / 8 % 2 =
        (if gtBalancedSourceDigits scalar 3 < 0 then 1 else 0) := by
  by_cases h0 : gtBalancedSourceDigits scalar 0 < 0 <;>
  by_cases h1 : gtBalancedSourceDigits scalar 1 < 0 <;>
  by_cases h2 : gtBalancedSourceDigits scalar 2 < 0 <;>
  by_cases h3 : gtBalancedSourceDigits scalar 3 < 0 <;>
  norm_num [gtBalancedSignMask, Fin.sum_univ_four, h0, h1, h2, h3]

theorem gt_balanced_signed_basis_reconstructs (scalar : Nat) :
    Finset.univ.sum (fun index : Fin 4 =>
      (gtBalancedMagnitudes scalar index : GtExponent) *
        gtBalancedSignedBasis scalar index) = (scalar : GtExponent) := by
  rw [Finset.sum_congr rfl (fun index _ => gt_magnitude_signed_basis scalar index)]
  exact gt_balanced_split_exponent_correct scalar

def gtEvalCompiledOp
    (table : Fin 16 -> GtExponent) (op : GtCompiledOp)
    (out : GtExponent) : GtExponent :=
  (2 ^ op.squares : GtExponent) * out + table op.index

def gtEvalCompiledOps
    (table : Fin 16 -> GtExponent) : List GtCompiledOp -> GtExponent -> GtExponent
  | [], out => out
  | op :: rest, out => gtEvalCompiledOps table rest (gtEvalCompiledOp table op out)

theorem gt_eval_compiled_append
    (table : Fin 16 -> GtExponent) (left right : List GtCompiledOp)
    (out : GtExponent) :
    gtEvalCompiledOps table (left ++ right) out =
      gtEvalCompiledOps table right (gtEvalCompiledOps table left out) := by
  induction left generalizing out with
  | nil => rfl
  | cons op rest ih =>
      simp only [List.cons_append, gtEvalCompiledOps]
      exact ih (gtEvalCompiledOp table op out)

theorem gt_eval_square_only_replicate
    (table : Fin 16 -> GtExponent) (tableZero : table 0 = 0)
    (count : Nat) (out : GtExponent) :
    gtEvalCompiledOps table (List.replicate count gtSquareOnlyOp) out =
      (2 ^ (15 * count) : GtExponent) * out := by
  induction count generalizing out with
  | zero => simp [gtEvalCompiledOps]
  | succ count ih =>
      rw [List.replicate_succ, gtEvalCompiledOps]
      rw [ih]
      simp only [gtEvalCompiledOp, gtSquareOnlyOp, tableZero, add_zero]
      rw [show 15 * (count + 1) = 15 * count + 15 by ring, pow_add]
      ring

theorem gt_eval_run_mul_ops
    (table : Fin 16 -> GtExponent) (tableZero : table 0 = 0)
    (run : Nat) (index : Fin 16) (runPositive : 0 < run)
    (out : GtExponent) :
    gtEvalCompiledOps table (gtRunMulOps run index) out =
      (2 ^ run : GtExponent) * out + table index := by
  rw [gtRunMulOps, gt_eval_compiled_append,
    gt_eval_square_only_replicate table tableZero]
  simp only [gtEvalCompiledOps, gtEvalCompiledOp]
  have arithmetic : 15 * ((run - 1) / 15) + ((run - 1) % 15 + 1) = run := by
    have division := Nat.mod_add_div (run - 1) 15
    omega
  rw [← mul_assoc, ← pow_add]
  have exponentEquality :
      (run - 1) % 15 + 1 + 15 * ((run - 1) / 15) = run := by
    omega
  rw [exponentEquality]

theorem gt_eval_trailing_square_ops
    (table : Fin 16 -> GtExponent) (tableZero : table 0 = 0)
    (count : Nat) (out : GtExponent) :
    gtEvalCompiledOps table (gtTrailingSquareOps count) out =
      (2 ^ count : GtExponent) * out := by
  rw [gtTrailingSquareOps, gt_eval_compiled_append,
    gt_eval_square_only_replicate table tableZero]
  by_cases remainderZero : count % 15 = 0
  · simp only [remainderZero, if_true, gtEvalCompiledOps]
    have division := Nat.mod_add_div count 15
    have arithmetic : 15 * (count / 15) = count := by omega
    rw [arithmetic]
  · simp only [remainderZero, if_false, gtEvalCompiledOps, gtEvalCompiledOp,
      tableZero, add_zero]
    have division := Nat.mod_add_div count 15
    rw [← mul_assoc, ← pow_add]
    rw [show count % 15 + 15 * (count / 15) = count by omega]

def gtEvalJointBitList
    (table : Fin 16 -> GtExponent) : List (Fin 16) -> GtExponent -> GtExponent
  | [], out => out
  | index :: rest, out =>
      gtEvalJointBitList table rest (out + out + table index)

theorem gt_eval_compile_tail
    (table : Fin 16 -> GtExponent) (tableZero : table 0 = 0)
    (bits : List (Fin 16)) (run : Nat) (runPositive : 0 < run)
    (out : GtExponent) :
    gtEvalCompiledOps table (gtCompileTail bits run) out =
      gtEvalJointBitList table bits ((2 ^ (run - 1) : GtExponent) * out) := by
  induction bits generalizing run out with
  | nil =>
      simp only [gtCompileTail, gtEvalJointBitList]
      exact gt_eval_trailing_square_ops table tableZero (run - 1) out
  | cons index rest ih =>
      rw [gtCompileTail]
      by_cases indexZero : index.val = 0
      · simp only [indexZero, if_true]
        rw [ih (run + 1) (by omega)]
        have tableIndexZero : table index = 0 := by
          have indexEq : index = 0 := Fin.ext indexZero
          rw [indexEq, tableZero]
        simp only [gtEvalJointBitList, tableIndexZero, add_zero]
        have powStep :
            (2 ^ run : GtExponent) = 2 * (2 ^ (run - 1) : GtExponent) := by
          obtain ⟨prior, rfl⟩ := Nat.exists_eq_succ_of_ne_zero (Nat.ne_of_gt runPositive)
          simp [pow_succ, mul_comm]
        rw [show run + 1 - 1 = run by omega, powStep]
        ring_nf
      · simp only [indexZero, if_false, gt_eval_compiled_append]
        rw [gt_eval_run_mul_ops table tableZero run index runPositive out,
          ih 1 (by omega)]
        simp only [Nat.reduceSub, pow_zero, one_mul,
          gtEvalJointBitList]
        have powStep :
            (2 ^ run : GtExponent) = 2 * (2 ^ (run - 1) : GtExponent) := by
          obtain ⟨prior, rfl⟩ := Nat.exists_eq_succ_of_ne_zero (Nat.ne_of_gt runPositive)
          simp [pow_succ, mul_comm]
        rw [powStep]
        ring_nf

theorem gt_eval_drop_leading_zero
    (table : Fin 16 -> GtExponent) (tableZero : table 0 = 0)
    (bits : List (Fin 16)) :
    gtEvalJointBitList table
        (bits.dropWhile (fun index => index.val = 0)) 0 =
      gtEvalJointBitList table bits 0 := by
  induction bits with
  | nil => rfl
  | cons index rest ih =>
      by_cases indexZero : index.val = 0
      · simp only [List.dropWhile, indexZero, decide_true, gtEvalJointBitList]
        have tableIndexZero : table index = 0 := by
          have indexEq : index = 0 := Fin.ext indexZero
          rw [indexEq, tableZero]
        rw [tableIndexZero]
        simpa using ih
      · simp [List.dropWhile, indexZero]

theorem gt_eval_compiled_schedule
    (table : Fin 16 -> GtExponent) (tableZero : table 0 = 0)
    (bits : List (Fin 16)) :
    gtEvalCompiledOps table (gtCompileSchedule bits) 0 =
      gtEvalJointBitList table bits 0 := by
  generalize trimmedEq : bits.dropWhile (fun index => index.val = 0) = trimmed
  have dropCorrect := gt_eval_drop_leading_zero table tableZero bits
  rw [trimmedEq] at dropCorrect
  rw [gtCompileSchedule, trimmedEq]
  cases trimmed with
  | nil => simpa [gtEvalCompiledOps, gtEvalJointBitList] using dropCorrect
  | cons first rest =>
      simp only [gtEvalCompiledOps, gtEvalCompiledOp, pow_zero,
        one_mul, zero_add]
      rw [gt_eval_compile_tail table tableZero rest 1 (by omega)]
      simp only [Nat.reduceSub, pow_zero, one_mul]
      simpa [gtEvalJointBitList] using dropCorrect

theorem gt_joint_bit_list_evaluation
    (table : Fin 16 -> GtExponent) (magnitudes : Fin 4 -> Nat)
    (bits : Nat) (out : GtExponent) :
    gtEvalJointBitList table
        ((List.range bits).reverse.map (gtJointIndex magnitudes)) out =
      gtJointCtLoop table magnitudes bits out := by
  induction bits generalizing out with
  | zero => rfl
  | succ bits ih =>
      rw [List.range_succ, List.reverse_append, List.map_append]
      simp only [List.reverse_singleton, List.map_singleton, List.singleton_append,
        gtEvalJointBitList, gtJointCtLoop]
      exact ih (out + out + table (gtJointIndex magnitudes bits))

theorem gt_compiled_glv4_executor_correct
    (scalar : Nat) (scalarBound : scalar < 2 ^ 256) :
    gtEvalCompiledOps
        (gtSourceTableExponent (gtBalancedSignedBasis scalar))
        (gtCompiledScheduleForScalar scalar) 0 =
      (scalar : GtExponent) := by
  let basis := gtBalancedSignedBasis scalar
  let magnitudes := gtBalancedMagnitudes scalar
  have tableZero : gtSourceTableExponent basis 0 = 0 :=
    gt_source_table_zero basis
  rw [gtCompiledScheduleForScalar, gt_eval_compiled_schedule _ tableZero,
    gtJointBitList, gt_joint_bit_list_evaluation]
  rw [gt_joint_ct_loop_correct]
  simp only [mul_zero, zero_add, gtJointExpected, Fin.sum_univ_four]
  have fit (index : Fin 4) : magnitudes index < 2 ^ 64 :=
    (gt_balanced_magnitude_bound scalar scalarBound index).trans (by omega)
  rw [Nat.mod_eq_of_lt (fit 0), Nat.mod_eq_of_lt (fit 1),
    Nat.mod_eq_of_lt (fit 2), Nat.mod_eq_of_lt (fit 3)]
  simpa [basis, magnitudes, Fin.sum_univ_four] using
    gt_balanced_signed_basis_reconstructs scalar

end CatBLST.AVX.MatrixNoK
