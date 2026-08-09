/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.GtBalancedSplit
import Mathlib.Tactic

set_option maxHeartbeats 800000

/-!
# Compiled GT GLV4 schedule

This models the high/low nibbles emitted by `avx_glv4_compile`.  A high nibble
is a run of at most fifteen cyclotomic squares; a low nibble is a joint table
index.  Runs longer than fifteen are emitted as square-only opcodes exactly as
in the source.
-/

namespace CatBLST.AVX.MatrixNoK

structure GtCompiledOp where
  squares : Nat
  index : Fin 16
  squaresBound : squares <= 15
  deriving DecidableEq

def GtCompiledOp.byte (op : GtCompiledOp) : Nat :=
  16 * op.squares + op.index.val

theorem gt_compiled_op_fits_byte (op : GtCompiledOp) : op.byte < 256 := by
  have indexBound := op.index.isLt
  have squareBound := op.squaresBound
  simp only [GtCompiledOp.byte]
  omega

theorem gt_compiled_op_high_nibble (op : GtCompiledOp) :
    op.byte / 16 = op.squares := by
  have indexBound := op.index.isLt
  rw [GtCompiledOp.byte, Nat.add_comm]
  rw [Nat.add_mul_div_left _ _ (by omega)]
  rw [Nat.div_eq_of_lt indexBound]
  omega

theorem gt_compiled_op_low_nibble (op : GtCompiledOp) :
    op.byte % 16 = op.index.val := by
  have indexBound := op.index.isLt
  simp [GtCompiledOp.byte, Nat.add_mod, Nat.mod_eq_of_lt indexBound]

def gtSquareOnlyOp : GtCompiledOp :=
  { squares := 15, index := 0, squaresBound := by omega }

def gtRunMulOps (run : Nat) (index : Fin 16) : List GtCompiledOp :=
  let full := (run - 1) / 15
  let last := (run - 1) % 15 + 1
  List.replicate full gtSquareOnlyOp ++
    [{ squares := last, index := index, squaresBound := by
        have remainderBound : (run - 1) % 15 < 15 := by omega
        omega }]

def gtTrailingSquareOps (count : Nat) : List GtCompiledOp :=
  let full := count / 15
  let remainder := count % 15
  List.replicate full gtSquareOnlyOp ++
    if remainder = 0 then [] else
      [{ squares := remainder, index := 0, squaresBound := by
          have remainderBound : count % 15 < 15 := by omega
          omega }]

/-- Exact source tail compiler with its pending `sqr_run`. -/
def gtCompileTail : List (Fin 16) -> Nat -> List GtCompiledOp
  | [], run => gtTrailingSquareOps (run - 1)
  | index :: rest, run =>
      if index.val = 0 then
        gtCompileTail rest (run + 1)
      else
        gtRunMulOps run index ++ gtCompileTail rest 1

def gtSquareSum (ops : List GtCompiledOp) : Nat :=
  (ops.map GtCompiledOp.squares).sum

theorem gt_run_mul_ops_square_sum (run : Nat) (index : Fin 16)
    (runPositive : 0 < run) :
    gtSquareSum (gtRunMulOps run index) = run := by
  have division := Nat.mod_add_div (run - 1) 15
  have arithmetic :
      ((run - 1) / 15) * 15 + ((run - 1) % 15 + 1) = run := by
    omega
  simpa [gtSquareSum, gtRunMulOps, gtSquareOnlyOp, nsmul_eq_mul,
    Nat.mul_comm] using arithmetic

theorem gt_trailing_square_ops_square_sum (count : Nat) :
    gtSquareSum (gtTrailingSquareOps count) = count := by
  have division := Nat.mod_add_div count 15
  by_cases remainderZero : count % 15 = 0
  · simp [gtSquareSum, gtTrailingSquareOps, gtSquareOnlyOp, remainderZero]
    omega
  · simp [gtSquareSum, gtTrailingSquareOps, gtSquareOnlyOp, remainderZero]
    omega

theorem gt_run_mul_ops_positive_squares
    (run : Nat) (index : Fin 16) (_runPositive : 0 < run) :
    (op : GtCompiledOp) -> op ∈ gtRunMulOps run index -> 0 < op.squares := by
  intro op membership
  simp only [gtRunMulOps, List.mem_append, List.mem_replicate,
    List.mem_singleton] at membership
  rcases membership with full | last
  · rcases full with ⟨_, rfl⟩
    simp [gtSquareOnlyOp]
  · subst op
    simp

theorem gt_trailing_square_ops_positive_squares (count : Nat) :
    (op : GtCompiledOp) -> op ∈ gtTrailingSquareOps count -> 0 < op.squares := by
  intro op membership
  by_cases remainderZero : count % 15 = 0
  · simp only [gtTrailingSquareOps, remainderZero, if_true, List.append_nil,
      List.mem_replicate] at membership
    rcases membership with ⟨_, rfl⟩
    simp [gtSquareOnlyOp]
  · simp only [gtTrailingSquareOps, remainderZero, if_false, List.mem_append,
      List.mem_replicate, List.mem_singleton] at membership
    rcases membership with full | last
    · rcases full with ⟨_, rfl⟩
      simp [gtSquareOnlyOp]
    · subst op
      have remainderPositive : 0 < count % 15 := Nat.pos_of_ne_zero remainderZero
      simpa using remainderPositive

theorem gt_compile_tail_square_sum
    (bits : List (Fin 16)) (run : Nat) (runPositive : 0 < run) :
    gtSquareSum (gtCompileTail bits run) = bits.length + run - 1 := by
  induction bits generalizing run with
  | nil =>
      simp only [gtCompileTail, List.length_nil, Nat.zero_add]
      rw [gt_trailing_square_ops_square_sum]
  | cons index rest ih =>
      rw [gtCompileTail]
      by_cases indexZero : index.val = 0
      · simp only [indexZero, if_true, List.length_cons]
        rw [ih (run + 1) (by omega)]
        omega
      · simp only [indexZero, if_false, gtSquareSum, List.map_append,
          List.sum_append]
        change gtSquareSum (gtRunMulOps run index) +
          gtSquareSum (gtCompileTail rest 1) = _
        rw [gt_run_mul_ops_square_sum run index runPositive,
          ih 1 (by omega)]
        simp
        omega

theorem gt_compile_tail_positive_squares
    (bits : List (Fin 16)) (run : Nat) (runPositive : 0 < run) :
    (op : GtCompiledOp) -> op ∈ gtCompileTail bits run -> 0 < op.squares := by
  induction bits generalizing run with
  | nil =>
      simpa [gtCompileTail] using gt_trailing_square_ops_positive_squares (run - 1)
  | cons index rest ih =>
      rw [gtCompileTail]
      by_cases indexZero : index.val = 0
      · simp only [indexZero, if_true]
        exact ih (run + 1) (by omega)
      · simp only [indexZero, if_false]
        intro op membership
        rcases List.mem_append.mp membership with inRun | inRest
        · exact gt_run_mul_ops_positive_squares run index runPositive op inRun
        · exact ih 1 (by omega) op inRest

theorem gt_list_length_le_square_sum
    (ops : List GtCompiledOp)
    (positive : (op : GtCompiledOp) -> op ∈ ops -> 0 < op.squares) :
    ops.length <= gtSquareSum ops := by
  induction ops with
  | nil => simp [gtSquareSum]
  | cons op rest ih =>
      have headPositive := positive op (by simp)
      have restPositive : (item : GtCompiledOp) ->
          item ∈ rest -> 0 < item.squares := by
        intro item membership
        exact positive item (by simp [membership])
      have tailBound := ih restPositive
      simp only [gtSquareSum] at tailBound
      simp only [gtSquareSum, List.map_cons, List.sum_cons, List.length_cons]
      omega

def gtCompileSchedule (bits : List (Fin 16)) : List GtCompiledOp :=
  match bits.dropWhile (fun index => index.val = 0) with
  | [] => []
  | first :: rest =>
      { squares := 0, index := first, squaresBound := by omega } ::
        gtCompileTail rest 1

theorem gt_dropWhile_length_le
    (predicate : Fin 16 -> Bool) (bits : List (Fin 16)) :
    (bits.dropWhile predicate).length <= bits.length := by
  induction bits with
  | nil => simp
  | cons head tail ih =>
      simp only [List.dropWhile]
      split
      all_goals simp_all
      all_goals omega

theorem gt_compile_schedule_length_le_input (bits : List (Fin 16)) :
    (gtCompileSchedule bits).length <= bits.length := by
  generalize trimmedEq : bits.dropWhile (fun index => index.val = 0) = trimmed
  have trimmedBound : trimmed.length <= bits.length := by
    rw [← trimmedEq]
    exact gt_dropWhile_length_le _ _
  rw [gtCompileSchedule, trimmedEq]
  cases trimmed with
  | nil => simp
  | cons first rest =>
      simp only [List.length_cons] at trimmedBound
      have tailPositive := gt_compile_tail_positive_squares rest 1 (by omega)
      have tailLength := gt_list_length_le_square_sum (gtCompileTail rest 1) tailPositive
      have tailSum := gt_compile_tail_square_sum rest 1 (by omega)
      simp only [List.length_cons]
      rw [tailSum] at tailLength
      simp at tailLength
      omega

theorem gt_compile_schedule_square_sum_le_input (bits : List (Fin 16)) :
    gtSquareSum (gtCompileSchedule bits) <= bits.length := by
  generalize trimmedEq : bits.dropWhile (fun index => index.val = 0) = trimmed
  have trimmedBound : trimmed.length <= bits.length := by
    rw [← trimmedEq]
    exact gt_dropWhile_length_le _ _
  rw [gtCompileSchedule, trimmedEq]
  cases trimmed with
  | nil => simp [gtSquareSum]
  | cons first rest =>
      simp only [gtSquareSum, List.map_cons, List.sum_cons]
      change 0 + gtSquareSum (gtCompileTail rest 1) <= bits.length
      rw [gt_compile_tail_square_sum rest 1 (by omega)]
      simp only [Nat.zero_add]
      simp only [List.length_cons] at trimmedBound
      omega

theorem gt_compile_schedule_square_sum_lt_input
    (bits : List (Fin 16)) (bitsNonempty : bits ≠ []) :
    gtSquareSum (gtCompileSchedule bits) < bits.length := by
  generalize trimmedEq : bits.dropWhile (fun index => index.val = 0) = trimmed
  have bitsPositive : 0 < bits.length := by
    apply Nat.pos_of_ne_zero
    simpa using bitsNonempty
  have trimmedBound : trimmed.length <= bits.length := by
    rw [← trimmedEq]
    exact gt_dropWhile_length_le _ _
  rw [gtCompileSchedule, trimmedEq]
  cases trimmed with
  | nil => simpa [gtSquareSum] using bitsPositive
  | cons first rest =>
      simp only [gtSquareSum, List.map_cons, List.sum_cons]
      change 0 + gtSquareSum (gtCompileTail rest 1) < bits.length
      rw [gt_compile_tail_square_sum rest 1 (by omega)]
      simp only [Nat.zero_add]
      simp only [List.length_cons] at trimmedBound
      omega

def gtBalancedMagnitudes (scalar : Nat) : Fin 4 -> Nat :=
  fun index => (gtBalancedSourceDigits scalar index).natAbs

def gtJointBitList (magnitudes : Fin 4 -> Nat) : List (Fin 16) :=
  (List.range 64).reverse.map (gtJointIndex magnitudes)

theorem gt_joint_bit_list_length (magnitudes : Fin 4 -> Nat) :
    (gtJointBitList magnitudes).length = 64 := by
  simp [gtJointBitList]

def gtCompiledScheduleForScalar (scalar : Nat) : List GtCompiledOp :=
  gtCompileSchedule (gtJointBitList (gtBalancedMagnitudes scalar))

theorem gt_compiled_schedule_fits_source_capacity (scalar : Nat) :
    (gtCompiledScheduleForScalar scalar).length <= 64 := by
  apply (gt_compile_schedule_length_le_input _).trans_eq
  exact gt_joint_bit_list_length _

theorem gt_compiled_schedule_fits_array_70 (scalar : Nat) :
    (gtCompiledScheduleForScalar scalar).length < 70 :=
  (gt_compiled_schedule_fits_source_capacity scalar).trans_lt (by omega)

theorem gt_compiled_schedule_square_bound (scalar : Nat) :
    gtSquareSum (gtCompiledScheduleForScalar scalar) <= 64 := by
  apply (gt_compile_schedule_square_sum_le_input _).trans_eq
  exact gt_joint_bit_list_length _

end CatBLST.AVX.MatrixNoK
