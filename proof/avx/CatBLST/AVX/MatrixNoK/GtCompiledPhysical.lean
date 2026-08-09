/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.GtCompiledExecutor
import Mathlib.Tactic

/-!
# Physical bounds for compiled GT exponentiation

This closes the source-level schedule boundary: exact byte nibbles, the
64-entry high-to-low discovery scan, the conservative 70-byte owner, and the
operation/storage bounds used by FFT root preparation.
-/

namespace CatBLST.AVX.MatrixNoK

theorem gt_balanced_joint_bit_63_zero
    (scalar : Nat) (scalarBound : scalar < 2 ^ 256) :
    gtJointIndex (gtBalancedMagnitudes scalar) 63 = 0 := by
  have b0 := gt_balanced_magnitude_bound scalar scalarBound 0
  have b1 := gt_balanced_magnitude_bound scalar scalarBound 1
  have b2 := gt_balanced_magnitude_bound scalar scalarBound 2
  have b3 := gt_balanced_magnitude_bound scalar scalarBound 3
  have q0 : gtBalancedMagnitudes scalar 0 / 9223372036854775808 = 0 :=
    Nat.div_eq_of_lt (by simpa using b0)
  have q1 : gtBalancedMagnitudes scalar 1 / 9223372036854775808 = 0 :=
    Nat.div_eq_of_lt (by simpa using b1)
  have q2 : gtBalancedMagnitudes scalar 2 / 9223372036854775808 = 0 :=
    Nat.div_eq_of_lt (by simpa using b2)
  have q3 : gtBalancedMagnitudes scalar 3 / 9223372036854775808 = 0 :=
    Nat.div_eq_of_lt (by simpa using b3)
  apply Fin.ext
  simp [gtJointIndex, gtBitValue, q0, q1, q2, q3]

theorem gt_joint_bit_list_top_cons (magnitudes : Fin 4 -> Nat) :
    gtJointBitList magnitudes =
      gtJointIndex magnitudes 63 ::
        (List.range 63).reverse.map (gtJointIndex magnitudes) := by
  norm_num [gtJointBitList, List.range_succ, List.reverse_append]

theorem gt_compiled_schedule_length_bound_63
    (scalar : Nat) (scalarBound : scalar < 2 ^ 256) :
    (gtCompiledScheduleForScalar scalar).length <= 63 := by
  let rest := (List.range 63).reverse.map
    (gtJointIndex (gtBalancedMagnitudes scalar))
  have decomposition := gt_joint_bit_list_top_cons (gtBalancedMagnitudes scalar)
  have topZero := gt_balanced_joint_bit_63_zero scalar scalarBound
  have strip : gtCompileSchedule (0 :: rest) = gtCompileSchedule rest := by
    simp [gtCompileSchedule, List.dropWhile]
  rw [gtCompiledScheduleForScalar, decomposition, topZero, strip]
  have bound := gt_compile_schedule_length_le_input rest
  have restLength : rest.length = 63 := by simp [rest]
  omega

theorem gt_compiled_schedule_square_bound_62
    (scalar : Nat) (scalarBound : scalar < 2 ^ 256) :
    gtSquareSum (gtCompiledScheduleForScalar scalar) <= 62 := by
  let rest := (List.range 63).reverse.map
    (gtJointIndex (gtBalancedMagnitudes scalar))
  have decomposition := gt_joint_bit_list_top_cons (gtBalancedMagnitudes scalar)
  have topZero := gt_balanced_joint_bit_63_zero scalar scalarBound
  have strip : gtCompileSchedule (0 :: rest) = gtCompileSchedule rest := by
    simp [gtCompileSchedule, List.dropWhile]
  rw [gtCompiledScheduleForScalar, decomposition, topZero, strip]
  have restNonempty : rest ≠ [] := by simp [rest]
  have bound := gt_compile_schedule_square_sum_lt_input rest restNonempty
  have restLength : rest.length = 63 := by simp [rest]
  omega

/-- Source executor multiplications after the table-selected initialization. -/
def gtExecutorMultiplicationCount (ops : List GtCompiledOp) : Nat :=
  ((ops.drop 1).filter fun op => op.index.val != 0).length

theorem gt_executor_multiplication_count_le_tail (ops : List GtCompiledOp) :
    gtExecutorMultiplicationCount ops <= ops.length - 1 := by
  calc
    gtExecutorMultiplicationCount ops <= (ops.drop 1).length := by
      simpa [gtExecutorMultiplicationCount] using
        (List.length_filter_le
          (fun op : GtCompiledOp => op.index.val != 0) (ops.drop 1))
    _ = ops.length - 1 := by simp

theorem gt_compiled_executor_multiplication_bound_62
    (scalar : Nat) (scalarBound : scalar < 2 ^ 256) :
    gtExecutorMultiplicationCount (gtCompiledScheduleForScalar scalar) <= 62 := by
  have countBound := gt_executor_multiplication_count_le_tail
    (gtCompiledScheduleForScalar scalar)
  have lengthBound := gt_compiled_schedule_length_bound_63 scalar scalarBound
  omega

/-- Eleven table products, at most 62 executor products, and normalization. -/
theorem gt_compiled_total_multiplication_bound_74
    (scalar : Nat) (scalarBound : scalar < 2 ^ 256) :
    11 + gtExecutorMultiplicationCount (gtCompiledScheduleForScalar scalar) + 1 <= 74 := by
  have executorBound :=
    gt_compiled_executor_multiplication_bound_62 scalar scalarBound
  omega

structure GtCompiledPhysicalBounds where
  scheduleCapacityBytes : Nat
  maximumScheduleLength : Nat
  tableEntries : Nat
  tableBuildMultiplications : Nat
  frobeniusMaps : Nat
  maximumConjugations : Nat
  maximumExecutorSquares : Nat
  maximumExecutorMultiplications : Nat
  normalizationMultiplications : Nat
  commonElementBound : Nat
  deriving DecidableEq

/-- Exact static sizes and tight source-level maxima for 256-bit inputs. -/
def gtCompiledPhysicalBounds : GtCompiledPhysicalBounds where
  scheduleCapacityBytes := 70
  maximumScheduleLength := 63
  tableEntries := 16
  tableBuildMultiplications := 11
  frobeniusMaps := 3
  maximumConjugations := 4
  maximumExecutorSquares := 62
  maximumExecutorMultiplications := 62
  normalizationMultiplications := 1
  commonElementBound := 6

def GtCompiledPhysicalBounds.maximumTotalMultiplications
    (bounds : GtCompiledPhysicalBounds) : Nat :=
  bounds.tableBuildMultiplications + bounds.maximumExecutorMultiplications +
    bounds.normalizationMultiplications

theorem gt_compiled_exact_physical_bounds :
    gtCompiledPhysicalBounds.scheduleCapacityBytes = 70 ∧
    gtCompiledPhysicalBounds.maximumScheduleLength = 63 ∧
    gtCompiledPhysicalBounds.maximumExecutorSquares = 62 ∧
    gtCompiledPhysicalBounds.maximumTotalMultiplications = 74 ∧
    gtCompiledPhysicalBounds.tableEntries = 16 ∧
    gtCompiledPhysicalBounds.commonElementBound = 6 := by
  decide

theorem gt_scalar_order_lt_two_pow_256 : blsScalarOrder < 2 ^ 256 := by
  native_decide

/-- Every canonical Fr root or inverse scale prepared by the FFT satisfies it. -/
theorem gt_compiled_fft_fr_scalar_correct
    (scalar : Nat) (canonicalFr : scalar < blsScalarOrder) :
    gtEvalCompiledOps
        (gtSourceTableExponent (gtBalancedSignedBasis scalar))
        (gtCompiledScheduleForScalar scalar) 0 =
      (scalar : GtExponent) :=
  gt_compiled_glv4_executor_correct scalar
    (canonicalFr.trans gt_scalar_order_lt_two_pow_256)

end CatBLST.AVX.MatrixNoK
