/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.MillerReduction
import Mathlib.Tactic

set_option maxHeartbeats 800000

/-!
# Literal-to-abstract MatrixNoK Miller composition

This final P4 layer keeps the source's literal flat Fp12 square and sparse
multiply schedules on one side and the accepted field operations on the
other.  The first doubling deliberately omits the source's no-op `one^2`.
-/

namespace CatBLST.AVX.MatrixNoK

@[ext] theorem SelectedMillerState.ext
    {left right : SelectedMillerState}
    (point : left.point = right.point)
    (accumulator : left.accumulator = right.accumulator) : left = right := by
  cases left
  cases right
  simp_all

def selectedSourceFp12Square (value : SelectedFp12) : SelectedFp12 :=
  packFlatFp12 (cppFp12LiteralSquare 0 (unpackFlatFp12 value))

theorem selected_source_fp12_square_correct (value : SelectedFp12) :
    selectedSourceFp12Square value = value ^ 2 := by
  rw [selectedSourceFp12Square, cpp_fp12_literal_square_packs,
    cpp_fp12_flat_square_correct]
  rw [pow_two]
  rfl

def selectedSourceSparseMul
    (accumulator : SelectedFp12) (line : MillerLine SelectedMillerFp2) :
    SelectedFp12 :=
  packFlatFp12
    (cppFp12LiteralMulXY00Z0 0 (unpackFlatFp12 accumulator)
      line.c line.x line.y)

theorem selected_source_sparse_mul_correct
    (accumulator : SelectedFp12) (line : MillerLine SelectedMillerFp2) :
    selectedSourceSparseMul accumulator line =
      accumulator * line.toSparseFp12 := by
  exact selected_cpp_flat_mulxy00z0_correct accumulator line

def selectedSourceMillerStep
    (b g1x g1y addendX addendY : SelectedMillerFp2)
    (state : SelectedMillerState) : MillerScheduleOp → SelectedMillerState
  | .double =>
      let result := cppMillerDouble b state.point
      { point := result.point
        accumulator := selectedSourceSparseMul
          (selectedSourceFp12Square state.accumulator)
          (result.line.scaleForG1 g1x g1y) }
  | .add =>
      let result := cppMillerAdd b state.point addendX addendY
      { point := result.point
        accumulator := selectedSourceSparseMul state.accumulator
          (result.line.scaleForG1 g1x g1y) }

theorem selected_source_miller_step_correct
    (b g1x g1y addendX addendY : SelectedMillerFp2)
    (state : SelectedMillerState) (operation : MillerScheduleOp) :
    selectedSourceMillerStep b g1x g1y addendX addendY state operation =
      selectedMillerStep b g1x g1y addendX addendY state operation := by
  cases operation <;>
    ext <;>
    simp only [selectedSourceMillerStep, selectedMillerStep] <;>
    rw [selected_source_sparse_mul_correct] <;>
    try rw [selected_source_fp12_square_correct]

def runSelectedSourceMillerSchedule
    (schedule : List MillerScheduleOp)
    (b g1x g1y addendX addendY : SelectedMillerFp2)
    (initial : SelectedMillerState) : SelectedMillerState :=
  schedule.foldl (selectedSourceMillerStep b g1x g1y addendX addendY) initial

theorem run_selected_source_miller_schedule_correct
    (schedule : List MillerScheduleOp)
    (b g1x g1y addendX addendY : SelectedMillerFp2)
    (initial : SelectedMillerState) :
    runSelectedSourceMillerSchedule schedule b g1x g1y addendX addendY initial =
      runSelectedMillerSchedule schedule b g1x g1y addendX addendY initial := by
  induction schedule generalizing initial with
  | nil => rfl
  | cons operation rest inductionHyp =>
      simp only [runSelectedSourceMillerSchedule, runSelectedMillerSchedule,
        List.foldl_cons]
      rw [selected_source_miller_step_correct]
      exact inductionHyp _

def initialSelectedMillerState
    (addendX addendY : SelectedMillerFp2) : SelectedMillerState :=
  { point := ⟨addendX, addendY, 1⟩
    accumulator := 1 }

/-- First source doubling, where squaring the identity is intentionally absent. -/
def selectedSourceFirstDouble
    (b g1x g1y : SelectedMillerFp2)
    (initial : SelectedMillerState) : SelectedMillerState :=
  let result := cppMillerDouble b initial.point
  { point := result.point
    accumulator := selectedSourceSparseMul initial.accumulator
      (result.line.scaleForG1 g1x g1y) }

theorem selected_source_first_double_matches_uniform_step
    (b g1x g1y addendX addendY : SelectedMillerFp2) :
    selectedSourceFirstDouble b g1x g1y
        (initialSelectedMillerState addendX addendY) =
      selectedSourceMillerStep b g1x g1y addendX addendY
        (initialSelectedMillerState addendX addendY) .double := by
  apply SelectedMillerState.ext
  · rfl
  · simp only [selectedSourceFirstDouble, selectedSourceMillerStep,
      initialSelectedMillerState]
    rw [selected_source_fp12_square_correct]
    simp

theorem run_cpp_miller_schedule_after_first
    (b g1x g1y addendX addendY : SelectedMillerFp2)
    (initial : SelectedMillerState) :
    runSelectedMillerSchedule cppMillerScheduleTail
        b g1x g1y addendX addendY
        (selectedMillerStep b g1x g1y addendX addendY initial .double) =
      runSelectedMillerSchedule cppMillerSchedule
        b g1x g1y addendX addendY initial := by
  simp only [runSelectedMillerSchedule, cppMillerSchedule, List.foldl_cons]

def selectedLiteralCppMillerLoop
    (b g1x g1y addendX addendY : SelectedMillerFp2) : SelectedFp12 :=
  let initial := initialSelectedMillerState addendX addendY
  let afterFirst := selectedSourceFirstDouble b g1x g1y initial
  cppFp12Conjugate
    (runSelectedSourceMillerSchedule cppMillerScheduleTail
      b g1x g1y addendX addendY afterFirst).accumulator

theorem selected_literal_cpp_miller_loop_correct
    (b g1x g1y addendX addendY : SelectedMillerFp2) :
    selectedLiteralCppMillerLoop b g1x g1y addendX addendY =
      selectedAbstractBlsMillerLoop b g1x g1y addendX addendY
        (initialSelectedMillerState addendX addendY) := by
  rw [selectedLiteralCppMillerLoop]
  rw [selected_source_first_double_matches_uniform_step]
  rw [run_selected_source_miller_schedule_correct]
  rw [selected_source_miller_step_correct]
  rw [run_cpp_miller_schedule_after_first]
  exact selected_cpp_miller_schedule_correct b g1x g1y addendX addendY _

theorem selected_literal_cpp_miller_point_invariant
    (b g1x g1y addendX addendY : SelectedMillerFp2)
    (addendOnCurve :
      (⟨addendX, addendY, 1⟩ : HomogeneousPoint SelectedMillerFp2).onCurve b) :
    (runSelectedMillerSchedule cppMillerSchedule b g1x g1y addendX addendY
      (initialSelectedMillerState addendX addendY)).point.onCurve b := by
  exact run_selected_miller_schedule_preserves_curve cppMillerSchedule
    b g1x g1y addendX addendY addendOnCurve
    (initialSelectedMillerState addendX addendY) addendOnCurve

end CatBLST.AVX.MatrixNoK
