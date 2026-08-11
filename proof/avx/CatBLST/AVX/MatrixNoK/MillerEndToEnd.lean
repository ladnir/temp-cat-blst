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

/-- CatBLST's affine boundary represents group infinity by two zero
coordinates. This predicate belongs to the API representation layer, not to
the finite affine Miller formulas above. -/
def selectedAffineInfinity (x y : SelectedMillerFp2) : Prop :=
  x = 0 ∧ y = 0

/-- Source-shaped public Miller dispatch, including the observed identity bits
that must be handled before finite affine points are converted to `(x,y,1)`.
The C++ representation observer establishing those bits remains a P0 boundary. -/
def selectedPublicCppMillerLoop
    (g1Infinity g2Infinity : Bool)
    (b g1x g1y addendX addendY : SelectedMillerFp2) : SelectedFp12 :=
  if g1Infinity || g2Infinity then
    1
  else
    selectedLiteralCppMillerLoop b g1x g1y addendX addendY

theorem selected_public_cpp_miller_loop_g1_infinity
    (g2Infinity : Bool)
    (b g1x g1y addendX addendY : SelectedMillerFp2) :
    selectedPublicCppMillerLoop true g2Infinity
      b g1x g1y addendX addendY = 1 := by
  simp [selectedPublicCppMillerLoop]

theorem selected_public_cpp_miller_loop_g2_infinity
    (g1Infinity : Bool)
    (b g1x g1y addendX addendY : SelectedMillerFp2) :
    selectedPublicCppMillerLoop g1Infinity true
      b g1x g1y addendX addendY = 1 := by
  cases g1Infinity <;> simp [selectedPublicCppMillerLoop]

theorem selected_public_cpp_miller_loop_finite
    (b g1x g1y addendX addendY : SelectedMillerFp2) :
    selectedPublicCppMillerLoop false false
      b g1x g1y addendX addendY =
      selectedLiteralCppMillerLoop b g1x g1y addendX addendY := by
  simp [selectedPublicCppMillerLoop]

/-- Sparse identity line stored in all 68 slots when G2 is infinity. -/
def selectedPreparedIdentityLine : MillerLine SelectedMillerFp2 :=
  ⟨1, 0, 0⟩

theorem selected_prepared_identity_line_scales_to_one
    (g1x g1y : SelectedMillerFp2) :
    (selectedPreparedIdentityLine.scaleForG1 g1x g1y).toSparseFp12 = 1 := by
  simp only [selectedPreparedIdentityLine, MillerLine.scaleForG1,
    MillerLine.toSparseFp12, zero_mul]
  change packFlatFp12 (sparseMillerLineFlat
    (1 : SelectedMillerFp2) (0 : SelectedMillerFp2)
    (0 : SelectedMillerFp2)) = Fp12Concrete.one
  rfl

def selectedPreparedIdentityStep
    (g1x g1y : SelectedMillerFp2)
    (accumulator : SelectedFp12) : MillerScheduleOp → SelectedFp12
  | .double => selectedSourceSparseMul
      (selectedSourceFp12Square accumulator)
      (selectedPreparedIdentityLine.scaleForG1 g1x g1y)
  | .add => selectedSourceSparseMul accumulator
      (selectedPreparedIdentityLine.scaleForG1 g1x g1y)

theorem selected_prepared_identity_step_one
    (g1x g1y : SelectedMillerFp2) (operation : MillerScheduleOp) :
    selectedPreparedIdentityStep g1x g1y 1 operation = 1 := by
  cases operation <;>
    simp [selectedPreparedIdentityStep, selected_source_sparse_mul_correct,
      selected_source_fp12_square_correct,
      selected_prepared_identity_line_scales_to_one]

theorem selected_prepared_identity_schedule_one
    (schedule : List MillerScheduleOp) (g1x g1y : SelectedMillerFp2) :
    schedule.foldl (selectedPreparedIdentityStep g1x g1y) 1 = 1 := by
  induction schedule with
  | nil => rfl
  | cons operation rest inductionHyp =>
      simp only [List.foldl_cons]
      rw [selected_prepared_identity_step_one]
      exact inductionHyp

/-- Exact prepared-table evaluator shape: the first identity line replaces
`one^2`, the 67 tail operations consume identity lines, and the result is
conjugated for the negative BLS parameter. -/
def selectedPreparedIdentityMillerLoop
    (g1x g1y : SelectedMillerFp2) : SelectedFp12 :=
  let afterFirst := selectedSourceSparseMul 1
    (selectedPreparedIdentityLine.scaleForG1 g1x g1y)
  cppFp12Conjugate
    (cppMillerScheduleTail.foldl
      (selectedPreparedIdentityStep g1x g1y) afterFirst)

theorem selected_prepared_identity_miller_loop_one
    (g1x g1y : SelectedMillerFp2) :
    selectedPreparedIdentityMillerLoop g1x g1y = 1 := by
  rw [selectedPreparedIdentityMillerLoop]
  simp only [selected_source_sparse_mul_correct,
    selected_prepared_identity_line_scales_to_one, mul_one]
  rw [selected_prepared_identity_schedule_one]
  rw [cpp_fp12_conjugate_correct]
  simp

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
