/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.MillerAdd
import CatBLST.AVX.MatrixNoK.Fp12Conjugate
import Mathlib.Tactic

/-!
# Fixed BLS12-381 Miller schedule

The selected source starts with one doubling line, then executes five
`add + n * double` segments with lengths `[2, 3, 9, 32, 16]`, followed by
conjugation for the negative BLS parameter.  This file proves that literal
schedule is the ordinary left-to-right binary schedule for
`|x| = 0xd201000000010000` and composes the accepted point, line, sparse-Fp12,
and conjugation semantics.
-/

namespace CatBLST.AVX.MatrixNoK

inductive MillerScheduleOp where
  | double
  | add
  deriving DecidableEq, Repr

def blsXAbs : Nat := 0xd201000000010000

def cppMillerSegments : List Nat := [2, 3, 9, 32, 16]

def cppMillerScheduleTail : List MillerScheduleOp :=
  cppMillerSegments.flatMap (fun count =>
    MillerScheduleOp.add :: List.replicate count MillerScheduleOp.double)

def cppMillerSchedule : List MillerScheduleOp :=
  MillerScheduleOp.double :: cppMillerScheduleTail

def binaryMillerSchedule : List MillerScheduleOp :=
  (List.range 63).reverse.flatMap (fun bit =>
    MillerScheduleOp.double ::
      if Nat.testBit blsXAbs bit then [MillerScheduleOp.add] else [])

theorem cpp_miller_schedule_matches_bls_x :
    cppMillerSchedule = binaryMillerSchedule := by
  native_decide

def millerPointScalarStep (value : Nat) : MillerScheduleOp -> Nat
  | .double => 2 * value
  | .add => value + 1

theorem cpp_miller_schedule_reaches_bls_x :
    cppMillerSchedule.foldl millerPointScalarStep 1 = blsXAbs := by
  native_decide

theorem cpp_miller_schedule_counts :
    cppMillerSchedule.length = 68 ∧
      cppMillerSchedule.count MillerScheduleOp.double = 63 ∧
      cppMillerSchedule.count MillerScheduleOp.add = 5 := by
  native_decide

structure SelectedMillerState where
  point : HomogeneousPoint SelectedMillerFp2
  accumulator : SelectedFp12

def selectedMillerStep
    (b g1x g1y addendX addendY : SelectedMillerFp2)
    (state : SelectedMillerState) : MillerScheduleOp -> SelectedMillerState
  | .double =>
      let result := cppMillerDouble b state.point
      { point := result.point
        accumulator := state.accumulator ^ 2 *
          (result.line.scaleForG1 g1x g1y).toSparseFp12 }
  | .add =>
      let result := cppMillerAdd b state.point addendX addendY
      { point := result.point
        accumulator := state.accumulator *
          (result.line.scaleForG1 g1x g1y).toSparseFp12 }

def runSelectedMillerSchedule
    (schedule : List MillerScheduleOp)
    (b g1x g1y addendX addendY : SelectedMillerFp2)
    (initial : SelectedMillerState) : SelectedMillerState :=
  schedule.foldl (selectedMillerStep b g1x g1y addendX addendY) initial

def selectedCppMillerLoop
    (b g1x g1y addendX addendY : SelectedMillerFp2)
    (initial : SelectedMillerState) : SelectedFp12 :=
  cppFp12Conjugate
    (runSelectedMillerSchedule cppMillerSchedule
      b g1x g1y addendX addendY initial).accumulator

/-- Abstract binary Miller loop for the fixed BLS12-381 parameter. -/
def selectedAbstractBlsMillerLoop
    (b g1x g1y addendX addendY : SelectedMillerFp2)
    (initial : SelectedMillerState) : SelectedFp12 :=
  cppFp12Conjugate
    (runSelectedMillerSchedule binaryMillerSchedule
      b g1x g1y addendX addendY initial).accumulator

theorem selected_cpp_miller_schedule_correct
    (b g1x g1y addendX addendY : SelectedMillerFp2)
    (initial : SelectedMillerState) :
    selectedCppMillerLoop b g1x g1y addendX addendY initial =
      selectedAbstractBlsMillerLoop b g1x g1y addendX addendY initial := by
  rw [selectedCppMillerLoop, selectedAbstractBlsMillerLoop,
    cpp_miller_schedule_matches_bls_x]

theorem selected_miller_step_preserves_curve
    (b g1x g1y addendX addendY : SelectedMillerFp2)
    (addendOnCurve :
      (⟨addendX, addendY, 1⟩ : HomogeneousPoint SelectedMillerFp2).onCurve b)
    (state : SelectedMillerState) (stateOnCurve : state.point.onCurve b)
    (operation : MillerScheduleOp) :
    (selectedMillerStep b g1x g1y addendX addendY state operation).point.onCurve b := by
  cases operation
  · exact cpp_miller_double_preserves_curve b state.point stateOnCurve
  · exact cpp_miller_add_preserves_curve b state.point addendX addendY
      stateOnCurve addendOnCurve

theorem run_selected_miller_schedule_preserves_curve
    (schedule : List MillerScheduleOp)
    (b g1x g1y addendX addendY : SelectedMillerFp2)
    (addendOnCurve :
      (⟨addendX, addendY, 1⟩ : HomogeneousPoint SelectedMillerFp2).onCurve b)
    (initial : SelectedMillerState) (initialOnCurve : initial.point.onCurve b) :
    (runSelectedMillerSchedule schedule b g1x g1y addendX addendY initial).point.onCurve b := by
  induction schedule generalizing initial with
  | nil => exact initialOnCurve
  | cons operation rest inductionHyp =>
      apply inductionHyp
      exact selected_miller_step_preserves_curve b g1x g1y addendX addendY
        addendOnCurve initial initialOnCurve operation

theorem selected_cpp_miller_final_conjugation
    (b g1x g1y addendX addendY : SelectedMillerFp2)
    (initial : SelectedMillerState) :
    selectedCppMillerLoop b g1x g1y addendX addendY initial =
      (runSelectedMillerSchedule cppMillerSchedule
        b g1x g1y addendX addendY initial).accumulator ^
          (Generated.target_modulus ^ 6) := by
  exact cpp_fp12_conjugate_correct _

end CatBLST.AVX.MatrixNoK
