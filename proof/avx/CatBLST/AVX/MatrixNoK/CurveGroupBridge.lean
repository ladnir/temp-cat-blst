/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.CurveSubgroup
import Mathlib.AlgebraicGeometry.EllipticCurve.Projective.Point
import Mathlib.Tactic

namespace CatBLST.AVX.MatrixNoK

open WeierstrassCurve
open scoped WeierstrassCurve.Projective

def shortWeierstrass [Zero K] (b : K) : WeierstrassCurve K :=
  ⟨0, 0, 0, 0, b⟩

def projectiveCoordinates (point : HomogeneousPoint K) : Fin 3 → K :=
  ![point.x, point.y, point.z]

@[simp] theorem projective_coordinates_zero
    (point : HomogeneousPoint K) : projectiveCoordinates point 0 = point.x :=
  (WeierstrassCurve.Projective.fin3_def_ext point.x point.y point.z).1

@[simp] theorem projective_coordinates_one
    (point : HomogeneousPoint K) : projectiveCoordinates point 1 = point.y :=
  (WeierstrassCurve.Projective.fin3_def_ext point.x point.y point.z).2.1

@[simp] theorem projective_coordinates_two
    (point : HomogeneousPoint K) : projectiveCoordinates point 2 = point.z :=
  (WeierstrassCurve.Projective.fin3_def_ext point.x point.y point.z).2.2

theorem projective_equation_iff [CommRing K] (b : K)
    (point : HomogeneousPoint K) :
    (shortWeierstrass b).toProjective.Equation
        (projectiveCoordinates point) ↔ point.onCurve b := by
  rw [WeierstrassCurve.Projective.equation_iff]
  have coordinates := WeierstrassCurve.Projective.fin3_def_ext
    point.x point.y point.z
  have coordinateX : projectiveCoordinates point 0 = point.x := coordinates.1
  have coordinateY : projectiveCoordinates point 1 = point.y := coordinates.2.1
  have coordinateZ : projectiveCoordinates point 2 = point.z := coordinates.2.2
  rw [coordinateX, coordinateY, coordinateZ]
  simp only [shortWeierstrass, HomogeneousPoint.onCurve]
  constructor <;> intro equation
  · linear_combination equation
  · linear_combination equation

def selectedG1Weierstrass : WeierstrassCurve SelectedCurveFp :=
  shortWeierstrass 4

theorem selected_g1_discriminant_nonzero :
    selectedG1Weierstrass.Δ ≠ 0 := by
  native_decide

instance selectedG1Elliptic : selectedG1Weierstrass.IsElliptic :=
  ⟨selected_g1_discriminant_nonzero.isUnit⟩

def HomogeneousPoint.projectivelyValid [CommRing K]
    (b : K) (point : HomogeneousPoint K) : Prop :=
  point.onCurve b ∧ (point.x ≠ 0 ∨ point.y ≠ 0 ∨ point.z ≠ 0)

theorem selected_g1_valid_nonsingular
    (point : HomogeneousPoint SelectedCurveFp)
    (valid : point.projectivelyValid 4) :
    selectedG1Weierstrass.toProjective.Nonsingular
      (projectiveCoordinates point) := by
  have equation : selectedG1Weierstrass.toProjective.Equation
      (projectiveCoordinates point) :=
    (projective_equation_iff 4 point).mpr valid.1
  by_cases zZero : point.z = 0
  · have coordinateZ : projectiveCoordinates point 2 = 0 := by
      simpa using zZero
    apply (WeierstrassCurve.Projective.nonsingular_of_Z_eq_zero
      coordinateZ).mpr
    constructor
    · exact equation
    · have xCubeZero : point.x ^ 3 = 0 := by
        have equationAtInfinity :=
          (WeierstrassCurve.Projective.equation_of_Z_eq_zero coordinateZ).mp equation
        simpa using equationAtInfinity
      have xZero : point.x = 0 := eq_zero_of_pow_eq_zero xCubeZero
      have yNonzero : point.y ≠ 0 := by
        rcases valid.2 with xNonzero | yNonzero | zNonzero
        · exact False.elim (xNonzero xZero)
        · exact yNonzero
        · exact False.elim (zNonzero zZero)
      right
      simpa [selectedG1Weierstrass, shortWeierstrass, xZero] using
        pow_ne_zero 2 yNonzero
  · have coordinateZ : projectiveCoordinates point 2 ≠ 0 := by
      simpa using zZero
    apply (WeierstrassCurve.Projective.nonsingular_of_Z_ne_zero
      coordinateZ).mpr
    apply (WeierstrassCurve.Affine.equation_iff_nonsingular).mp
    exact (WeierstrassCurve.Projective.equation_of_Z_ne_zero
      coordinateZ).mp equation

def toSelectedG1Point
    (point : HomogeneousPoint SelectedCurveFp)
    (valid : point.projectivelyValid 4) :
    selectedG1Weierstrass.toProjective.Point :=
  ⟨(WeierstrassCurve.Projective.nonsingularLift_iff
    (projectiveCoordinates point)).mpr
      (selected_g1_valid_nonsingular point valid)⟩

theorem mathlib_equiv_of_homogeneous_equivalent
    (left right : HomogeneousPoint SelectedCurveFp)
    (leftValid : left.projectivelyValid 4)
    (rightValid : right.projectivelyValid 4)
    (equivalent : left.equivalent right) :
    projectiveCoordinates left ≈ projectiveCoordinates right := by
  have leftNonsingular := selected_g1_valid_nonsingular left leftValid
  have rightNonsingular := selected_g1_valid_nonsingular right rightValid
  rcases equivalent with atInfinity | finite
  · apply WeierstrassCurve.Projective.equiv_of_Z_eq_zero
      leftNonsingular rightNonsingular
    · simpa using atInfinity.1
    · simpa using atInfinity.2
  · by_cases leftZZero : left.z = 0
    · have leftYNonzero : left.y ≠ 0 := by
        apply WeierstrassCurve.Projective.Y_ne_zero_of_Z_eq_zero leftNonsingular
        simpa using leftZZero
      have rightZZero : right.z = 0 := by
        apply (mul_eq_zero.mp ?_).resolve_left leftYNonzero
        calc
          left.y * right.z = right.y * left.z := finite.2
          _ = 0 := by rw [leftZZero, mul_zero]
      apply WeierstrassCurve.Projective.equiv_of_Z_eq_zero
        leftNonsingular rightNonsingular
      · simpa using leftZZero
      · simpa using rightZZero
    · have rightZNonzero : right.z ≠ 0 := by
        intro rightZZero
        have rightYNonzero : right.y ≠ 0 := by
          apply WeierstrassCurve.Projective.Y_ne_zero_of_Z_eq_zero rightNonsingular
          simpa using rightZZero
        have rightYZero : right.y = 0 := by
          apply (mul_eq_zero.mp ?_).resolve_right leftZZero
          calc
            right.y * left.z = left.y * right.z := finite.2.symm
            _ = 0 := by rw [rightZZero, mul_zero]
        exact rightYNonzero rightYZero
      apply WeierstrassCurve.Projective.equiv_of_X_eq_of_Y_eq
        (by simpa using leftZZero) (by simpa using rightZNonzero)
      · simpa using finite.1
      · simpa using finite.2

theorem homogeneous_equivalent_of_mathlib_equiv
    (left right : HomogeneousPoint SelectedCurveFp)
    (mathlibEquivalent :
      projectiveCoordinates left ≈ projectiveCoordinates right) :
    left.equivalent right := by
  by_cases leftZZero : left.z = 0
  · left
    constructor
    · exact leftZZero
    · have zeroEquivalent :=
        WeierstrassCurve.Projective.Z_eq_zero_of_equiv mathlibEquivalent
      have rightZZero : right.z = 0 := by
        apply (zeroEquivalent.mp ?_)
        simpa using leftZZero
      exact rightZZero
  · right
    constructor
    · simpa using WeierstrassCurve.Projective.X_eq_of_equiv mathlibEquivalent
    · simpa using WeierstrassCurve.Projective.Y_eq_of_equiv mathlibEquivalent

theorem to_selected_g1_point_eq_iff
    (left right : HomogeneousPoint SelectedCurveFp)
    (leftValid : left.projectivelyValid 4)
    (rightValid : right.projectivelyValid 4) :
    toSelectedG1Point left leftValid =
        toSelectedG1Point right rightValid ↔
      left.equivalent right := by
  constructor
  · intro pointEqual
    apply homogeneous_equivalent_of_mathlib_equiv left right
    have classEqual := WeierstrassCurve.Projective.Point.ext_iff.mp pointEqual
    change (⟦projectiveCoordinates left⟧ :
      WeierstrassCurve.Projective.PointClass SelectedCurveFp) =
        ⟦projectiveCoordinates right⟧ at classEqual
    exact Quotient.eq.mp classEqual
  · intro equivalent
    apply WeierstrassCurve.Projective.Point.ext
    change (⟦projectiveCoordinates left⟧ :
      WeierstrassCurve.Projective.PointClass SelectedCurveFp) =
        ⟦projectiveCoordinates right⟧
    exact Quotient.eq.mpr
      (mathlib_equiv_of_homogeneous_equivalent left right
        leftValid rightValid equivalent)

theorem selected_g1_zero_valid :
    (cppProjectiveZero : HomogeneousPoint SelectedCurveFp).projectivelyValid 4 := by
  constructor
  · exact cpp_projective_zero_on_curve 4
  · right
    left
    simp [cppProjectiveZero]

theorem to_selected_g1_point_zero :
    toSelectedG1Point cppProjectiveZero selected_g1_zero_valid = 0 := by
  apply WeierstrassCurve.Projective.Point.ext
  rw [WeierstrassCurve.Projective.Point.zero_point]
  rfl

theorem negate_projectively_valid
    (point : HomogeneousPoint SelectedCurveFp)
    (valid : point.projectivelyValid 4) :
    point.negate.projectivelyValid 4 := by
  constructor
  · exact cpp_negate_preserves_curve 4 point valid.1
  · rcases valid.2 with xNonzero | yNonzero | zNonzero
    · exact Or.inl xNonzero
    · exact Or.inr (Or.inl (neg_ne_zero.mpr yNonzero))
    · exact Or.inr (Or.inr zNonzero)

theorem projective_coordinates_negate
    (point : HomogeneousPoint SelectedCurveFp) :
    projectiveCoordinates point.negate =
      selectedG1Weierstrass.toProjective.neg (projectiveCoordinates point) := by
  funext coordinate
  fin_cases coordinate <;>
    simp [HomogeneousPoint.negate, WeierstrassCurve.Projective.neg,
      WeierstrassCurve.Projective.negY, selectedG1Weierstrass,
      shortWeierstrass]

theorem to_selected_g1_point_negate
    (point : HomogeneousPoint SelectedCurveFp)
    (valid : point.projectivelyValid 4) :
    toSelectedG1Point point.negate (negate_projectively_valid point valid) =
      -toSelectedG1Point point valid := by
  apply WeierstrassCurve.Projective.Point.ext
  rw [WeierstrassCurve.Projective.Point.neg_point]
  change (⟦projectiveCoordinates point.negate⟧ :
      WeierstrassCurve.Projective.PointClass SelectedCurveFp) =
    selectedG1Weierstrass.toProjective.negMap
      ⟦projectiveCoordinates point⟧
  rw [WeierstrassCurve.Projective.negMap_eq]
  rw [projective_coordinates_negate]

end CatBLST.AVX.MatrixNoK
