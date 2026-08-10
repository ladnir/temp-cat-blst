/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.CurveGroupBridge
import Mathlib.Tactic

/-!
# Complete addition on the projective group carrier

This file connects the source-shaped complete-add formula to Mathlib's
projective elliptic-curve group law. The first layer proves the formula's exact
bihomogeneous scaling and reduces arbitrary valid representatives to either the
point at infinity or an affine representative.
-/

namespace CatBLST.AVX.MatrixNoK

open WeierstrassCurve
open scoped WeierstrassCurve.Projective

theorem cpp_complete_add_homogeneous_scale [CommRing K] (b leftScale rightScale : K)
    (left right : HomogeneousPoint K) :
    cppCompleteAdd b (homogeneousScale leftScale left)
        (homogeneousScale rightScale right) =
      homogeneousScale ((leftScale * rightScale) ^ 2)
        (cppCompleteAdd b left right) := by
  rw [cpp_complete_add_coordinates, cpp_complete_add_coordinates]
  ext <;> simp only [homogeneousScale] <;> ring

theorem cpp_complete_add_homogeneous_scale_equivalent [Field K]
    (b leftScale rightScale : K)
    (leftScaleNonzero : leftScale ≠ 0)
    (rightScaleNonzero : rightScale ≠ 0)
    (left right : HomogeneousPoint K) :
    (cppCompleteAdd b (homogeneousScale leftScale left)
      (homogeneousScale rightScale right)).equivalent
        (cppCompleteAdd b left right) := by
  rw [cpp_complete_add_homogeneous_scale]
  apply homogeneous_scale_equivalent
  exact pow_ne_zero 2 (mul_ne_zero leftScaleNonzero rightScaleNonzero)

def affineRepresentative [Field K] (point : HomogeneousPoint K) :
    HomogeneousPoint K :=
  ⟨point.x / point.z, point.y / point.z, 1⟩

theorem homogeneous_scale_affine_representative [Field K]
    (point : HomogeneousPoint K) (zNonzero : point.z ≠ 0) :
    homogeneousScale point.z (affineRepresentative point) = point := by
  ext <;> simp only [homogeneousScale, affineRepresentative]
  · calc
      point.z * (point.x / point.z) = point.z * point.x / point.z := by ring
      _ = point.x := mul_div_cancel_left₀ point.x zNonzero
  · calc
      point.z * (point.y / point.z) = point.z * point.y / point.z := by ring
      _ = point.y := mul_div_cancel_left₀ point.y zNonzero
  · simp

theorem affine_representative_on_curve [Field K] (b : K)
    (point : HomogeneousPoint K) (onCurve : point.onCurve b)
    (zNonzero : point.z ≠ 0) :
    (affineRepresentative point).onCurve b := by
  have scaledOnCurve :
      (homogeneousScale point.z (affineRepresentative point)).onCurve b := by
    rwa [homogeneous_scale_affine_representative point zNonzero]
  simp only [HomogeneousPoint.onCurve, homogeneousScale,
    affineRepresentative] at scaledOnCurve ⊢
  field_simp [zNonzero] at scaledOnCurve ⊢
  exact scaledOnCurve

theorem projectively_valid_at_infinity [Field K] (b : K)
    (point : HomogeneousPoint K) (valid : point.projectivelyValid b)
    (zZero : point.z = 0) :
    point.x = 0 ∧ point.y ≠ 0 := by
  have xCubeZero : point.x ^ 3 = 0 := by
    have onCurve := valid.1
    simp only [HomogeneousPoint.onCurve] at onCurve
    rw [zZero] at onCurve
    simpa using onCurve.symm
  have xZero : point.x = 0 := eq_zero_of_pow_eq_zero xCubeZero
  refine ⟨xZero, ?_⟩
  rcases valid.2 with xNonzero | yNonzero | zNonzero
  · exact False.elim (xNonzero xZero)
  · exact yNonzero
  · exact False.elim (zNonzero zZero)

theorem homogeneous_scale_projective_zero [Field K] (b : K)
    (point : HomogeneousPoint K) (valid : point.projectivelyValid b)
    (zZero : point.z = 0) :
    homogeneousScale point.y cppProjectiveZero = point := by
  obtain ⟨xZero, yNonzero⟩ :=
    projectively_valid_at_infinity b point valid zZero
  ext <;> simp [homogeneousScale, cppProjectiveZero, xZero, zZero]

theorem homogeneous_scale_projectively_valid [Field K] (b factor : K)
    (point : HomogeneousPoint K) (factorNonzero : factor ≠ 0)
    (valid : point.projectivelyValid b) :
    (homogeneousScale factor point).projectivelyValid b := by
  constructor
  · have onCurve := valid.1
    simp only [HomogeneousPoint.onCurve, homogeneousScale] at onCurve ⊢
    linear_combination factor ^ 3 * onCurve
  · rcases valid.2 with xNonzero | yNonzero | zNonzero
    · left
      exact mul_ne_zero factorNonzero xNonzero
    · right
      left
      exact mul_ne_zero factorNonzero yNonzero
    · right
      right
      exact mul_ne_zero factorNonzero zNonzero

theorem to_projective_point_homogeneous_scale [Field K] (b factor : K)
    [elliptic : (shortWeierstrass b).IsElliptic]
    (point : HomogeneousPoint K) (factorNonzero : factor ≠ 0)
    (valid : point.projectivelyValid b) :
    toProjectivePoint b (homogeneousScale factor point)
        (homogeneous_scale_projectively_valid b factor point factorNonzero valid) =
      toProjectivePoint b point valid := by
  apply (to_projective_point_eq_iff b _ _ _ _).mpr
  exact homogeneous_scale_equivalent factor point factorNonzero

theorem to_projective_point_eq_zero_of_z_eq_zero [Field K] (b : K)
    [elliptic : (shortWeierstrass b).IsElliptic]
    (point : HomogeneousPoint K) (valid : point.projectivelyValid b)
    (zZero : point.z = 0) :
    toProjectivePoint b point valid = 0 := by
  have pointEqualsScaledZero :=
    (homogeneous_scale_projective_zero b point valid zZero).symm
  have yNonzero := (projectively_valid_at_infinity b point valid zZero).2
  have pointEquivalentZero : point.equivalent cppProjectiveZero := by
    rw [pointEqualsScaledZero]
    exact homogeneous_scale_equivalent point.y cppProjectiveZero yNonzero
  have pointEqualsZeroRepresentative := (to_projective_point_eq_iff b
    point cppProjectiveZero valid (projective_zero_valid b)).mpr
      pointEquivalentZero
  exact pointEqualsZeroRepresentative.trans (to_projective_point_zero b)

theorem cpp_complete_add_zero_left_eq_scale_y [CommRing K] (b : K)
    (point : HomogeneousPoint K) :
    cppCompleteAdd b cppProjectiveZero point =
      homogeneousScale point.y point := by
  rw [cpp_complete_add_coordinates]
  ext
  · simp [cppProjectiveZero, homogeneousScale]
  · simp [cppProjectiveZero, homogeneousScale]
  · simp [cppProjectiveZero, homogeneousScale]
    ring

theorem cpp_complete_add_zero_right_eq_scale_y [CommRing K] (b : K)
    (point : HomogeneousPoint K) :
    cppCompleteAdd b point cppProjectiveZero =
      homogeneousScale point.y point := by
  rw [cpp_complete_add_commutative]
  exact cpp_complete_add_zero_left_eq_scale_y b point

theorem cpp_complete_add_infinity_left_eq_scale [Field K] (b : K)
    (left right : HomogeneousPoint K) (leftValid : left.projectivelyValid b)
    (leftAtInfinity : left.z = 0) :
    cppCompleteAdd b left right =
      homogeneousScale (left.y ^ 2 * right.y) right := by
  obtain ⟨leftXZero, _⟩ :=
    projectively_valid_at_infinity b left leftValid leftAtInfinity
  rw [cpp_complete_add_coordinates]
  ext <;>
    simp [leftXZero, leftAtInfinity, homogeneousScale] <;> ring

theorem cpp_complete_add_infinity_right_eq_scale [Field K] (b : K)
    (left right : HomogeneousPoint K) (rightValid : right.projectivelyValid b)
    (rightAtInfinity : right.z = 0) :
    cppCompleteAdd b left right =
      homogeneousScale (right.y ^ 2 * left.y) left := by
  rw [cpp_complete_add_commutative]
  exact cpp_complete_add_infinity_left_eq_scale b right left
    rightValid rightAtInfinity

theorem affine_representative_valid [Field K] (b : K)
    (point : HomogeneousPoint K) (valid : point.projectivelyValid b)
    (zNonzero : point.z ≠ 0) :
    (affineRepresentative point).projectivelyValid b := by
  constructor
  · exact affine_representative_on_curve b point valid.1 zNonzero
  · right
    right
    simp [affineRepresentative]

def affinePoint [One K] (x y : K) : HomogeneousPoint K := ⟨x, y, 1⟩

theorem affine_point_valid [Field K] (b x y : K)
    (onCurve : y ^ 2 = x ^ 3 + b) :
    (affinePoint x y).projectivelyValid b := by
  constructor
  · simpa [affinePoint, HomogeneousPoint.onCurve] using onCurve
  · right
    right
    simp [affinePoint]

def shortSecantX [Field K] (x₁ y₁ x₂ y₂ : K) : K :=
  ((y₁ - y₂) / (x₁ - x₂)) ^ 2 - x₁ - x₂

def shortSecantY [Field K] (x₁ y₁ x₂ y₂ : K) : K :=
  -((y₁ - y₂) / (x₁ - x₂) *
      (shortSecantX x₁ y₁ x₂ y₂ - x₁) + y₁)

def shortDifferenceX [Field K] (x₁ y₁ x₂ y₂ : K) : K :=
  ((y₁ + y₂) / (x₁ - x₂)) ^ 2 - x₁ - x₂

def shortDifferenceY [Field K] (x₁ y₁ x₂ y₂ : K) : K :=
  -((y₁ + y₂) / (x₁ - x₂) *
      (shortDifferenceX x₁ y₁ x₂ y₂ - x₁) + y₁)

def shortTangentX [Field K] (x y : K) : K :=
  ((3 * x ^ 2) / (2 * y)) ^ 2 - 2 * x

def shortTangentY [Field K] (x y : K) : K :=
  -((3 * x ^ 2) / (2 * y) * (shortTangentX x y - x) + y)

/-- The short Weierstrass curve `y² = x³ + b` has no affine point of order two.
For this curve shape, such a point would have `y = 0`, hence an `x` satisfying
`x³ + b = 0`. -/
def NoAffineTwoTorsion [Field K] (b : K) : Prop :=
  ∀ x : K, x ^ 3 + b ≠ 0

theorem affine_y_ne_zero_of_no_two_torsion [Field K] (b x y : K)
    (noTwoTorsion : NoAffineTwoTorsion b)
    (onCurve : y ^ 2 = x ^ 3 + b) : y ≠ 0 := by
  intro yZero
  apply noTwoTorsion x
  rw [← onCurve, yZero]
  simp

theorem projectively_valid_y_ne_zero [Field K] (b : K)
    (noTwoTorsion : NoAffineTwoTorsion b)
    (point : HomogeneousPoint K) (valid : point.projectivelyValid b) :
    point.y ≠ 0 := by
  by_cases zZero : point.z = 0
  · exact (projectively_valid_at_infinity b point valid zZero).2
  · have affineValid := affine_representative_valid b point valid zZero
    have affineYNonzero := affine_y_ne_zero_of_no_two_torsion b
      (affineRepresentative point).x (affineRepresentative point).y
      noTwoTorsion (by
        simpa [affinePoint, affineRepresentative, HomogeneousPoint.onCurve] using
          affineValid.1)
    intro yZero
    apply affineYNonzero
    simp [affineRepresentative, yZero]

theorem affine_representative_equation [Field K] (b : K)
    (point : HomogeneousPoint K) (valid : point.projectivelyValid b)
    (zNonzero : point.z ≠ 0) :
    (affineRepresentative point).y ^ 2 =
      (affineRepresentative point).x ^ 3 + b := by
  have onCurve := affine_representative_on_curve b point valid.1 zNonzero
  simpa [HomogeneousPoint.onCurve, affineRepresentative] using onCurve

theorem y_ne_neg_of_two_ne_zero [Field K] (y : K)
    (twoNonzero : (2 : K) ≠ 0) (yNonzero : y ≠ 0) : y ≠ -y := by
  intro equalNeg
  apply yNonzero
  apply (mul_eq_zero.mp
    (show (2 : K) * y = 0 by linear_combination equalNeg)).resolve_left
  exact twoNonzero

theorem short_tangent_slope [Field K] [DecidableEq K] (b x y : K)
    (twoNonzero : (2 : K) ≠ 0) (yNonzero : y ≠ 0) :
    (shortWeierstrass b).toAffine.slope x x y y =
      (3 * x ^ 2) / (2 * y) := by
  rw [WeierstrassCurve.Affine.slope_of_Y_ne rfl]
  · simp [shortWeierstrass]
    ring
  · simpa [shortWeierstrass] using
      y_ne_neg_of_two_ne_zero y twoNonzero yNonzero

theorem short_weierstrass_affine_equation_iff [Field K] (b x y : K) :
    (shortWeierstrass b).toAffine.Equation x y ↔ y ^ 2 = x ^ 3 + b := by
  simp [shortWeierstrass, WeierstrassCurve.Affine.equation_iff]

theorem short_difference_on_curve [Field K] [DecidableEq K]
    (b x₁ y₁ x₂ y₂ : K)
    (leftOnCurve : y₁ ^ 2 = x₁ ^ 3 + b)
    (rightOnCurve : y₂ ^ 2 = x₂ ^ 3 + b) (xDistinct : x₁ ≠ x₂) :
    shortDifferenceY x₁ y₁ x₂ y₂ ^ 2 =
      shortDifferenceX x₁ y₁ x₂ y₂ ^ 3 + b := by
  have leftEquation : (shortWeierstrass b).toAffine.Equation x₁ y₁ :=
    (short_weierstrass_affine_equation_iff b x₁ y₁).mpr leftOnCurve
  have rightEquation : (shortWeierstrass b).toAffine.Equation x₂ (-y₂) := by
    apply (short_weierstrass_affine_equation_iff b x₂ (-y₂)).mpr
    simpa using rightOnCurve
  have differenceEquation := WeierstrassCurve.Affine.equation_add
    leftEquation rightEquation (by
      intro exceptional
      exact xDistinct exceptional.1)
  rw [WeierstrassCurve.Affine.slope_of_X_ne xDistinct] at differenceEquation
  apply (short_weierstrass_affine_equation_iff b
    (shortDifferenceX x₁ y₁ x₂ y₂)
    (shortDifferenceY x₁ y₁ x₂ y₂)).mp
  simpa [shortWeierstrass, shortDifferenceX, shortDifferenceY,
    WeierstrassCurve.Affine.addY, WeierstrassCurve.Affine.negY] using
      differenceEquation

theorem short_secant_on_curve [Field K] [DecidableEq K]
    (b x₁ y₁ x₂ y₂ : K)
    (leftOnCurve : y₁ ^ 2 = x₁ ^ 3 + b)
    (rightOnCurve : y₂ ^ 2 = x₂ ^ 3 + b) (xDistinct : x₁ ≠ x₂) :
    shortSecantY x₁ y₁ x₂ y₂ ^ 2 =
      shortSecantX x₁ y₁ x₂ y₂ ^ 3 + b := by
  have leftEquation : (shortWeierstrass b).toAffine.Equation x₁ y₁ :=
    (short_weierstrass_affine_equation_iff b x₁ y₁).mpr leftOnCurve
  have rightEquation : (shortWeierstrass b).toAffine.Equation x₂ y₂ :=
    (short_weierstrass_affine_equation_iff b x₂ y₂).mpr rightOnCurve
  have sumEquation := WeierstrassCurve.Affine.equation_add
    leftEquation rightEquation (by
      intro exceptional
      exact xDistinct exceptional.1)
  rw [WeierstrassCurve.Affine.slope_of_X_ne xDistinct] at sumEquation
  apply (short_weierstrass_affine_equation_iff b
    (shortSecantX x₁ y₁ x₂ y₂)
    (shortSecantY x₁ y₁ x₂ y₂)).mp
  simpa [shortWeierstrass, shortSecantX, shortSecantY,
    WeierstrassCurve.Affine.addY, WeierstrassCurve.Affine.negY] using
      sumEquation

theorem short_tangent_on_curve [Field K] [DecidableEq K]
    (b x y : K) (twoNonzero : (2 : K) ≠ 0)
    (onCurve : y ^ 2 = x ^ 3 + b) (yNonzero : y ≠ 0) :
    shortTangentY x y ^ 2 = shortTangentX x y ^ 3 + b := by
  have equation : (shortWeierstrass b).toAffine.Equation x y :=
    (short_weierstrass_affine_equation_iff b x y).mpr onCurve
  have sumEquation := WeierstrassCurve.Affine.equation_add equation equation (by
    intro exceptional
    exact (y_ne_neg_of_two_ne_zero y twoNonzero yNonzero) (by
      simpa [shortWeierstrass] using exceptional.2))
  rw [short_tangent_slope b x y twoNonzero yNonzero] at sumEquation
  apply (short_weierstrass_affine_equation_iff b
    (shortTangentX x y) (shortTangentY x y)).mp
  convert sumEquation using 1 <;>
    simp only [shortWeierstrass, shortTangentX, shortTangentY,
      WeierstrassCurve.Affine.addX, WeierstrassCurve.Affine.addY,
      WeierstrassCurve.Affine.negAddY, WeierstrassCurve.Affine.negY,
      zero_mul, sub_zero] <;> ring

def affineXLeftCertificate [CommRing K] (b x₁ y₁ x₂ y₂ : K) : K :=
  -3 * b * y₁ + 2 * b * y₂ - 3 * x₁ ^ 2 * x₂ * y₂ -
    3 * x₁ * x₂ ^ 2 * y₁ + 3 * x₁ * x₂ ^ 2 * y₂ +
    2 * x₂ ^ 3 * y₂ - y₁ ^ 2 * y₂ + y₁ * y₂ ^ 2 + y₂ ^ 3

def affineXRightCertificate [CommRing K] (b x₁ y₁ x₂ y₂ : K) : K :=
  3 * b * y₁ - 2 * b * y₂ + 3 * x₁ ^ 3 * y₁ + x₁ ^ 3 * y₂ +
    3 * x₁ ^ 2 * x₂ * y₁ - 3 * x₁ ^ 2 * x₂ * y₂ -
    3 * x₁ * x₂ ^ 2 * y₁ - y₁ * y₂ ^ 2

def affineYLeftCertificate [CommRing K] (b x₁ y₁ x₂ y₂ : K) : K :=
  3 * b ^ 2 + 12 * b * x₁ * x₂ ^ 2 - 6 * b * x₂ ^ 3 +
    3 * b * y₁ ^ 2 - 5 * b * y₁ * y₂ - 2 * b * y₂ ^ 2 +
    9 * x₁ ^ 2 * x₂ ^ 4 + 3 * x₁ ^ 2 * x₂ * y₁ * y₂ -
    15 * x₁ ^ 2 * x₂ * y₂ ^ 2 - 6 * x₁ * x₂ ^ 5 +
    3 * x₁ * x₂ ^ 2 * y₁ ^ 2 - 6 * x₁ * x₂ ^ 2 * y₁ * y₂ +
    15 * x₁ * x₂ ^ 2 * y₂ ^ 2 - 2 * x₂ ^ 3 * y₁ * y₂ -
    2 * x₂ ^ 3 * y₂ ^ 2 + y₁ ^ 3 * y₂ - 2 * y₁ ^ 2 * y₂ ^ 2 +
    2 * y₂ ^ 4

def affineYRightCertificate [CommRing K] (b x₁ y₁ x₂ y₂ : K) : K :=
  -3 * b ^ 2 + 6 * b * x₁ ^ 3 - 27 * b * x₁ ^ 2 * x₂ +
    15 * b * x₁ * x₂ ^ 2 + 5 * b * y₁ * y₂ - b * y₂ ^ 2 -
    9 * x₁ ^ 5 * x₂ + 6 * x₁ ^ 4 * x₂ ^ 2 +
    2 * x₁ ^ 3 * y₁ * y₂ + 2 * x₁ ^ 3 * y₂ ^ 2 +
    6 * x₁ ^ 2 * x₂ * y₁ * y₂ - 3 * x₁ ^ 2 * x₂ * y₂ ^ 2 -
    3 * x₁ * x₂ ^ 2 * y₁ * y₂ - y₁ * y₂ ^ 3

theorem cpp_complete_add_affine_x_cross [CommRing K] (b x₁ y₁ x₂ y₂ : K)
    (leftOnCurve : y₁ ^ 2 = x₁ ^ 3 + b)
    (rightOnCurve : y₂ ^ 2 = x₂ ^ 3 + b) :
    (x₁ - x₂) ^ 2 * (cppCompleteAdd b (affinePoint x₁ y₁)
      (affinePoint x₂ y₂)).x =
      ((y₁ - y₂) ^ 2 - (x₁ + x₂) * (x₁ - x₂) ^ 2) *
        (cppCompleteAdd b (affinePoint x₁ y₁) (affinePoint x₂ y₂)).z := by
  rw [cpp_complete_add_coordinates]
  simp only [affinePoint, one_mul, mul_one]
  have leftEquation : y₁ ^ 2 - x₁ ^ 3 - b = 0 := by
    linear_combination leftOnCurve
  have rightEquation : y₂ ^ 2 - x₂ ^ 3 - b = 0 := by
    linear_combination rightOnCurve
  linear_combination (norm :=
    (simp only [affineXLeftCertificate, affineXRightCertificate]; ring))
    affineXLeftCertificate b x₁ y₁ x₂ y₂ * leftEquation +
      affineXRightCertificate b x₁ y₁ x₂ y₂ * rightEquation

theorem cpp_complete_add_affine_y_cross [CommRing K] (b x₁ y₁ x₂ y₂ : K)
    (leftOnCurve : y₁ ^ 2 = x₁ ^ 3 + b)
    (rightOnCurve : y₂ ^ 2 = x₂ ^ 3 + b) :
    (x₁ - x₂) ^ 3 * (cppCompleteAdd b (affinePoint x₁ y₁)
      (affinePoint x₂ y₂)).y =
      (-(y₁ * (x₁ - x₂) ^ 3) - (y₁ - y₂) *
        ((y₁ - y₂) ^ 2 - (x₁ + x₂) * (x₁ - x₂) ^ 2 -
          x₁ * (x₁ - x₂) ^ 2)) *
        (cppCompleteAdd b (affinePoint x₁ y₁) (affinePoint x₂ y₂)).z := by
  rw [cpp_complete_add_coordinates]
  simp only [affinePoint, one_mul, mul_one]
  have leftEquation : y₁ ^ 2 - x₁ ^ 3 - b = 0 := by
    linear_combination leftOnCurve
  have rightEquation : y₂ ^ 2 - x₂ ^ 3 - b = 0 := by
    linear_combination rightOnCurve
  linear_combination (norm :=
    (simp only [affineYLeftCertificate, affineYRightCertificate]; ring))
    affineYLeftCertificate b x₁ y₁ x₂ y₂ * leftEquation +
      affineYRightCertificate b x₁ y₁ x₂ y₂ * rightEquation

theorem cpp_complete_add_affine_z_difference [CommRing K] (b x₁ y₁ x₂ y₂ : K)
    (leftOnCurve : y₁ ^ 2 = x₁ ^ 3 + b)
    (rightOnCurve : y₂ ^ 2 = x₂ ^ 3 + b) :
    (cppCompleteAdd b (affinePoint x₁ y₁) (affinePoint x₂ y₂)).z +
      (-(y₁ * (x₁ - x₂) ^ 3) - (y₁ + y₂) *
        ((y₁ + y₂) ^ 2 - (x₁ + x₂) * (x₁ - x₂) ^ 2 -
          x₁ * (x₁ - x₂) ^ 2)) = 0 := by
  rw [cpp_complete_add_coordinates]
  simp only [affinePoint, one_mul, mul_one]
  have leftEquation : y₁ ^ 2 - x₁ ^ 3 - b = 0 := by
    linear_combination leftOnCurve
  have rightEquation : y₂ ^ 2 - x₂ ^ 3 - b = 0 := by
    linear_combination rightOnCurve
  linear_combination
    (-y₁ - 2 * y₂) * leftEquation + (-2 * y₁ - y₂) * rightEquation

theorem cpp_complete_add_affine_z_ne_zero [Field K] [DecidableEq K]
    (b x₁ y₁ x₂ y₂ : K) (noTwoTorsion : NoAffineTwoTorsion b)
    (leftOnCurve : y₁ ^ 2 = x₁ ^ 3 + b)
    (rightOnCurve : y₂ ^ 2 = x₂ ^ 3 + b) (xDistinct : x₁ ≠ x₂) :
    (cppCompleteAdd b (affinePoint x₁ y₁) (affinePoint x₂ y₂)).z ≠ 0 := by
  let differenceY := shortDifferenceY x₁ y₁ x₂ y₂
  have differenceOnCurve := short_difference_on_curve b x₁ y₁ x₂ y₂
    leftOnCurve rightOnCurve xDistinct
  have zRelation := cpp_complete_add_affine_z_difference b x₁ y₁ x₂ y₂
    leftOnCurve rightOnCurve
  have numeratorRelation :
      (x₁ - x₂) ^ 3 * differenceY =
        -(y₁ * (x₁ - x₂) ^ 3) - (y₁ + y₂) *
          ((y₁ + y₂) ^ 2 - (x₁ + x₂) * (x₁ - x₂) ^ 2 -
            x₁ * (x₁ - x₂) ^ 2) := by
    dsimp [differenceY, shortDifferenceY, shortDifferenceX]
    field_simp [sub_ne_zero.mpr xDistinct]
    ring
  intro zZero
  have scaledYZero : (x₁ - x₂) ^ 3 * differenceY = 0 := by
    rw [numeratorRelation]
    linear_combination zRelation - zZero
  have differenceYZero : differenceY = 0 :=
    (mul_eq_zero.mp scaledYZero).resolve_left
      (pow_ne_zero 3 (sub_ne_zero.mpr xDistinct))
  apply noTwoTorsion (shortDifferenceX x₁ y₁ x₂ y₂)
  rw [← differenceOnCurve]
  have shortDifferenceYZero : shortDifferenceY x₁ y₁ x₂ y₂ = 0 := by
    simpa [differenceY] using differenceYZero
  rw [shortDifferenceYZero]
  simp

theorem cpp_complete_add_affine_eq_scaled_secant [Field K] (b x₁ y₁ x₂ y₂ : K)
    (leftOnCurve : y₁ ^ 2 = x₁ ^ 3 + b)
    (rightOnCurve : y₂ ^ 2 = x₂ ^ 3 + b) (xDistinct : x₁ ≠ x₂) :
    cppCompleteAdd b (affinePoint x₁ y₁) (affinePoint x₂ y₂) =
      homogeneousScale
        (cppCompleteAdd b (affinePoint x₁ y₁) (affinePoint x₂ y₂)).z
        (affinePoint (shortSecantX x₁ y₁ x₂ y₂)
          (shortSecantY x₁ y₁ x₂ y₂)) := by
  apply HomogeneousPoint.ext
  · have cross := cpp_complete_add_affine_x_cross b x₁ y₁ x₂ y₂
      leftOnCurve rightOnCurve
    change (cppCompleteAdd b (affinePoint x₁ y₁) (affinePoint x₂ y₂)).x =
      (cppCompleteAdd b (affinePoint x₁ y₁) (affinePoint x₂ y₂)).z *
        shortSecantX x₁ y₁ x₂ y₂
    rw [shortSecantX]
    field_simp [sub_ne_zero.mpr xDistinct]
    convert cross using 1 <;> ring
  · have cross := cpp_complete_add_affine_y_cross b x₁ y₁ x₂ y₂
      leftOnCurve rightOnCurve
    change (cppCompleteAdd b (affinePoint x₁ y₁) (affinePoint x₂ y₂)).y =
      (cppCompleteAdd b (affinePoint x₁ y₁) (affinePoint x₂ y₂)).z *
        shortSecantY x₁ y₁ x₂ y₂
    rw [shortSecantY, shortSecantX]
    field_simp [sub_ne_zero.mpr xDistinct]
    convert cross using 1 <;> ring
  · simp [homogeneousScale, affinePoint]

theorem cpp_complete_add_affine_equivalent_secant [Field K] [DecidableEq K]
    (b x₁ y₁ x₂ y₂ : K) (noTwoTorsion : NoAffineTwoTorsion b)
    (leftOnCurve : y₁ ^ 2 = x₁ ^ 3 + b)
    (rightOnCurve : y₂ ^ 2 = x₂ ^ 3 + b) (xDistinct : x₁ ≠ x₂) :
    (cppCompleteAdd b (affinePoint x₁ y₁) (affinePoint x₂ y₂)).equivalent
      (affinePoint (shortSecantX x₁ y₁ x₂ y₂)
        (shortSecantY x₁ y₁ x₂ y₂)) := by
  rw [cpp_complete_add_affine_eq_scaled_secant b x₁ y₁ x₂ y₂
    leftOnCurve rightOnCurve xDistinct]
  exact homogeneous_scale_equivalent
    (cppCompleteAdd b (affinePoint x₁ y₁) (affinePoint x₂ y₂)).z
    (affinePoint (shortSecantX x₁ y₁ x₂ y₂)
      (shortSecantY x₁ y₁ x₂ y₂))
    (cpp_complete_add_affine_z_ne_zero b x₁ y₁ x₂ y₂ noTwoTorsion
      leftOnCurve rightOnCurve xDistinct)

theorem cpp_public_double_affine_eq_scaled_tangent [Field K]
    (b x y : K) (onCurve : y ^ 2 = x ^ 3 + b)
    (twoNonzero : (2 : K) ≠ 0) (yNonzero : y ≠ 0) :
    cppPublicDouble b (affinePoint x y) =
      homogeneousScale (cppPublicDouble b (affinePoint x y)).z
        (affinePoint (shortTangentX x y) (shortTangentY x y)) := by
  have curveEquation : y ^ 2 - x ^ 3 - b = 0 := by
    linear_combination onCurve
  apply HomogeneousPoint.ext
  · simp only [cppPublicDouble, affinePoint, one_pow, mul_one]
    dsimp [shortTangentX, homogeneousScale]
    field_simp [twoNonzero, yNonzero]
    linear_combination (72 * x) * curveEquation
  · simp only [cppPublicDouble, affinePoint, one_pow, mul_one]
    dsimp [shortTangentY, shortTangentX, homogeneousScale]
    field_simp [twoNonzero, yNonzero]
    linear_combination
      (-72 * (-3 * b + 3 * x ^ 3 - y ^ 2)) * curveEquation
  · simp [cppPublicDouble, affinePoint, homogeneousScale]

theorem cpp_public_double_affine_z_ne_zero [Field K] (b x y : K)
    (twoNonzero : (2 : K) ≠ 0) (yNonzero : y ≠ 0) :
    (cppPublicDouble b (affinePoint x y)).z ≠ 0 := by
  simp only [cppPublicDouble, affinePoint, one_pow, mul_one]
  have eightNonzero : (8 : K) ≠ 0 := by
    have eightAsPower : (8 : K) = 2 ^ 3 := by norm_num
    rw [eightAsPower]
    exact pow_ne_zero 3 twoNonzero
  exact mul_ne_zero (mul_ne_zero eightNonzero (pow_ne_zero 2 yNonzero)) yNonzero

theorem cpp_public_double_affine_valid [Field K] (b x y : K)
    (onCurve : y ^ 2 = x ^ 3 + b)
    (twoNonzero : (2 : K) ≠ 0) (yNonzero : y ≠ 0) :
    (cppPublicDouble b (affinePoint x y)).projectivelyValid b := by
  constructor
  · exact cpp_public_double_preserves_curve b (affinePoint x y)
      (affine_point_valid b x y onCurve).1
  · right
    right
    exact cpp_public_double_affine_z_ne_zero b x y twoNonzero yNonzero

theorem cpp_complete_add_affine_inverse_y_eq_double_y [Field K]
    (b x y : K) (onCurve : y ^ 2 = x ^ 3 + b) :
    (cppCompleteAdd b (affinePoint x y) (affinePoint x (-y))).y =
      (cppPublicDouble b (affinePoint x y)).y := by
  rw [cpp_complete_add_coordinates]
  simp only [affinePoint, cppPublicDouble, one_mul, mul_one, one_pow]
  have curveEquation : y ^ 2 - x ^ 3 - b = 0 := by
    linear_combination onCurve
  linear_combination (-18 * b) * curveEquation

theorem cpp_public_double_affine_y_ne_zero [Field K] [DecidableEq K]
    (b x y : K) (noTwoTorsion : NoAffineTwoTorsion b)
    (twoNonzero : (2 : K) ≠ 0) (onCurve : y ^ 2 = x ^ 3 + b) :
    (cppPublicDouble b (affinePoint x y)).y ≠ 0 := by
  let yNonzero := affine_y_ne_zero_of_no_two_torsion b x y
    noTwoTorsion onCurve
  have tangentOnCurve := short_tangent_on_curve b x y twoNonzero
    onCurve yNonzero
  have tangentYNonzero := affine_y_ne_zero_of_no_two_torsion b
    (shortTangentX x y) (shortTangentY x y) noTwoTorsion tangentOnCurve
  have scaleNonzero := cpp_public_double_affine_z_ne_zero b x y
    twoNonzero yNonzero
  have coordinateEquality := congrArg HomogeneousPoint.y
    (cpp_public_double_affine_eq_scaled_tangent b x y onCurve
      twoNonzero yNonzero)
  simp only [homogeneousScale, affinePoint] at coordinateEquality
  have coordinateEquality' :
      (cppPublicDouble b (affinePoint x y)).y =
        (cppPublicDouble b (affinePoint x y)).z * shortTangentY x y := by
    simpa [affinePoint] using coordinateEquality
  rw [coordinateEquality']
  exact mul_ne_zero scaleNonzero tangentYNonzero

theorem cpp_complete_add_affine_inverse_valid [Field K] [DecidableEq K]
    (b x y : K) (noTwoTorsion : NoAffineTwoTorsion b)
    (twoNonzero : (2 : K) ≠ 0) (onCurve : y ^ 2 = x ^ 3 + b) :
    (cppCompleteAdd b (affinePoint x y)
      (affinePoint x (-y))).projectivelyValid b := by
  constructor
  · apply cpp_complete_add_preserves_curve b
    · exact (affine_point_valid b x y onCurve).1
    · exact (affine_point_valid b x (-y) (by simpa using onCurve)).1
  · right
    left
    rw [cpp_complete_add_affine_inverse_y_eq_double_y b x y onCurve]
    exact cpp_public_double_affine_y_ne_zero b x y noTwoTorsion
      twoNonzero onCurve

theorem cpp_complete_add_affine_valid [Field K] [DecidableEq K]
    (b x₁ y₁ x₂ y₂ : K) (noTwoTorsion : NoAffineTwoTorsion b)
    (leftOnCurve : y₁ ^ 2 = x₁ ^ 3 + b)
    (rightOnCurve : y₂ ^ 2 = x₂ ^ 3 + b) (xDistinct : x₁ ≠ x₂) :
    (cppCompleteAdd b (affinePoint x₁ y₁)
      (affinePoint x₂ y₂)).projectivelyValid b := by
  constructor
  · exact cpp_complete_add_preserves_curve b
      (affinePoint x₁ y₁) (affinePoint x₂ y₂)
      (affine_point_valid b x₁ y₁ leftOnCurve).1
      (affine_point_valid b x₂ y₂ rightOnCurve).1
  · right
    right
    exact cpp_complete_add_affine_z_ne_zero b x₁ y₁ x₂ y₂ noTwoTorsion
      leftOnCurve rightOnCurve xDistinct

theorem projective_affine_secant_eq_add [Field K] [DecidableEq K]
    (b x₁ y₁ x₂ y₂ : K) [elliptic : (shortWeierstrass b).IsElliptic]
    (leftOnCurve : y₁ ^ 2 = x₁ ^ 3 + b)
    (rightOnCurve : y₂ ^ 2 = x₂ ^ 3 + b) (xDistinct : x₁ ≠ x₂) :
    toProjectivePoint b
        (affinePoint (shortSecantX x₁ y₁ x₂ y₂)
          (shortSecantY x₁ y₁ x₂ y₂))
        (affine_point_valid b _ _
          (short_secant_on_curve b x₁ y₁ x₂ y₂
            leftOnCurve rightOnCurve xDistinct)) =
      toProjectivePoint b (affinePoint x₁ y₁)
          (affine_point_valid b x₁ y₁ leftOnCurve) +
        toProjectivePoint b (affinePoint x₂ y₂)
          (affine_point_valid b x₂ y₂ rightOnCurve) := by
  apply WeierstrassCurve.Projective.Point.ext
  rw [WeierstrassCurve.Projective.Point.add_point]
  change
    (⟦projectiveCoordinates (affinePoint (shortSecantX x₁ y₁ x₂ y₂)
      (shortSecantY x₁ y₁ x₂ y₂))⟧ :
      WeierstrassCurve.Projective.PointClass K) =
      (shortWeierstrass b).toProjective.addMap
        ⟦projectiveCoordinates (affinePoint x₁ y₁)⟧
        ⟦projectiveCoordinates (affinePoint x₂ y₂)⟧
  rw [WeierstrassCurve.Projective.addMap_of_Z_ne_zero
    ((projective_equation_iff b _).mpr
      (affine_point_valid b x₁ y₁ leftOnCurve).1)
    ((projective_equation_iff b _).mpr
      (affine_point_valid b x₂ y₂ rightOnCurve).1)
    (by simp [projectiveCoordinates, affinePoint])
    (by simp [projectiveCoordinates, affinePoint])
    (by
      intro exceptional
      apply xDistinct
      simpa [projectiveCoordinates, affinePoint] using exceptional.1)]
  congr 1
  funext i
  fin_cases i <;>
    simp [projectiveCoordinates, affinePoint, shortWeierstrass,
      shortSecantX, shortSecantY,
      WeierstrassCurve.Affine.slope_of_X_ne xDistinct,
      WeierstrassCurve.Affine.addY, WeierstrassCurve.Affine.negY]

theorem projective_affine_tangent_eq_add [Field K] [DecidableEq K]
    (b x y : K) [elliptic : (shortWeierstrass b).IsElliptic]
    (twoNonzero : (2 : K) ≠ 0) (onCurve : y ^ 2 = x ^ 3 + b)
    (yNonzero : y ≠ 0) :
    toProjectivePoint b
        (affinePoint (shortTangentX x y) (shortTangentY x y))
        (affine_point_valid b _ _
          (short_tangent_on_curve b x y twoNonzero onCurve yNonzero)) =
      toProjectivePoint b (affinePoint x y)
          (affine_point_valid b x y onCurve) +
        toProjectivePoint b (affinePoint x y)
          (affine_point_valid b x y onCurve) := by
  apply WeierstrassCurve.Projective.Point.ext
  rw [WeierstrassCurve.Projective.Point.add_point]
  change
    (⟦projectiveCoordinates (affinePoint (shortTangentX x y)
      (shortTangentY x y))⟧ : WeierstrassCurve.Projective.PointClass K) =
      (shortWeierstrass b).toProjective.addMap
        ⟦projectiveCoordinates (affinePoint x y)⟧
        ⟦projectiveCoordinates (affinePoint x y)⟧
  rw [WeierstrassCurve.Projective.addMap_of_Z_ne_zero
    ((projective_equation_iff b _).mpr
      (affine_point_valid b x y onCurve).1)
    ((projective_equation_iff b _).mpr
      (affine_point_valid b x y onCurve).1)
    (by simp [projectiveCoordinates, affinePoint])
    (by simp [projectiveCoordinates, affinePoint])
    (by
      intro exceptional
      exact (y_ne_neg_of_two_ne_zero y twoNonzero yNonzero) (by
        simpa [projectiveCoordinates, affinePoint, shortWeierstrass,
          WeierstrassCurve.Projective.negY] using exceptional.2))]
  congr 1
  funext i
  fin_cases i <;>
    simp [projectiveCoordinates, affinePoint,
      shortWeierstrass, shortTangentX, shortTangentY,
      WeierstrassCurve.Affine.slope,
      y_ne_neg_of_two_ne_zero y twoNonzero yNonzero,
      WeierstrassCurve.Affine.addY, WeierstrassCurve.Affine.negY] <;> ring

theorem projective_public_double_affine_eq_add [Field K] [DecidableEq K]
    (b x y : K) [elliptic : (shortWeierstrass b).IsElliptic]
    (noTwoTorsion : NoAffineTwoTorsion b) (twoNonzero : (2 : K) ≠ 0)
    (onCurve : y ^ 2 = x ^ 3 + b) :
    toProjectivePoint b (cppPublicDouble b (affinePoint x y))
        (cpp_public_double_affine_valid b x y onCurve twoNonzero
          (affine_y_ne_zero_of_no_two_torsion b x y noTwoTorsion onCurve)) =
      toProjectivePoint b (affinePoint x y)
          (affine_point_valid b x y onCurve) +
        toProjectivePoint b (affinePoint x y)
          (affine_point_valid b x y onCurve) := by
  let yNonzero := affine_y_ne_zero_of_no_two_torsion b x y
    noTwoTorsion onCurve
  have doubleEquivalent :
      (cppPublicDouble b (affinePoint x y)).equivalent
        (affinePoint (shortTangentX x y) (shortTangentY x y)) := by
    rw [cpp_public_double_affine_eq_scaled_tangent b x y onCurve
      twoNonzero yNonzero]
    exact homogeneous_scale_equivalent
      (cppPublicDouble b (affinePoint x y)).z
      (affinePoint (shortTangentX x y) (shortTangentY x y))
      (cpp_public_double_affine_z_ne_zero b x y twoNonzero yNonzero)
  have doubleEqualsTangent := (to_projective_point_eq_iff b
    (cppPublicDouble b (affinePoint x y))
    (affinePoint (shortTangentX x y) (shortTangentY x y))
    (cpp_public_double_affine_valid b x y onCurve twoNonzero yNonzero)
    (affine_point_valid b _ _
      (short_tangent_on_curve b x y twoNonzero onCurve yNonzero))).mpr
      doubleEquivalent
  exact doubleEqualsTangent.trans
    (projective_affine_tangent_eq_add b x y twoNonzero onCurve yNonzero)

theorem projective_complete_add_affine_inverse_eq_add
    [Field K] [DecidableEq K]
    (b x y : K) [elliptic : (shortWeierstrass b).IsElliptic]
    (noTwoTorsion : NoAffineTwoTorsion b) (twoNonzero : (2 : K) ≠ 0)
    (onCurve : y ^ 2 = x ^ 3 + b) :
    toProjectivePoint b
        (cppCompleteAdd b (affinePoint x y) (affinePoint x (-y)))
        (cpp_complete_add_affine_inverse_valid b x y noTwoTorsion
          twoNonzero onCurve) =
      toProjectivePoint b (affinePoint x y)
          (affine_point_valid b x y onCurve) +
        toProjectivePoint b (affinePoint x (-y))
          (affine_point_valid b x (-y) (by simpa using onCurve)) := by
  have outputEquivalentZero :
      (cppCompleteAdd b (affinePoint x y)
        (affinePoint x (-y))).equivalent cppProjectiveZero := by
    simpa [affinePoint, HomogeneousPoint.negate] using
      cpp_complete_add_inverse b (affinePoint x y)
  have outputEqualsZero := (to_projective_point_eq_iff b
    (cppCompleteAdd b (affinePoint x y) (affinePoint x (-y)))
    cppProjectiveZero
    (cpp_complete_add_affine_inverse_valid b x y noTwoTorsion
      twoNonzero onCurve)
    (projective_zero_valid b)).mpr outputEquivalentZero
  have rightEqualsNeg :
      toProjectivePoint b (affinePoint x (-y))
          (affine_point_valid b x (-y) (by simpa using onCurve)) =
        -toProjectivePoint b (affinePoint x y)
          (affine_point_valid b x y onCurve) := by
    simpa [affinePoint, HomogeneousPoint.negate] using
      to_projective_point_negate b (affinePoint x y)
        (affine_point_valid b x y onCurve)
  rw [rightEqualsNeg, add_neg_cancel]
  exact outputEqualsZero.trans (to_projective_point_zero b)

theorem projective_complete_add_affine_eq_add [Field K] [DecidableEq K]
    (b x₁ y₁ x₂ y₂ : K) [elliptic : (shortWeierstrass b).IsElliptic]
    (noTwoTorsion : NoAffineTwoTorsion b)
    (leftOnCurve : y₁ ^ 2 = x₁ ^ 3 + b)
    (rightOnCurve : y₂ ^ 2 = x₂ ^ 3 + b) (xDistinct : x₁ ≠ x₂) :
    toProjectivePoint b
        (cppCompleteAdd b (affinePoint x₁ y₁) (affinePoint x₂ y₂))
        (cpp_complete_add_affine_valid b x₁ y₁ x₂ y₂ noTwoTorsion
          leftOnCurve rightOnCurve xDistinct) =
      toProjectivePoint b (affinePoint x₁ y₁)
          (affine_point_valid b x₁ y₁ leftOnCurve) +
        toProjectivePoint b (affinePoint x₂ y₂)
          (affine_point_valid b x₂ y₂ rightOnCurve) := by
  have outputEquivalent := cpp_complete_add_affine_equivalent_secant
    b x₁ y₁ x₂ y₂ noTwoTorsion leftOnCurve rightOnCurve xDistinct
  have outputEqualsSecant := (to_projective_point_eq_iff b
    (cppCompleteAdd b (affinePoint x₁ y₁) (affinePoint x₂ y₂))
    (affinePoint (shortSecantX x₁ y₁ x₂ y₂)
      (shortSecantY x₁ y₁ x₂ y₂))
    (cpp_complete_add_affine_valid b x₁ y₁ x₂ y₂ noTwoTorsion
      leftOnCurve rightOnCurve xDistinct)
    (affine_point_valid b _ _
      (short_secant_on_curve b x₁ y₁ x₂ y₂
        leftOnCurve rightOnCurve xDistinct))).mpr outputEquivalent
  exact outputEqualsSecant.trans
    (projective_affine_secant_eq_add b x₁ y₁ x₂ y₂
      leftOnCurve rightOnCurve xDistinct)

theorem affine_y_eq_or_eq_neg [Field K] (b x₁ y₁ x₂ y₂ : K)
    (leftOnCurve : y₁ ^ 2 = x₁ ^ 3 + b)
    (rightOnCurve : y₂ ^ 2 = x₂ ^ 3 + b) (xEqual : x₁ = x₂) :
    y₁ = y₂ ∨ y₂ = -y₁ := by
  have squareEqual : y₁ ^ 2 = y₂ ^ 2 := by
    calc
      y₁ ^ 2 = x₁ ^ 3 + b := leftOnCurve
      _ = x₂ ^ 3 + b := by rw [xEqual]
      _ = y₂ ^ 2 := rightOnCurve.symm
  have productZero : (y₁ - y₂) * (y₁ + y₂) = 0 := by
    linear_combination squareEqual
  rcases mul_eq_zero.mp productZero with differenceZero | sumZero
  · left
    exact sub_eq_zero.mp differenceZero
  · right
    linear_combination sumZero

theorem cpp_complete_add_affine_valid_all [Field K] [DecidableEq K]
    (b x₁ y₁ x₂ y₂ : K) (noTwoTorsion : NoAffineTwoTorsion b)
    (twoNonzero : (2 : K) ≠ 0)
    (leftOnCurve : y₁ ^ 2 = x₁ ^ 3 + b)
    (rightOnCurve : y₂ ^ 2 = x₂ ^ 3 + b) :
    (cppCompleteAdd b (affinePoint x₁ y₁)
      (affinePoint x₂ y₂)).projectivelyValid b := by
  by_cases xDistinct : x₁ ≠ x₂
  · exact cpp_complete_add_affine_valid b x₁ y₁ x₂ y₂ noTwoTorsion
      leftOnCurve rightOnCurve xDistinct
  · have xEqual : x₁ = x₂ := not_ne_iff.mp xDistinct
    rcases affine_y_eq_or_eq_neg b x₁ y₁ x₂ y₂
      leftOnCurve rightOnCurve xEqual with yEqual | yNeg
    · subst x₂
      subst y₂
      rw [cpp_complete_add_self_eq_double b (affinePoint x₁ y₁)
        (affine_point_valid b x₁ y₁ leftOnCurve).1]
      exact cpp_public_double_affine_valid b x₁ y₁ leftOnCurve twoNonzero
        (affine_y_ne_zero_of_no_two_torsion b x₁ y₁
          noTwoTorsion leftOnCurve)
    · subst x₂
      rw [yNeg]
      exact cpp_complete_add_affine_inverse_valid b x₁ y₁ noTwoTorsion
        twoNonzero leftOnCurve

theorem cpp_complete_add_projectively_valid [Field K] [DecidableEq K]
    (b : K) (noTwoTorsion : NoAffineTwoTorsion b)
    (twoNonzero : (2 : K) ≠ 0)
    (left right : HomogeneousPoint K)
    (leftValid : left.projectivelyValid b)
    (rightValid : right.projectivelyValid b) :
    (cppCompleteAdd b left right).projectivelyValid b := by
  by_cases leftAtInfinity : left.z = 0
  · rw [cpp_complete_add_infinity_left_eq_scale b left right
      leftValid leftAtInfinity]
    apply homogeneous_scale_projectively_valid b
    · exact mul_ne_zero
        (pow_ne_zero 2
          (projectively_valid_y_ne_zero b noTwoTorsion left leftValid))
        (projectively_valid_y_ne_zero b noTwoTorsion right rightValid)
    · exact rightValid
  · by_cases rightAtInfinity : right.z = 0
    · rw [cpp_complete_add_infinity_right_eq_scale b left right
        rightValid rightAtInfinity]
      apply homogeneous_scale_projectively_valid b
      · exact mul_ne_zero
          (pow_ne_zero 2
            (projectively_valid_y_ne_zero b noTwoTorsion right rightValid))
          (projectively_valid_y_ne_zero b noTwoTorsion left leftValid)
      · exact leftValid
    · rw [← homogeneous_scale_affine_representative left leftAtInfinity,
        ← homogeneous_scale_affine_representative right rightAtInfinity,
        cpp_complete_add_homogeneous_scale]
      apply homogeneous_scale_projectively_valid b
      · exact pow_ne_zero 2 (mul_ne_zero leftAtInfinity rightAtInfinity)
      · exact cpp_complete_add_affine_valid_all b
          (affineRepresentative left).x (affineRepresentative left).y
          (affineRepresentative right).x (affineRepresentative right).y
          noTwoTorsion twoNonzero
          (affine_representative_equation b left leftValid leftAtInfinity)
          (affine_representative_equation b right rightValid rightAtInfinity)

theorem projective_complete_add_affine_eq_add_all
    [Field K] [DecidableEq K]
    (b x₁ y₁ x₂ y₂ : K) [elliptic : (shortWeierstrass b).IsElliptic]
    (noTwoTorsion : NoAffineTwoTorsion b) (twoNonzero : (2 : K) ≠ 0)
    (leftOnCurve : y₁ ^ 2 = x₁ ^ 3 + b)
    (rightOnCurve : y₂ ^ 2 = x₂ ^ 3 + b) :
    toProjectivePoint b
        (cppCompleteAdd b (affinePoint x₁ y₁) (affinePoint x₂ y₂))
        (cpp_complete_add_affine_valid_all b x₁ y₁ x₂ y₂ noTwoTorsion
          twoNonzero leftOnCurve rightOnCurve) =
      toProjectivePoint b (affinePoint x₁ y₁)
          (affine_point_valid b x₁ y₁ leftOnCurve) +
        toProjectivePoint b (affinePoint x₂ y₂)
          (affine_point_valid b x₂ y₂ rightOnCurve) := by
  by_cases xDistinct : x₁ ≠ x₂
  · exact projective_complete_add_affine_eq_add b x₁ y₁ x₂ y₂
      noTwoTorsion leftOnCurve rightOnCurve xDistinct
  · have xEqual : x₁ = x₂ := not_ne_iff.mp xDistinct
    rcases affine_y_eq_or_eq_neg b x₁ y₁ x₂ y₂
      leftOnCurve rightOnCurve xEqual with yEqual | yNeg
    · subst x₂
      subst y₂
      simpa only [cpp_complete_add_self_eq_double b (affinePoint x₁ y₁)
        (affine_point_valid b x₁ y₁ leftOnCurve).1] using
          projective_public_double_affine_eq_add b x₁ y₁
            noTwoTorsion twoNonzero leftOnCurve
    · subst x₂
      simpa only [yNeg] using
        projective_complete_add_affine_inverse_eq_add b x₁ y₁
          noTwoTorsion twoNonzero leftOnCurve

theorem projective_complete_add_eq_add [Field K] [DecidableEq K]
    (b : K) [elliptic : (shortWeierstrass b).IsElliptic]
    (noTwoTorsion : NoAffineTwoTorsion b) (twoNonzero : (2 : K) ≠ 0)
    (left right : HomogeneousPoint K)
    (leftValid : left.projectivelyValid b)
    (rightValid : right.projectivelyValid b) :
    toProjectivePoint b (cppCompleteAdd b left right)
        (cpp_complete_add_projectively_valid b noTwoTorsion twoNonzero
          left right leftValid rightValid) =
      toProjectivePoint b left leftValid +
        toProjectivePoint b right rightValid := by
  by_cases leftAtInfinity : left.z = 0
  · have leftIsZero := to_projective_point_eq_zero_of_z_eq_zero b
      left leftValid leftAtInfinity
    have leftYNonzero := projectively_valid_y_ne_zero b noTwoTorsion
      left leftValid
    have rightYNonzero := projectively_valid_y_ne_zero b noTwoTorsion
      right rightValid
    have factorNonzero : left.y ^ 2 * right.y ≠ 0 :=
      mul_ne_zero (pow_ne_zero 2 leftYNonzero) rightYNonzero
    have outputScale := to_projective_point_homogeneous_scale b
      (left.y ^ 2 * right.y) right factorNonzero rightValid
    have outputIsRight :
        toProjectivePoint b (cppCompleteAdd b left right)
            (cpp_complete_add_projectively_valid b noTwoTorsion twoNonzero
              left right leftValid rightValid) =
          toProjectivePoint b right rightValid := by
      simpa only [cpp_complete_add_infinity_left_eq_scale b left right
        leftValid leftAtInfinity] using outputScale
    rw [leftIsZero, zero_add]
    exact outputIsRight
  · by_cases rightAtInfinity : right.z = 0
    · have rightIsZero := to_projective_point_eq_zero_of_z_eq_zero b
        right rightValid rightAtInfinity
      have rightYNonzero := projectively_valid_y_ne_zero b noTwoTorsion
        right rightValid
      have leftYNonzero := projectively_valid_y_ne_zero b noTwoTorsion
        left leftValid
      have factorNonzero : right.y ^ 2 * left.y ≠ 0 :=
        mul_ne_zero (pow_ne_zero 2 rightYNonzero) leftYNonzero
      have outputScale := to_projective_point_homogeneous_scale b
        (right.y ^ 2 * left.y) left factorNonzero leftValid
      have outputIsLeft :
          toProjectivePoint b (cppCompleteAdd b left right)
              (cpp_complete_add_projectively_valid b noTwoTorsion twoNonzero
                left right leftValid rightValid) =
            toProjectivePoint b left leftValid := by
        simpa only [cpp_complete_add_infinity_right_eq_scale b left right
          rightValid rightAtInfinity] using outputScale
      rw [rightIsZero, add_zero]
      exact outputIsLeft
    · let leftAffine := affineRepresentative left
      let rightAffine := affineRepresentative right
      have leftAffineValid := affine_representative_valid b left
        leftValid leftAtInfinity
      have rightAffineValid := affine_representative_valid b right
        rightValid rightAtInfinity
      have leftEquation := affine_representative_equation b left
        leftValid leftAtInfinity
      have rightEquation := affine_representative_equation b right
        rightValid rightAtInfinity
      have affineOutputValid := cpp_complete_add_affine_valid_all b
        leftAffine.x leftAffine.y rightAffine.x rightAffine.y
        noTwoTorsion twoNonzero leftEquation rightEquation
      have scaleFactorNonzero : (left.z * right.z) ^ 2 ≠ 0 :=
        pow_ne_zero 2 (mul_ne_zero leftAtInfinity rightAtInfinity)
      have leftRepresentation : homogeneousScale left.z leftAffine = left := by
        simpa [leftAffine] using
          homogeneous_scale_affine_representative left leftAtInfinity
      have rightRepresentation : homogeneousScale right.z rightAffine = right := by
        simpa [rightAffine] using
          homogeneous_scale_affine_representative right rightAtInfinity
      have outputEquality :
          cppCompleteAdd b left right =
            homogeneousScale ((left.z * right.z) ^ 2)
              (cppCompleteAdd b leftAffine rightAffine) := by
        calc
          cppCompleteAdd b left right =
              cppCompleteAdd b (homogeneousScale left.z leftAffine)
                (homogeneousScale right.z rightAffine) := by
            rw [leftRepresentation, rightRepresentation]
          _ = homogeneousScale ((left.z * right.z) ^ 2)
              (cppCompleteAdd b leftAffine rightAffine) :=
            cpp_complete_add_homogeneous_scale b left.z right.z
              leftAffine rightAffine
      have outputEquivalent :
          (cppCompleteAdd b left right).equivalent
            (cppCompleteAdd b leftAffine rightAffine) := by
        rw [outputEquality]
        exact homogeneous_scale_equivalent ((left.z * right.z) ^ 2)
          (cppCompleteAdd b leftAffine rightAffine) scaleFactorNonzero
      have outputEqualsAffine := (to_projective_point_eq_iff b
        (cppCompleteAdd b left right) (cppCompleteAdd b leftAffine rightAffine)
        (cpp_complete_add_projectively_valid b noTwoTorsion twoNonzero
          left right leftValid rightValid)
        (by simpa [leftAffine, rightAffine, affineRepresentative, affinePoint]
          using affineOutputValid)).mpr outputEquivalent
      have leftEquivalent : left.equivalent leftAffine := by
        rw [← homogeneous_scale_affine_representative left leftAtInfinity]
        exact homogeneous_scale_equivalent left.z leftAffine leftAtInfinity
      have rightEquivalent : right.equivalent rightAffine := by
        rw [← homogeneous_scale_affine_representative right rightAtInfinity]
        exact homogeneous_scale_equivalent right.z rightAffine rightAtInfinity
      have leftEqualsAffine := (to_projective_point_eq_iff b left leftAffine
        leftValid leftAffineValid).mpr leftEquivalent
      have rightEqualsAffine := (to_projective_point_eq_iff b right rightAffine
        rightValid rightAffineValid).mpr rightEquivalent
      rw [outputEqualsAffine, leftEqualsAffine, rightEqualsAffine]
      simpa [leftAffine, rightAffine, affineRepresentative, affinePoint] using
        projective_complete_add_affine_eq_add_all b
          leftAffine.x leftAffine.y rightAffine.x rightAffine.y
          noTwoTorsion twoNonzero leftEquation rightEquation

/-! ## Selected BLS12-381 fields -/

set_option maxRecDepth 100000 in
theorem selected_g1_no_affine_two_torsion :
    NoAffineTwoTorsion (4 : SelectedCurveFp) := by
  intro x root
  have xNonzero : x ≠ 0 := by
    intro xZero
    rw [xZero] at root
    have fourNonzero : (4 : SelectedCurveFp) ≠ 0 := by native_decide
    exact fourNonzero (by simpa using root)
  have cube : x ^ 3 = -(4 : SelectedCurveFp) := by
    linear_combination root
  let exponent := (Generated.target_modulus - 1) / 3
  have exponentIdentity : 3 * exponent = Generated.target_modulus - 1 := by
    native_decide
  have powered := congrArg (fun value : SelectedCurveFp => value ^ exponent) cube
  have leftOne : (x ^ 3) ^ exponent = 1 := by
    rw [← pow_mul, exponentIdentity]
    exact ZMod.pow_card_sub_one_eq_one xNonzero
  have rightNotOne : (-(4 : SelectedCurveFp)) ^ exponent ≠ 1 := by
    native_decide
  exact rightNotOne (powered.symm.trans leftOne)

def fp2ValueEquivProd (R : Type*) : Fp2Value R ≃ R × R where
  toFun value := (value.x, value.y)
  invFun value := ⟨value.1, value.2⟩
  left_inv value := by cases value; rfl
  right_inv value := by cases value; rfl

noncomputable instance fp2ValueFintype {R : Type*} [Fintype R] :
    Fintype (Fp2Value R) :=
  Fintype.ofEquiv (R × R) (fp2ValueEquivProd R).symm

theorem selected_fp2_card :
    Fintype.card SelectedMillerFp2 = Generated.target_modulus ^ 2 := by
  rw [Fintype.card_congr (fp2ValueEquivProd SelectedCurveFp)]
  simp [ZMod.card, pow_two]

/-- A proof-only computable view of exponentiation. The selected field instance
has a noncomputable inverse, although exponentiation itself uses only the
already-computable `Fp2Value` commutative monoid. -/
def computableSelectedFp2Pow (value : SelectedMillerFp2) (exponent : Nat) :
    SelectedMillerFp2 := by
  letI : CommMonoid SelectedMillerFp2 := Fp2Concrete.fp2CommMonoid
  exact npowRecAuto exponent value

theorem computable_selected_fp2_pow_eq (value : SelectedMillerFp2)
    (exponent : Nat) :
    computableSelectedFp2Pow value exponent = value ^ exponent := by
  rfl

set_option maxRecDepth 100000 in
theorem selected_g2_no_affine_two_torsion :
    NoAffineTwoTorsion selectedG2CurveB := by
  intro x root
  have xNonzero : x ≠ 0 := by
    intro xZero
    rw [xZero] at root
    have bNonzero : selectedG2CurveB ≠ 0 := by native_decide
    exact bNonzero (by simpa using root)
  have cube : x ^ 3 = -selectedG2CurveB := by
    linear_combination root
  let fieldCard := Generated.target_modulus ^ 2
  let exponent := (fieldCard - 1) / 3
  have exponentIdentity : 3 * exponent = fieldCard - 1 := by native_decide
  have powered := congrArg
    (fun value : SelectedMillerFp2 => value ^ exponent) cube
  have leftOne : (x ^ 3) ^ exponent = 1 := by
    rw [← pow_mul, exponentIdentity]
    dsimp [fieldCard]
    rw [← selected_fp2_card]
    exact FiniteField.pow_card_sub_one_eq_one x xNonzero
  have computableRightNotOne :
      computableSelectedFp2Pow (-selectedG2CurveB) exponent ≠ 1 := by
    native_decide
  have rightNotOne : (-selectedG2CurveB) ^ exponent ≠ 1 := by
    rwa [computable_selected_fp2_pow_eq] at computableRightNotOne
  exact rightNotOne (powered.symm.trans leftOne)

theorem selected_g1_two_ne_zero : (2 : SelectedCurveFp) ≠ 0 := by
  native_decide

theorem selected_g2_two_ne_zero : (2 : SelectedMillerFp2) ≠ 0 := by
  intro twoZero
  have realZero := congrArg Fp2Value.x twoZero
  change (2 : SelectedCurveFp) = 0 at realZero
  exact selected_g1_two_ne_zero realZero

theorem selected_g1_complete_add_eq_add
    (left right : HomogeneousPoint SelectedCurveFp)
    (leftValid : left.projectivelyValid 4)
    (rightValid : right.projectivelyValid 4) :
    toSelectedG1Point (cppCompleteAdd 4 left right)
        (cpp_complete_add_projectively_valid 4
          selected_g1_no_affine_two_torsion selected_g1_two_ne_zero
          left right leftValid rightValid) =
      toSelectedG1Point left leftValid + toSelectedG1Point right rightValid := by
  exact projective_complete_add_eq_add 4
    selected_g1_no_affine_two_torsion selected_g1_two_ne_zero
    left right leftValid rightValid

theorem selected_g2_complete_add_eq_add
    (left right : HomogeneousPoint SelectedMillerFp2)
    (leftValid : left.projectivelyValid selectedG2CurveB)
    (rightValid : right.projectivelyValid selectedG2CurveB) :
    toSelectedG2Point (cppCompleteAdd selectedG2CurveB left right)
        (cpp_complete_add_projectively_valid selectedG2CurveB
          selected_g2_no_affine_two_torsion selected_g2_two_ne_zero
          left right leftValid rightValid) =
      toSelectedG2Point left leftValid + toSelectedG2Point right rightValid := by
  exact projective_complete_add_eq_add selectedG2CurveB
    selected_g2_no_affine_two_torsion selected_g2_two_ne_zero
    left right leftValid rightValid

end CatBLST.AVX.MatrixNoK
