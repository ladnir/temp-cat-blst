/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Fp2Algebra
import Mathlib.Tactic

set_option maxHeartbeats 800000

/-!
# Source-shaped Miller doubling geometry

This is the first P4 packet. It models the exact homogeneous point and line
formulas shared by `Miller::double_step` and the prepared `line_dbl` path.
The output line is the tangent at the input and contains the negation of the
doubled point, as required by the Miller divisor. Sparse Fp12 multiplication
and the fixed Miller schedule remain separate correspondence obligations.
-/

namespace CatBLST.AVX.MatrixNoK

structure HomogeneousPoint (K : Type*) where
  x : K
  y : K
  z : K
  deriving DecidableEq

@[ext] theorem HomogeneousPoint.ext {left right : HomogeneousPoint K}
    (x : left.x = right.x) (y : left.y = right.y) (z : left.z = right.z) :
    left = right := by
  cases left
  cases right
  simp_all

structure MillerLine (K : Type*) where
  c : K
  x : K
  y : K
  deriving DecidableEq

@[ext] theorem MillerLine.ext {left right : MillerLine K}
    (c : left.c = right.c) (x : left.x = right.x) (y : left.y = right.y) :
    left = right := by
  cases left
  cases right
  simp_all

structure MillerDoubleResult (K : Type*) where
  point : HomogeneousPoint K
  line : MillerLine K
  deriving DecidableEq

def HomogeneousPoint.onCurve [CommRing K]
    (b : K) (point : HomogeneousPoint K) : Prop :=
  point.y ^ 2 * point.z = point.x ^ 3 + b * point.z ^ 3

def HomogeneousPoint.negate [Neg K]
    (point : HomogeneousPoint K) : HomogeneousPoint K :=
  ⟨point.x, -point.y, point.z⟩

def MillerLine.evaluateHomogeneous [CommRing K]
    (line : MillerLine K) (point : HomogeneousPoint K) : K :=
  line.c * point.z + line.x * point.x + line.y * point.y

/-- Prepared-line scaling performed later by `scale_line_xy`. -/
def MillerLine.scaleForG1 [CommRing K]
    (line : MillerLine K) (g1x g1y : K) : MillerLine K :=
  ⟨line.c, line.x * g1x, line.y * g1y⟩

/--
Literal algebra after the accepted Fp2 reductions in both C++ doubling paths.
The names and grouping retain the source's `yy`, `zz`, `xx`, `xy`, `yz`,
`b3_zz`, `yy_minus_b9_zz`, and second multiplication round.
-/
def cppMillerDouble [CommRing K]
    (b : K) (point : HomogeneousPoint K) : MillerDoubleResult K :=
  let yy := point.y ^ 2
  let zz := point.z ^ 2
  let xx := point.x ^ 2
  let xy := point.x * point.y
  let yz := point.y * point.z
  let eightYY := 8 * yy
  let b3ZZ := 3 * b * zz
  let yyPlusB3ZZ := yy + b3ZZ
  let b9ZZ := 3 * b3ZZ
  let yyMinusB9ZZ := yy - b9ZZ
  let twoXY := 2 * xy
  let pointOut : HomogeneousPoint K :=
    ⟨yyMinusB9ZZ * twoXY,
      yyMinusB9ZZ * yyPlusB3ZZ + eightYY * b3ZZ,
      eightYY * yz⟩
  let lineOut : MillerLine K :=
    ⟨b3ZZ - yy, 3 * xx, -(2 * yz)⟩
  ⟨pointOut, lineOut⟩

/-- Line coefficients formed directly inside `Miller::double_step`. -/
def cppMillerDoubleDirectLine [CommRing K]
    (b : K) (point : HomogeneousPoint K) (g1x g1y : K) : MillerLine K :=
  let yy := point.y ^ 2
  let zz := point.z ^ 2
  let xx := point.x ^ 2
  let yz := point.y * point.z
  ⟨3 * b * zz - yy, xx * (3 * g1x), yz * (-2 * g1y)⟩

theorem cpp_miller_double_point_coordinates [CommRing K]
    (b : K) (point : HomogeneousPoint K) :
    (cppMillerDouble b point).point =
      ⟨(point.y ^ 2 - 9 * b * point.z ^ 2) *
          (2 * point.x * point.y),
        (point.y ^ 2 - 9 * b * point.z ^ 2) *
            (point.y ^ 2 + 3 * b * point.z ^ 2) +
          24 * b * point.y ^ 2 * point.z ^ 2,
        8 * point.y ^ 3 * point.z⟩ := by
  ext
  all_goals simp only [cppMillerDouble]
  all_goals ring

theorem cpp_miller_double_line_coordinates [CommRing K]
    (b : K) (point : HomogeneousPoint K) :
    (cppMillerDouble b point).line =
      ⟨3 * b * point.z ^ 2 - point.y ^ 2,
        3 * point.x ^ 2,
        -(2 * point.y * point.z)⟩ := by
  ext
  all_goals simp only [cppMillerDouble]
  all_goals ring

theorem cpp_miller_double_prepared_line_matches_direct [CommRing K]
    (b : K) (point : HomogeneousPoint K) (g1x g1y : K) :
    (cppMillerDouble b point).line.scaleForG1 g1x g1y =
      cppMillerDoubleDirectLine b point g1x g1y := by
  rw [cpp_miller_double_line_coordinates]
  ext
  all_goals simp only [MillerLine.scaleForG1, cppMillerDoubleDirectLine]
  all_goals ring

theorem cpp_miller_double_line_contains_input [CommRing K]
    (b : K) (point : HomogeneousPoint K)
    (onCurve : point.onCurve b) :
    (cppMillerDouble b point).line.evaluateHomogeneous point = 0 := by
  simp only [HomogeneousPoint.onCurve] at onCurve
  simp only [MillerLine.evaluateHomogeneous]
  rw [cpp_miller_double_line_coordinates]
  simp only
  linear_combination -3 * onCurve

theorem cpp_miller_double_line_contains_negative_output [CommRing K]
    (b : K) (point : HomogeneousPoint K)
    (onCurve : point.onCurve b) :
    (cppMillerDouble b point).line.evaluateHomogeneous
        (cppMillerDouble b point).point.negate = 0 := by
  simp only [HomogeneousPoint.onCurve] at onCurve
  rw [cpp_miller_double_line_coordinates,
    cpp_miller_double_point_coordinates]
  simp only [MillerLine.evaluateHomogeneous, HomogeneousPoint.negate]
  linear_combination
    -6 * point.y * (point.y ^ 2 - 9 * b * point.z ^ 2) * onCurve

theorem cpp_miller_double_preserves_curve [CommRing K]
    (b : K) (point : HomogeneousPoint K)
    (onCurve : point.onCurve b) :
    (cppMillerDouble b point).point.onCurve b := by
  simp only [HomogeneousPoint.onCurve] at onCurve
  rw [cpp_miller_double_point_coordinates]
  simp only [HomogeneousPoint.onCurve]
  linear_combination
    8 * point.y ^ 3 * (point.y ^ 2 - 9 * b * point.z ^ 2) ^ 3 * onCurve

/-- The exact selected P3 Fp2 carrier used by the MatrixNoK G2/Miller path. -/
abbrev SelectedMillerFp2 := Fp2Value (ZMod Generated.target_modulus)

theorem selected_cpp_miller_double_correct
    (b : SelectedMillerFp2) (point : HomogeneousPoint SelectedMillerFp2)
    (onCurve : point.onCurve b) :
    (cppMillerDouble b point).point.onCurve b ∧
      (cppMillerDouble b point).line.evaluateHomogeneous point = 0 ∧
      (cppMillerDouble b point).line.evaluateHomogeneous
          (cppMillerDouble b point).point.negate = 0 := by
  exact ⟨cpp_miller_double_preserves_curve b point onCurve,
    cpp_miller_double_line_contains_input b point onCurve,
    cpp_miller_double_line_contains_negative_output b point onCurve⟩

end CatBLST.AVX.MatrixNoK
