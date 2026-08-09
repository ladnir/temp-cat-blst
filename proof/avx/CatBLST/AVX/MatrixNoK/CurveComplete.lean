/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.MillerDouble
import Mathlib.Tactic

set_option maxHeartbeats 1600000

/-!
# Complete homogeneous G1/G2 operations

The MatrixNoK G1 and G2 owners share the same complete-addition schedule in
`matrixnok_curve.h` and the same Algorithm-9 doubling schedule in `ec.h`.
This file transcribes the reduced field formulas exactly.  The generic
commutative-ring theorems apply to both the selected Fp and Fp2 carriers.
-/

namespace CatBLST.AVX.MatrixNoK

/-- Exact reduced-field meaning of `matrixnok::detail::add_point`. -/
def cppCompleteAdd [CommRing K] (b : K)
    (left right : HomogeneousPoint K) : HomogeneousPoint K :=
  let xx := left.x * right.x
  let yy := left.y * right.y
  let zz := left.z * right.z
  let xyCross := (left.x + left.y) * (right.x + right.y) - xx - yy
  let yzCross := (left.y + left.z) * (right.y + right.z) - yy - zz
  let xzCross := (left.x + left.z) * (right.x + right.z) - xx - zz
  let b3ZZ := 3 * b * zz
  let b3XZCross := 3 * b * xzCross
  let yyMinusB3ZZ := yy - b3ZZ
  let threeXX := 3 * xx
  let yyPlusB3ZZ := yy + b3ZZ
  ⟨yyMinusB3ZZ * xyCross - yzCross * b3XZCross,
    yyMinusB3ZZ * yyPlusB3ZZ + threeXX * b3XZCross,
    yzCross * yyPlusB3ZZ + threeXX * xyCross⟩

theorem cpp_complete_add_coordinates [CommRing K] (b : K)
    (left right : HomogeneousPoint K) :
    cppCompleteAdd b left right =
      ⟨(left.y * right.y - 3 * b * left.z * right.z) *
            (left.x * right.y + left.y * right.x) -
          (left.y * right.z + left.z * right.y) *
            (3 * b * (left.x * right.z + left.z * right.x)),
        (left.y * right.y - 3 * b * left.z * right.z) *
            (left.y * right.y + 3 * b * left.z * right.z) +
          (3 * left.x * right.x) *
            (3 * b * (left.x * right.z + left.z * right.x)),
        (left.y * right.z + left.z * right.y) *
            (left.y * right.y + 3 * b * left.z * right.z) +
          (3 * left.x * right.x) *
            (left.x * right.y + left.y * right.x)⟩ := by
  ext <;> simp only [cppCompleteAdd] <;> ring

/-!
Offline polynomial certificates for curve preservation.  They were obtained
by reducing `Y²Z - X³ - bZ³` by the two input curve equations; Lean's
`linear_combination` rechecks the complete identity below.
-/
def completeAddLeftCertificate [CommRing K] (b lx ly lz rx ry rz : K) : K :=
  -27 * b ^ 3 * rz ^ 2 *
      (9 * lx * ly * lz * rx ^ 2 * rz ^ 2 -
        9 * lx * lz ^ 2 * rx ^ 2 * ry * rz + ly ^ 3 * rz ^ 4 +
        6 * ly ^ 2 * lz * ry * rz ^ 3 +
        12 * ly * lz ^ 2 * ry ^ 2 * rz ^ 2 -
        9 * lz ^ 3 * rx ^ 3 * ry + 8 * lz ^ 3 * ry ^ 3 * rz) -
    27 * b ^ 2 *
      (9 * lx ^ 3 * rx ^ 3 * ry * rz ^ 2 +
        9 * lx ^ 2 * ly * rx ^ 4 * rz ^ 2 +
        18 * lx ^ 2 * lz * rx ^ 4 * ry * rz +
        3 * lx * ly ^ 2 * rx ^ 2 * ry * rz ^ 3 +
        18 * lx * ly * lz * rx ^ 5 * rz +
        9 * lx * lz ^ 2 * rx ^ 5 * ry -
        3 * lx * lz ^ 2 * rx ^ 2 * ry ^ 3 * rz -
        ly ^ 3 * ry ^ 2 * rz ^ 4 +
        9 * ly ^ 2 * lz * rx ^ 3 * ry * rz ^ 2 -
        4 * ly ^ 2 * lz * ry ^ 3 * rz ^ 3 +
        9 * ly * lz ^ 2 * rx ^ 6 -
        4 * ly * lz ^ 2 * ry ^ 4 * rz ^ 2 +
        lz ^ 3 * rx ^ 3 * ry ^ 3) +
    9 * b *
      (3 * lx ^ 3 * rx ^ 3 * ry ^ 3 +
        9 * lx ^ 2 * ly * rx ^ 4 * ry ^ 2 +
        9 * lx * ly ^ 2 * rx ^ 5 * ry -
        3 * lx * ly ^ 2 * rx ^ 2 * ry ^ 3 * rz +
        3 * lx * ly * lz * rx ^ 2 * ry ^ 4 +
        3 * ly ^ 3 * rx ^ 6 - ly ^ 3 * ry ^ 4 * rz ^ 2 +
        3 * ly ^ 2 * lz * rx ^ 3 * ry ^ 3 -
        2 * ly ^ 2 * lz * ry ^ 5 * rz) +
    ly ^ 3 * ry ^ 6

def completeAddRightCertificate [CommRing K] (b lx ly lz rx ry rz : K) : K :=
  -243 * b ^ 4 * lz ^ 5 * rz ^ 2 * (ly * rz + lz * ry) +
    27 * b ^ 3 * lz ^ 3 *
      (18 * lx ^ 2 * ly * rx * rz ^ 2 +
        9 * lx ^ 2 * lz * rx * ry * rz +
        18 * lx * ly * lz * rx ^ 2 * rz +
        9 * lx * lz ^ 2 * rx ^ 2 * ry + 2 * ly ^ 3 * rz ^ 3 +
        6 * ly ^ 2 * lz * ry * rz ^ 2 +
        9 * ly * lz ^ 2 * rx ^ 3 + 3 * ly * lz ^ 2 * ry ^ 2 * rz) -
    27 * b ^ 2 * ly * lz *
      (6 * lx ^ 2 * ly ^ 2 * rx * rz ^ 2 +
        18 * lx ^ 2 * ly * lz * rx * ry * rz +
        6 * lx ^ 2 * lz ^ 2 * rx * ry ^ 2 +
        18 * lx * ly ^ 2 * lz * rx ^ 2 * rz +
        12 * lx * ly * lz ^ 2 * rx ^ 2 * ry + ly ^ 4 * rz ^ 3 +
        5 * ly ^ 3 * lz * ry * rz ^ 2 +
        10 * ly ^ 2 * lz ^ 2 * rx ^ 3 +
        6 * ly ^ 2 * lz ^ 2 * ry ^ 2 * rz + ly * lz ^ 3 * ry ^ 3) +
    9 * b * ly ^ 3 *
      (3 * lx ^ 2 * ly * rx * ry * rz +
        6 * lx ^ 2 * lz * rx * ry ^ 2 +
        9 * lx * ly * lz * rx ^ 2 * ry + 3 * ly ^ 2 * lz * rx ^ 3 +
        ly ^ 2 * lz * ry ^ 2 * rz + 2 * ly * lz ^ 2 * ry ^ 3) +
    ly ^ 6 * ry ^ 3

theorem cpp_complete_add_preserves_curve [CommRing K] (b : K)
    (left right : HomogeneousPoint K)
    (leftOnCurve : left.onCurve b) (rightOnCurve : right.onCurve b) :
    (cppCompleteAdd b left right).onCurve b := by
  simp only [HomogeneousPoint.onCurve] at leftOnCurve rightOnCurve ⊢
  rw [cpp_complete_add_coordinates]
  simp only
  have leftMultiple := congrArg
    (fun value => completeAddLeftCertificate b left.x left.y left.z
      right.x right.y right.z * value) leftOnCurve
  have rightMultiple := congrArg
    (fun value => completeAddRightCertificate b left.x left.y left.z
      right.x right.y right.z * value) rightOnCurve
  simp only [completeAddLeftCertificate] at leftMultiple
  simp only [completeAddRightCertificate] at rightMultiple
  linear_combination leftMultiple + rightMultiple

/-!
The selected public `double_point` is byte-for-byte the same reduced formula
already used by the proved Miller doubling point.  This theorem closes that
source-level connection rather than maintaining a second algebraic model.
-/
def cppPublicDouble [CommRing K]
    (b : K) (point : HomogeneousPoint K) : HomogeneousPoint K :=
  let yy := point.y ^ 2
  let zz := point.z ^ 2
  let xy := point.x * point.y
  let yz := point.y * point.z
  let eightYY := 8 * yy
  let b3ZZ := 3 * b * zz
  let yyPlusB3ZZ := yy + b3ZZ
  let yyMinusB9ZZ := yy - 3 * b3ZZ
  ⟨yyMinusB9ZZ * (2 * xy),
    yyMinusB9ZZ * yyPlusB3ZZ + eightYY * b3ZZ,
    eightYY * yz⟩

theorem cpp_public_double_eq_miller_point [CommRing K]
    (b : K) (point : HomogeneousPoint K) :
    cppPublicDouble b point = (cppMillerDouble b point).point := by
  rfl

theorem cpp_public_double_preserves_curve [CommRing K]
    (b : K) (point : HomogeneousPoint K) (onCurve : point.onCurve b) :
    (cppPublicDouble b point).onCurve b := by
  rw [cpp_public_double_eq_miller_point]
  exact cpp_miller_double_preserves_curve b point onCurve

def cppProjectiveZero [Zero K] [One K] : HomogeneousPoint K := ⟨0, 1, 0⟩

theorem cpp_projective_zero_on_curve [CommRing K] (b : K) :
    (cppProjectiveZero : HomogeneousPoint K).onCurve b := by
  simp [cppProjectiveZero, HomogeneousPoint.onCurve]

theorem cpp_negate_preserves_curve [CommRing K] (b : K)
    (point : HomogeneousPoint K) (onCurve : point.onCurve b) :
    point.negate.onCurve b := by
  simpa [HomogeneousPoint.negate, HomogeneousPoint.onCurve] using onCurve

/-- Equality of the affine values represented by homogeneous coordinates. -/
def HomogeneousPoint.equivalent [CommRing K]
    (left right : HomogeneousPoint K) : Prop :=
  (left.z = 0 ∧ right.z = 0) ∨
    (left.x * right.z = right.x * left.z ∧
      left.y * right.z = right.y * left.z)

theorem cpp_complete_add_zero_left [CommRing K] (b : K)
    (point : HomogeneousPoint K) :
    (cppCompleteAdd b cppProjectiveZero point).equivalent point := by
  rw [cpp_complete_add_coordinates]
  simp only [cppProjectiveZero, HomogeneousPoint.equivalent]
  by_cases atInfinity : point.z = 0
  · exact Or.inl ⟨by simp [atInfinity], atInfinity⟩
  · right
    constructor <;> ring

theorem cpp_complete_add_zero_right [CommRing K] (b : K)
    (point : HomogeneousPoint K) :
    (cppCompleteAdd b point cppProjectiveZero).equivalent point := by
  rw [cpp_complete_add_coordinates]
  simp only [cppProjectiveZero, HomogeneousPoint.equivalent]
  by_cases atInfinity : point.z = 0
  · exact Or.inl ⟨by simp [atInfinity], atInfinity⟩
  · right
    constructor <;> ring

theorem cpp_complete_add_commutative [CommRing K] (b : K)
    (left right : HomogeneousPoint K) :
    cppCompleteAdd b left right = cppCompleteAdd b right left := by
  rw [cpp_complete_add_coordinates, cpp_complete_add_coordinates]
  ext <;> simp only <;> ring

theorem cpp_complete_add_inverse [CommRing K] (b : K)
    (point : HomogeneousPoint K) :
    (cppCompleteAdd b point point.negate).equivalent cppProjectiveZero := by
  rw [cpp_complete_add_coordinates]
  simp only [HomogeneousPoint.negate, HomogeneousPoint.equivalent,
    cppProjectiveZero]
  left
  constructor
  · ring
  · trivial

theorem cpp_complete_add_self_eq_double [CommRing K] (b : K)
    (point : HomogeneousPoint K) (onCurve : point.onCurve b) :
    cppCompleteAdd b point point = cppPublicDouble b point := by
  simp only [HomogeneousPoint.onCurve] at onCurve
  rw [cpp_complete_add_coordinates]
  ext
  · simp only [cppPublicDouble]
    ring
  · simp only [cppPublicDouble]
    linear_combination -18 * b * point.z * onCurve
  · simp only [cppPublicDouble]
    linear_combination -6 * point.y * onCurve

abbrev SelectedCurveFp := ZMod Generated.target_modulus

theorem selected_g1_complete_operations
    (b : SelectedCurveFp) (left right : HomogeneousPoint SelectedCurveFp)
    (leftOnCurve : left.onCurve b) (rightOnCurve : right.onCurve b) :
    (cppCompleteAdd b left right).onCurve b ∧
      (cppPublicDouble b left).onCurve b := by
  exact ⟨cpp_complete_add_preserves_curve b left right
      leftOnCurve rightOnCurve,
    cpp_public_double_preserves_curve b left leftOnCurve⟩

theorem selected_g2_complete_operations
    (b : SelectedMillerFp2)
    (left right : HomogeneousPoint SelectedMillerFp2)
    (leftOnCurve : left.onCurve b) (rightOnCurve : right.onCurve b) :
    (cppCompleteAdd b left right).onCurve b ∧
      (cppPublicDouble b left).onCurve b := by
  exact ⟨cpp_complete_add_preserves_curve b left right
      leftOnCurve rightOnCurve,
    cpp_public_double_preserves_curve b left leftOnCurve⟩

end CatBLST.AVX.MatrixNoK
