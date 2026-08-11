/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Fp12
import Mathlib.Tactic

set_option maxHeartbeats 800000

/-!
# Source-shaped Fp12 inversion decomposition

This module models the exact algebraic stages in `FP12Ring::fp12_flat_inverse`:

1. take the quadratic norm from Fp12 to Fp6;
2. form the three cubic-extension cofactors and their Fp2 norm;
3. take the Fp2 norm down to the base field;
4. consume one caller-supplied base-field reciprocal; and
5. lift the result back through Fp2, Fp6, and Fp12.

The theorem assumes exactly the external inverter contract used by C++: its
returned base-field value multiplies the computed nonzero norm to one. It does
not assume an unproved generic `Field (Fp12Value R)` instance or replace the
source decomposition with a library inverse.
-/

namespace CatBLST.AVX.MatrixNoK

def fp6Sub {R : Type*} [Sub R]
    (x y : Fp6Value R) : Fp6Value R :=
  { c0 := Fp2Value.sub x.c0 y.c0
    c1 := Fp2Value.sub x.c1 y.c1
    c2 := Fp2Value.sub x.c2 y.c2 }

def fp6Neg {R : Type*} [Neg R] (x : Fp6Value R) : Fp6Value R :=
  { c0 := Fp2Value.neg x.c0
    c1 := Fp2Value.neg x.c1
    c2 := Fp2Value.neg x.c2 }

def fp6One {R : Type*} [Zero R] [One R] : Fp6Value R :=
  { c0 := { x := 1, y := 0 }
    c1 := { x := 0, y := 0 }
    c2 := { x := 0, y := 0 } }

def fp6Zero {R : Type*} [Zero R] : Fp6Value R :=
  { c0 := { x := 0, y := 0 }
    c1 := { x := 0, y := 0 }
    c2 := { x := 0, y := 0 } }

def fp12One {R : Type*} [Zero R] [One R] : Fp12Value R :=
  { a := fp6One
    b := fp6Zero }

/-- C++ `fp12_flat_inverse_to_fp6`: `a^2 - v*b^2`. -/
def cppFp12InverseToFp6 {R : Type*} [CommRing R]
    (x : Fp12Value R) : Fp6Value R :=
  fp6Sub (Fp6Value.mul x.a x.a)
    (Fp6Value.mulByV (Fp6Value.mul x.b x.b))

/-- C++ `fp6_flat_inverse_to_fp2_c0_c1_c2`. -/
def cppFp6InverseCofactors {R : Type*} [CommRing R]
    (x : Fp6Value R) : Fp6Value R :=
  { c0 := Fp2Value.sub (Fp2Value.square x.c0)
      (Fp2Value.mulByNonresidue (Fp2Value.mul x.c1 x.c2))
    c1 := Fp2Value.sub
      (Fp2Value.mulByNonresidue (Fp2Value.square x.c2))
      (Fp2Value.mul x.c0 x.c1)
    c2 := Fp2Value.sub (Fp2Value.square x.c1)
      (Fp2Value.mul x.c0 x.c2) }

/-- C++ `fp6_flat_inverse_to_fp2_ret`. -/
def cppFp6InverseNorm {R : Type*} [CommRing R]
    (cofactors x : Fp6Value R) : Fp2Value R :=
  Fp2Value.add (Fp2Value.mul cofactors.c0 x.c0)
    (Fp2Value.mulByNonresidue
      (Fp2Value.add (Fp2Value.mul cofactors.c1 x.c2)
        (Fp2Value.mul cofactors.c2 x.c1)))

/-- C++ `fp2_flat_inverse_to_fp`. -/
def cppFp2InverseNorm {R : Type*} [CommRing R]
    (x : Fp2Value R) : R := x.x * x.x + x.y * x.y

/-- C++ `fp_flat_inverse_to_fp2`. -/
def cppFp2InverseFromBase {R : Type*} [CommRing R]
    (x : Fp2Value R) (baseInverse : R) : Fp2Value R :=
  { x := x.x * baseInverse
    y := -x.y * baseInverse }

/-- C++ `fp2_flat_inverse_to_fp6`. -/
def cppFp6Scale {R : Type*} [CommRing R]
    (x : Fp6Value R) (scale : Fp2Value R) : Fp6Value R :=
  { c0 := Fp2Value.mul scale x.c0
    c1 := Fp2Value.mul scale x.c1
    c2 := Fp2Value.mul scale x.c2 }

/-- C++ `fp6_flat_inverse_to_fp12`. -/
def cppFp12LiftInverse {R : Type*} [CommRing R]
    (x : Fp12Value R) (inverseNorm : Fp6Value R) : Fp12Value R :=
  { a := Fp6Value.mul x.a inverseNorm
    b := fp6Neg (Fp6Value.mul x.b inverseNorm) }

def cppFp12InversePipeline {R : Type*} [CommRing R]
    (x : Fp12Value R) (baseInverse : R) : Fp12Value R :=
  let fp6Norm := cppFp12InverseToFp6 x
  let cofactors := cppFp6InverseCofactors fp6Norm
  let fp2Norm := cppFp6InverseNorm cofactors fp6Norm
  let fp2Inverse := cppFp2InverseFromBase fp2Norm baseInverse
  let fp6Inverse := cppFp6Scale cofactors fp2Inverse
  cppFp12LiftInverse x fp6Inverse

theorem cpp_fp6_cofactor_product {R : Type*} [CommRing R]
    (x : Fp6Value R) :
    Fp6Value.mul x (cppFp6InverseCofactors x) =
      { c0 := cppFp6InverseNorm (cppFp6InverseCofactors x) x
        c1 := { x := 0, y := 0 }
        c2 := { x := 0, y := 0 } } := by
  ext <;>
    simp [cppFp6InverseCofactors, cppFp6InverseNorm, Fp6Value.mul,
      Fp2Value.add, Fp2Value.sub, Fp2Value.mul, Fp2Value.square,
      Fp2Value.mulByNonresidue] <;>
    ring

theorem cpp_fp2_inverse_from_base_correct {R : Type*} [CommRing R]
    (x : Fp2Value R) (baseInverse : R)
    (inverseCorrect : cppFp2InverseNorm x * baseInverse = 1) :
    Fp2Value.mul x (cppFp2InverseFromBase x baseInverse) =
      { x := 1, y := 0 } := by
  ext
  · calc
      x.x * (x.x * baseInverse) - x.y * (-x.y * baseInverse) =
          cppFp2InverseNorm x * baseInverse := by
            simp [cppFp2InverseNorm]
            ring
      _ = 1 := inverseCorrect
  · simp [cppFp2InverseFromBase, Fp2Value.mul]
    ring

theorem cpp_fp6_mul_scale_right {R : Type*} [CommRing R]
    (x y : Fp6Value R) (scale : Fp2Value R) :
    Fp6Value.mul x (cppFp6Scale y scale) =
      cppFp6Scale (Fp6Value.mul x y) scale := by
  ext <;>
    simp [cppFp6Scale, Fp6Value.mul, Fp2Value.add, Fp2Value.mul,
      Fp2Value.mulByNonresidue] <;>
    ring

theorem fp2_mul_comm {R : Type*} [CommRing R]
    (x y : Fp2Value R) : Fp2Value.mul x y = Fp2Value.mul y x := by
  ext <;> simp [Fp2Value.mul] <;> ring

theorem cpp_fp6_inverse_pipeline_correct {R : Type*} [CommRing R]
    (x : Fp6Value R) (baseInverse : R)
    (inverseCorrect :
      cppFp2InverseNorm
        (cppFp6InverseNorm (cppFp6InverseCofactors x) x) *
          baseInverse = 1) :
    Fp6Value.mul x
      (cppFp6Scale (cppFp6InverseCofactors x)
        (cppFp2InverseFromBase
          (cppFp6InverseNorm (cppFp6InverseCofactors x) x)
      baseInverse)) = fp6One := by
  rw [cpp_fp6_mul_scale_right, cpp_fp6_cofactor_product]
  have fp2Correct := cpp_fp2_inverse_from_base_correct
    (cppFp6InverseNorm (cppFp6InverseCofactors x) x)
    baseInverse inverseCorrect
  unfold cppFp6Scale
  rw [fp2_mul_comm, fp2Correct]
  ext <;> simp [fp6One, Fp2Value.mul]

theorem cpp_fp12_lift_product {R : Type*} [CommRing R]
    (x : Fp12Value R) (inverseNorm : Fp6Value R) :
    Fp12Value.mul x (cppFp12LiftInverse x inverseNorm) =
      { a := Fp6Value.mul (cppFp12InverseToFp6 x) inverseNorm
        b := fp6Zero } := by
  ext <;>
    simp [cppFp12LiftInverse, cppFp12InverseToFp6, fp6Sub, fp6Neg,
      fp6Zero, Fp12Value.mul, Fp6Value.add, Fp6Value.mul,
      Fp6Value.mulByV, Fp2Value.add, Fp2Value.sub, Fp2Value.neg,
      Fp2Value.mul, Fp2Value.mulByNonresidue] <;>
    ring

theorem cpp_fp12_inverse_pipeline_correct {R : Type*} [CommRing R]
    (x : Fp12Value R) (baseInverse : R)
    (inverseCorrect :
      let fp6Norm := cppFp12InverseToFp6 x
      let cofactors := cppFp6InverseCofactors fp6Norm
      let fp2Norm := cppFp6InverseNorm cofactors fp6Norm
      cppFp2InverseNorm fp2Norm * baseInverse = 1) :
    Fp12Value.mul x (cppFp12InversePipeline x baseInverse) = fp12One := by
  unfold cppFp12InversePipeline
  rw [cpp_fp12_lift_product]
  rw [cpp_fp6_inverse_pipeline_correct
    (cppFp12InverseToFp6 x) baseInverse inverseCorrect]
  rfl

end CatBLST.AVX.MatrixNoK
