/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Fp12Coordinates
import Mathlib.Tactic

/-!
The exact selected `p²` schedule factored into a coefficient-parametric linear
template.  Keeping the coefficients opaque prevents algebraic proofs from
normalizing their 381-bit concrete values.
-/

namespace CatBLST.AVX.MatrixNoK

def frobenius2FlatTemplate {R : Type*} [Neg R] [Mul R]
    (coefficient4 coefficient2 coefficient23 coefficient3 : R)
    (x : FlatFp12 R) : FlatFp12 R :=
  ![x 0, x 1,
    x 2 * coefficient4, x 3 * coefficient4,
    x 4 * coefficient2, x 5 * coefficient2,
    x 6 * coefficient23, x 7 * coefficient23,
    -x 8, -x 9,
    x 10 * coefficient3, x 11 * coefficient3]

theorem frobenius2_flat_template_zero {R : Type*} [CommRing R]
    (coefficient4 coefficient2 coefficient23 coefficient3 : R) :
    frobenius2FlatTemplate coefficient4 coefficient2 coefficient23 coefficient3 0 = 0 := by
  funext index
  fin_cases index <;> simp [frobenius2FlatTemplate]

theorem frobenius2_flat_template_add {R : Type*} [CommRing R]
    (coefficient4 coefficient2 coefficient23 coefficient3 : R)
    (x y : FlatFp12 R) :
    frobenius2FlatTemplate coefficient4 coefficient2 coefficient23 coefficient3 (x + y) =
      frobenius2FlatTemplate coefficient4 coefficient2 coefficient23 coefficient3 x +
        frobenius2FlatTemplate coefficient4 coefficient2 coefficient23 coefficient3 y := by
  funext index
  fin_cases index <;> simp [frobenius2FlatTemplate, add_mul, add_comm]

theorem frobenius2_flat_template_scale {R : Type*} [CommRing R]
    (coefficient4 coefficient2 coefficient23 coefficient3 scalarValue : R)
    (x : FlatFp12 R) :
    frobenius2FlatTemplate coefficient4 coefficient2 coefficient23 coefficient3
        (fun index => scalarValue * x index) =
      fun index => scalarValue *
        frobenius2FlatTemplate coefficient4 coefficient2 coefficient23 coefficient3 x index := by
  funext index
  fin_cases index <;>
    simp [frobenius2FlatTemplate, mul_left_comm, mul_comm]

theorem cpp_frobenius2_flat_eq_template (x : FlatFp12 SelectedBaseField) :
    cppFp12Frobenius2Flat x =
      frobenius2FlatTemplate frob1_4 frob1_2 frob2_3 frob1_3 x := rfl

@[simp] theorem unpack_cpp_fp12_frobenius2
    (x : SelectedFp12) (index : Fin 12) :
    unpackFlatFp12 (cppFp12Frobenius2 x) index =
      cppFp12Frobenius2Flat (unpackFlatFp12 x) index := by
  exact congrFun (unpack_pack_fp12 _) index

end CatBLST.AVX.MatrixNoK
