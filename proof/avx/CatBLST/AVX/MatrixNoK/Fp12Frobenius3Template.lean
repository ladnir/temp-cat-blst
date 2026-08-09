/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Fp12Coordinates
import Mathlib.Tactic

/-!
The exact selected `p³` schedule factored into a coefficient-parametric linear
template, keeping the concrete 381-bit constants opaque to algebraic proofs.
-/

namespace CatBLST.AVX.MatrixNoK

def frobenius3FlatTemplate {R : Type*} [Neg R] [Add R] [Mul R]
    (coefficient7 coefficient8 : R) (x : FlatFp12 R) : FlatFp12 R :=
  ![x 0, -x 1, x 3, x 2, -x 4, x 5,
    x 7 * coefficient7 + x 6 * coefficient8,
    x 7 * coefficient7 + x 6 * coefficient7,
    x 8 * coefficient8 + x 9 * coefficient8,
    x 8 * coefficient8 + x 9 * coefficient7,
    x 11 * coefficient8 + x 10 * coefficient7,
    x 11 * coefficient8 + x 10 * coefficient8]

theorem frobenius3_flat_template_zero {R : Type*} [CommRing R]
    (coefficient7 coefficient8 : R) :
    frobenius3FlatTemplate coefficient7 coefficient8 0 = 0 := by
  funext index
  fin_cases index <;> simp [frobenius3FlatTemplate]

theorem frobenius3_flat_template_add {R : Type*} [CommRing R]
    (coefficient7 coefficient8 : R) (x y : FlatFp12 R) :
    frobenius3FlatTemplate coefficient7 coefficient8 (x + y) =
      frobenius3FlatTemplate coefficient7 coefficient8 x +
        frobenius3FlatTemplate coefficient7 coefficient8 y := by
  funext index
  fin_cases index <;>
    simp [frobenius3FlatTemplate, add_mul, add_assoc, add_left_comm, add_comm]

theorem frobenius3_flat_template_scale {R : Type*} [CommRing R]
    (coefficient7 coefficient8 scalarValue : R) (x : FlatFp12 R) :
    frobenius3FlatTemplate coefficient7 coefficient8
        (fun index => scalarValue * x index) =
      fun index => scalarValue *
        frobenius3FlatTemplate coefficient7 coefficient8 x index := by
  funext index
  fin_cases index <;>
    simp [frobenius3FlatTemplate, mul_add, mul_left_comm, mul_comm]

theorem cpp_frobenius3_flat_eq_template (x : FlatFp12 SelectedBaseField) :
    cppFp12Frobenius3Flat x =
      frobenius3FlatTemplate frob1_7 frob1_8 x := rfl

@[simp] theorem unpack_cpp_fp12_frobenius3
    (x : SelectedFp12) (index : Fin 12) :
    unpackFlatFp12 (cppFp12Frobenius3 x) index =
      cppFp12Frobenius3Flat (unpackFlatFp12 x) index := by
  exact congrFun (unpack_pack_fp12 _) index

end CatBLST.AVX.MatrixNoK
