/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Fp12Frobenius3Template

/-! Additivity of the selected `p³` Fp12 source schedule. -/

namespace CatBLST.AVX.MatrixNoK

def cppFp12Frobenius3AddHom : SelectedFp12 →+ SelectedFp12 where
  toFun := cppFp12Frobenius3
  map_zero' := by
    apply unpack_flat_fp12_injective
    funext index
    simp only [unpack_cpp_fp12_frobenius3, unpack_flat_fp12_zero,
      unpack_flat_fp12_zero_fun, cpp_frobenius3_flat_eq_template,
      frobenius3_flat_template_zero, Pi.zero_apply]
  map_add' := by
    intro x y
    apply unpack_flat_fp12_injective
    funext index
    simp only [unpack_cpp_fp12_frobenius3, unpack_flat_fp12_add,
      unpack_flat_fp12_add_fun, cpp_frobenius3_flat_eq_template,
      frobenius3_flat_template_add, Pi.add_apply]

@[simp] theorem cpp_frobenius3_add_hom_apply (x : SelectedFp12) :
    cppFp12Frobenius3AddHom x = cppFp12Frobenius3 x := rfl

end CatBLST.AVX.MatrixNoK
