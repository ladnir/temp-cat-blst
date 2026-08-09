/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Fp12Frobenius3Add
import CatBLST.AVX.MatrixNoK.Fp12Frobenius3Scalar

/-! Base-field linearity of the selected `p³` Fp12 source schedule. -/

set_option maxHeartbeats 800000

namespace CatBLST.AVX.MatrixNoK

open Fp12Concrete

theorem cpp_frobenius3_basis_expansion (x : SelectedFp12) :
    cppFp12Frobenius3 x =
      ∑ index : Fin 12,
        scalar (unpackFlatFp12 x index) *
          cppFp12Frobenius3 (fp12Basis index) := by
  calc
    cppFp12Frobenius3 x =
        cppFp12Frobenius3 (fp12BasisExpansion x) := by
      rw [fp12_basis_expansion]
    _ = ∑ index : Fin 12,
          cppFp12Frobenius3
            (scalar (unpackFlatFp12 x index) * fp12Basis index) := by
      change cppFp12Frobenius3AddHom
          (∑ index : Fin 12,
            scalar (unpackFlatFp12 x index) * fp12Basis index) = _
      rw [map_sum]
      rfl
    _ = ∑ index : Fin 12,
          scalar (unpackFlatFp12 x index) *
            cppFp12Frobenius3 (fp12Basis index) := by
      apply Finset.sum_congr rfl
      intro index _
      exact cpp_frobenius3_scalar_mul _ _

end CatBLST.AVX.MatrixNoK
