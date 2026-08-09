/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Fp12Frobenius2Add
import CatBLST.AVX.MatrixNoK.Fp12Frobenius2Scalar

/-! Base-field linearity of the selected `p²` Fp12 source schedule. -/

set_option maxHeartbeats 800000

namespace CatBLST.AVX.MatrixNoK

open Fp12Concrete

theorem cpp_frobenius2_basis_expansion (x : SelectedFp12) :
    cppFp12Frobenius2 x =
      ∑ index : Fin 12,
        scalar (unpackFlatFp12 x index) *
          cppFp12Frobenius2 (fp12Basis index) := by
  calc
    cppFp12Frobenius2 x =
        cppFp12Frobenius2 (fp12BasisExpansion x) := by
      rw [fp12_basis_expansion]
    _ = ∑ index : Fin 12,
          cppFp12Frobenius2
            (scalar (unpackFlatFp12 x index) * fp12Basis index) := by
      change cppFp12Frobenius2AddHom
          (∑ index : Fin 12,
            scalar (unpackFlatFp12 x index) * fp12Basis index) = _
      rw [map_sum]
      rfl
    _ = ∑ index : Fin 12,
          scalar (unpackFlatFp12 x index) *
            cppFp12Frobenius2 (fp12Basis index) := by
      apply Finset.sum_congr rfl
      intro index _
      exact cpp_frobenius2_scalar_mul _ _

end CatBLST.AVX.MatrixNoK
