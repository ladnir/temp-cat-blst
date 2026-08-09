/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Fp12Frobenius2Linear
import CatBLST.AVX.MatrixNoK.Fp12Frobenius2Schedule

/-! The selected `p²` Fp12 source schedule, composed from the accepted `p` map. -/

set_option maxHeartbeats 800000

namespace CatBLST.AVX.MatrixNoK

open Fp12Concrete

theorem selected_frobenius2_basis (index : Fin 12) :
    frobenius SelectedFp12 Generated.target_modulus
        (frobenius SelectedFp12 Generated.target_modulus
          (fp12Basis index)) =
      cppFp12Frobenius2 (fp12Basis index) := by
  rw [← cpp_selected_fp12_frobenius1_correct,
    ← cpp_selected_fp12_frobenius1_correct]
  exact (cpp_frobenius2_basis_schedule index).symm

theorem cpp_selected_fp12_frobenius2_iterate (x : SelectedFp12) :
    cppFp12Frobenius2 x =
      frobenius SelectedFp12 Generated.target_modulus
        (frobenius SelectedFp12 Generated.target_modulus x) := by
  rw [cpp_frobenius2_basis_expansion]
  have mappedExpansion :
      frobenius SelectedFp12 Generated.target_modulus
          (frobenius SelectedFp12 Generated.target_modulus
            (fp12BasisExpansion x)) =
        ∑ index : Fin 12,
          scalar (unpackFlatFp12 x index) *
            cppFp12Frobenius2 (fp12Basis index) := by
    simp only [fp12BasisExpansion, map_sum, map_mul,
      selected_frobenius1_scalar, selected_frobenius2_basis]
  rw [fp12_basis_expansion] at mappedExpansion
  exact mappedExpansion.symm

theorem cpp_selected_fp12_frobenius2_correct (x : SelectedFp12) :
    cppFp12Frobenius2 x = x ^ (Generated.target_modulus ^ 2) := by
  rw [cpp_selected_fp12_frobenius2_iterate]
  simp [frobenius_def, pow_mul, pow_two]

end CatBLST.AVX.MatrixNoK
