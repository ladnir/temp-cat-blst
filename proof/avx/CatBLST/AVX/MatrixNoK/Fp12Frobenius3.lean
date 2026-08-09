/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Fp12Frobenius2
import CatBLST.AVX.MatrixNoK.Fp12Frobenius3Linear
import CatBLST.AVX.MatrixNoK.Fp12Frobenius3Schedule

/-! The selected `p³` Fp12 source schedule, composed from the accepted maps. -/

set_option maxHeartbeats 800000

namespace CatBLST.AVX.MatrixNoK

open Fp12Concrete

theorem selected_frobenius3_basis (index : Fin 12) :
    frobenius SelectedFp12 Generated.target_modulus
        (frobenius SelectedFp12 Generated.target_modulus
          (frobenius SelectedFp12 Generated.target_modulus
            (fp12Basis index))) =
      cppFp12Frobenius3 (fp12Basis index) := by
  rw [← cpp_selected_fp12_frobenius1_correct,
    ← cpp_selected_fp12_frobenius1_correct,
    ← cpp_selected_fp12_frobenius1_correct]
  rw [← cpp_frobenius2_basis_schedule index]
  exact (cpp_frobenius3_basis_schedule index).symm

theorem cpp_selected_fp12_frobenius3_iterate (x : SelectedFp12) :
    cppFp12Frobenius3 x =
      frobenius SelectedFp12 Generated.target_modulus
        (frobenius SelectedFp12 Generated.target_modulus
          (frobenius SelectedFp12 Generated.target_modulus x)) := by
  rw [cpp_frobenius3_basis_expansion]
  have mappedExpansion :
      frobenius SelectedFp12 Generated.target_modulus
          (frobenius SelectedFp12 Generated.target_modulus
            (frobenius SelectedFp12 Generated.target_modulus
              (fp12BasisExpansion x))) =
        ∑ index : Fin 12,
          scalar (unpackFlatFp12 x index) *
            cppFp12Frobenius3 (fp12Basis index) := by
    simp only [fp12BasisExpansion, map_sum, map_mul,
      selected_frobenius1_scalar, selected_frobenius3_basis]
  rw [fp12_basis_expansion] at mappedExpansion
  exact mappedExpansion.symm

theorem cpp_selected_fp12_frobenius3_correct (x : SelectedFp12) :
    cppFp12Frobenius3 x = x ^ (Generated.target_modulus ^ 3) := by
  rw [cpp_selected_fp12_frobenius3_iterate]
  simp [frobenius_def, pow_mul, pow_succ]

end CatBLST.AVX.MatrixNoK
