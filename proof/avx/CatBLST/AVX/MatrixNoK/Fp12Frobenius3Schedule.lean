/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Fp12Frobenius

/-! Finite certificate for the selected `p³` Fp12 source schedule. -/

namespace CatBLST.AVX.MatrixNoK

theorem cpp_frobenius3_basis_schedule : ∀ index : Fin 12,
    cppFp12Frobenius3 (fp12Basis index) =
      cppFp12Frobenius1 (cppFp12Frobenius2 (fp12Basis index)) := by
  native_decide

end CatBLST.AVX.MatrixNoK
