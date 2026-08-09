/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Fp12Frobenius2Template

/-! Scalar compatibility of the selected `p²` Fp12 source schedule. -/

namespace CatBLST.AVX.MatrixNoK

open Fp12Concrete

theorem cpp_frobenius2_scalar_mul
    (coefficient : SelectedBaseField) (x : SelectedFp12) :
    cppFp12Frobenius2 (scalar coefficient * x) =
      scalar coefficient * cppFp12Frobenius2 x := by
  apply unpack_flat_fp12_injective
  funext index
  simp only [unpack_cpp_fp12_frobenius2, unpack_flat_fp12_scalar_mul,
    unpack_flat_fp12_scalar_mul_fun, cpp_frobenius2_flat_eq_template,
    frobenius2_flat_template_scale]

end CatBLST.AVX.MatrixNoK
