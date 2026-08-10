/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Fp12FinalExponentiationSchedule
import CatBLST.AVX.MatrixNoK.TargetPrimeCheck

namespace CatBLST.AVX.MatrixNoK.ScalarOrderPrimeChecks

def checks : List LucasCheck :=
  [⟨63690073, 7, [2, 2, 2, 3, 2653753]⟩,
   ⟨254760293, 2, [2, 2, 63690073]⟩,
   ⟨blsScalarOrder, 7,
    [2, 3, 11, 19, 10177, 125527, 859267, 906349, 2508409,
      2529403, 52437899, 254760293]⟩]

theorem checks_valid : checks.all lucasCheckValid = true := by
  native_decide

private theorem check_of_mem (check : LucasCheck) (member : check ∈ checks) :
    lucasCheckValid check = true :=
  (List.all_eq_true.mp checks_valid) check member

theorem check_63690073 :
    lucasCheckValid ⟨63690073, 7, [2, 2, 2, 3, 2653753]⟩ = true :=
  check_of_mem _ (by decide)

theorem check_254760293 :
    lucasCheckValid ⟨254760293, 2, [2, 2, 63690073]⟩ = true :=
  check_of_mem _ (by decide)

theorem check_scalar_order :
    lucasCheckValid
      ⟨blsScalarOrder, 7,
        [2, 3, 11, 19, 10177, 125527, 859267, 906349, 2508409,
          2529403, 52437899, 254760293]⟩ = true :=
  check_of_mem _ (by decide)

end CatBLST.AVX.MatrixNoK.ScalarOrderPrimeChecks
