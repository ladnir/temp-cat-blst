/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.TargetPrimeCheck

namespace CatBLST.AVX.MatrixNoK.PrimeChecksMedium

def checks : List LucasCheck :=
  [⟨64881703735777, 5,
      [2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 927093389]⟩,
   ⟨1928745244171409, 3, [2, 2, 2, 2, 13, 9272813673901]⟩,
   ⟨7259797099061183477, 2, [2, 2, 941, 1928745244171409]⟩,
   ⟨2584487767265781317813, 2, [2, 2, 89, 7259797099061183477]⟩,
   ⟨3819663927398918131021, 6,
      [2, 2, 3, 3, 5, 19, 113, 755057, 13090036741]⟩]

theorem checks_valid : checks.all lucasCheckValid = true := by
  native_decide

private theorem check_of_mem (check : LucasCheck) (member : check ∈ checks) :
    lucasCheckValid check = true :=
  (List.all_eq_true.mp checks_valid) check member

theorem check_64881703735777 :
    lucasCheckValid ⟨64881703735777, 5,
      [2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 927093389]⟩ = true :=
  check_of_mem _ (by decide)
theorem check_1928745244171409 :
    lucasCheckValid ⟨1928745244171409, 3,
      [2, 2, 2, 2, 13, 9272813673901]⟩ = true :=
  check_of_mem _ (by decide)
theorem check_7259797099061183477 :
    lucasCheckValid ⟨7259797099061183477, 2,
      [2, 2, 941, 1928745244171409]⟩ = true :=
  check_of_mem _ (by decide)
theorem check_2584487767265781317813 :
    lucasCheckValid ⟨2584487767265781317813, 2,
      [2, 2, 89, 7259797099061183477]⟩ = true :=
  check_of_mem _ (by decide)
theorem check_3819663927398918131021 :
    lucasCheckValid ⟨3819663927398918131021, 6,
      [2, 2, 3, 3, 5, 19, 113, 755057, 13090036741]⟩ = true :=
  check_of_mem _ (by decide)

end CatBLST.AVX.MatrixNoK.PrimeChecksMedium
