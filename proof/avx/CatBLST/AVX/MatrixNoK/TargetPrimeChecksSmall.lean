/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.TargetPrimeCheck

namespace CatBLST.AVX.MatrixNoK.PrimeChecksSmall

def checks : List LucasCheck :=
  [⟨1686913, 10, [2, 2, 2, 2, 2, 2, 2, 3, 23, 191]⟩,
   ⟨51376543, 3, [2, 3, 7, 151, 8101]⟩,
   ⟨52437899, 2, [2, 43, 609743]⟩,
   ⟨475709467, 2, [2, 3, 47, 1686913]⟩,
   ⟨927093389, 3, [2, 2, 13, 409, 43591]⟩,
   ⟨13090036741, 10, [2, 2, 3, 5, 11, 47, 421987]⟩,
   ⟨43670061551, 7, [2, 5, 5, 17, 51376543]⟩,
   ⟨9272813673901, 2, [2, 2, 3, 5, 5, 7, 7577, 582767]⟩]

theorem checks_valid : checks.all lucasCheckValid = true := by
  native_decide

private theorem check_of_mem (check : LucasCheck) (member : check ∈ checks) :
    lucasCheckValid check = true :=
  (List.all_eq_true.mp checks_valid) check member

theorem check_1686913 :
    lucasCheckValid ⟨1686913, 10, [2, 2, 2, 2, 2, 2, 2, 3, 23, 191]⟩ = true :=
  check_of_mem _ (by decide)
theorem check_51376543 :
    lucasCheckValid ⟨51376543, 3, [2, 3, 7, 151, 8101]⟩ = true :=
  check_of_mem _ (by decide)
theorem check_52437899 :
    lucasCheckValid ⟨52437899, 2, [2, 43, 609743]⟩ = true :=
  check_of_mem _ (by decide)
theorem check_475709467 :
    lucasCheckValid ⟨475709467, 2, [2, 3, 47, 1686913]⟩ = true :=
  check_of_mem _ (by decide)
theorem check_927093389 :
    lucasCheckValid ⟨927093389, 3, [2, 2, 13, 409, 43591]⟩ = true :=
  check_of_mem _ (by decide)
theorem check_13090036741 :
    lucasCheckValid ⟨13090036741, 10, [2, 2, 3, 5, 11, 47, 421987]⟩ = true :=
  check_of_mem _ (by decide)
theorem check_43670061551 :
    lucasCheckValid ⟨43670061551, 7, [2, 5, 5, 17, 51376543]⟩ = true :=
  check_of_mem _ (by decide)
theorem check_9272813673901 :
    lucasCheckValid ⟨9272813673901, 2,
      [2, 2, 3, 5, 5, 7, 7577, 582767]⟩ = true :=
  check_of_mem _ (by decide)

end CatBLST.AVX.MatrixNoK.PrimeChecksSmall
