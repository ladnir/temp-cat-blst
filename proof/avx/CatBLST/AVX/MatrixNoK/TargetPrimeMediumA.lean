/- Copyright 2026 Category Labs, Inc. SPDX-License-Identifier: Apache-2.0 -/

import CatBLST.AVX.MatrixNoK.TargetPrimeSmall
import CatBLST.AVX.MatrixNoK.TargetPrimeChecksMedium
import CatBLST.AVX.MatrixNoK.TargetPrimeLeaves
import Mathlib.Tactic

namespace CatBLST.AVX.MatrixNoK

theorem prime_64881703735777 : Nat.Prime 64881703735777 := by
  apply lucas_of_factor_list 64881703735777 5
    [2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 927093389]
  · norm_num
  · intro q hq
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hq
    rcases hq with rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl |
      rfl | rfl | rfl | rfl | rfl
    all_goals first
      | exact prime_927093389
      | exact target_prime_leaf_factor_prime (by decide)
  · exact PrimeChecksMedium.check_64881703735777

theorem prime_1928745244171409 : Nat.Prime 1928745244171409 := by
  apply lucas_of_factor_list 1928745244171409 3
    [2, 2, 2, 2, 13, 9272813673901]
  · norm_num
  · intro q hq
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hq
    rcases hq with rfl | rfl | rfl | rfl | rfl | rfl
    all_goals first
      | exact prime_9272813673901
      | exact target_prime_leaf_factor_prime (by decide)
  · exact PrimeChecksMedium.check_1928745244171409

end CatBLST.AVX.MatrixNoK
