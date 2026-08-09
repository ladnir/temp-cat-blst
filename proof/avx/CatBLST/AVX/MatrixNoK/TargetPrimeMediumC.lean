/- Copyright 2026 Category Labs, Inc. SPDX-License-Identifier: Apache-2.0 -/

import CatBLST.AVX.MatrixNoK.TargetPrimeMediumB
import CatBLST.AVX.MatrixNoK.TargetPrimeLeaves
import Mathlib.Tactic

namespace CatBLST.AVX.MatrixNoK

theorem prime_3819663927398918131021 :
    Nat.Prime 3819663927398918131021 := by
  apply lucas_of_factor_list 3819663927398918131021 6
    [2, 2, 3, 3, 5, 19, 113, 755057, 13090036741]
  · norm_num
  · intro q hq
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hq
    rcases hq with rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl
    all_goals first
      | exact prime_13090036741
      | exact target_prime_leaf_factor_prime (by decide)
  · exact PrimeChecksMedium.check_3819663927398918131021

end CatBLST.AVX.MatrixNoK
