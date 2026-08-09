/- Copyright 2026 Category Labs, Inc. SPDX-License-Identifier: Apache-2.0 -/

import CatBLST.AVX.MatrixNoK.TargetPrimeMediumA
import CatBLST.AVX.MatrixNoK.TargetPrimeLeaves
import Mathlib.Tactic

namespace CatBLST.AVX.MatrixNoK

theorem prime_7259797099061183477 : Nat.Prime 7259797099061183477 := by
  apply lucas_of_factor_list 7259797099061183477 2
    [2, 2, 941, 1928745244171409]
  · norm_num
  · intro q hq
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hq
    rcases hq with rfl | rfl | rfl | rfl
    all_goals first
      | exact prime_1928745244171409
      | exact target_prime_leaf_factor_prime (by decide)
  · exact PrimeChecksMedium.check_7259797099061183477

theorem prime_2584487767265781317813 :
    Nat.Prime 2584487767265781317813 := by
  apply lucas_of_factor_list 2584487767265781317813 2
    [2, 2, 89, 7259797099061183477]
  · norm_num
  · intro q hq
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hq
    rcases hq with rfl | rfl | rfl | rfl
    all_goals first
      | exact prime_7259797099061183477
      | exact target_prime_leaf_factor_prime (by decide)
  · exact PrimeChecksMedium.check_2584487767265781317813

end CatBLST.AVX.MatrixNoK
