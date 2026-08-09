/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import Mathlib.Tactic

namespace CatBLST.AVX.MatrixNoK

def targetPrimeLeafFactors : List Nat :=
  [2, 3, 5, 7, 11, 13, 17, 19, 23, 31, 43, 47, 53, 89, 113,
    151, 191, 409, 467, 941, 3373, 7577, 8101, 10177, 16447, 43591,
    421987, 582767, 609743, 755057, 859267]

theorem target_prime_leaf_factor_prime {factor : Nat}
    (member : factor ∈ targetPrimeLeafFactors) : Nat.Prime factor := by
  simp only [targetPrimeLeafFactors, List.mem_cons, List.not_mem_nil,
    or_false] at member
  rcases member with rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl |
    rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl |
    rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl |
    rfl
  all_goals norm_num

end CatBLST.AVX.MatrixNoK
