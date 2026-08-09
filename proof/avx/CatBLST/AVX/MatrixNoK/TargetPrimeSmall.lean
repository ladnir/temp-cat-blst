/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.TargetPrimeLucas
import CatBLST.AVX.MatrixNoK.TargetPrimeChecksSmall
import Mathlib.Tactic

namespace CatBLST.AVX.MatrixNoK

theorem prime_1686913 : Nat.Prime 1686913 := by
  apply lucas_of_factor_list 1686913 10
    [2, 2, 2, 2, 2, 2, 2, 3, 23, 191]
  · norm_num
  · intro q hq
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hq
    rcases hq with rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl
    all_goals norm_num
  · exact PrimeChecksSmall.check_1686913

theorem prime_51376543 : Nat.Prime 51376543 := by
  apply lucas_of_factor_list 51376543 3 [2, 3, 7, 151, 8101]
  · norm_num
  · intro q hq
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hq
    rcases hq with rfl | rfl | rfl | rfl | rfl
    all_goals norm_num
  · exact PrimeChecksSmall.check_51376543

theorem prime_52437899 : Nat.Prime 52437899 := by
  apply lucas_of_factor_list 52437899 2 [2, 43, 609743]
  · norm_num
  · intro q hq
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hq
    rcases hq with rfl | rfl | rfl
    all_goals norm_num
  · exact PrimeChecksSmall.check_52437899

theorem prime_475709467 : Nat.Prime 475709467 := by
  apply lucas_of_factor_list 475709467 2 [2, 3, 47, 1686913]
  · norm_num
  · intro q hq
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hq
    rcases hq with rfl | rfl | rfl | rfl
    all_goals norm_num [prime_1686913]
  · exact PrimeChecksSmall.check_475709467

theorem prime_927093389 : Nat.Prime 927093389 := by
  apply lucas_of_factor_list 927093389 3 [2, 2, 13, 409, 43591]
  · norm_num
  · intro q hq
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hq
    rcases hq with rfl | rfl | rfl | rfl | rfl
    all_goals norm_num
  · exact PrimeChecksSmall.check_927093389

theorem prime_13090036741 : Nat.Prime 13090036741 := by
  apply lucas_of_factor_list 13090036741 10
    [2, 2, 3, 5, 11, 47, 421987]
  · norm_num
  · intro q hq
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hq
    rcases hq with rfl | rfl | rfl | rfl | rfl | rfl | rfl
    all_goals norm_num
  · exact PrimeChecksSmall.check_13090036741

theorem prime_43670061551 : Nat.Prime 43670061551 := by
  apply lucas_of_factor_list 43670061551 7 [2, 5, 5, 17, 51376543]
  · norm_num
  · intro q hq
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hq
    rcases hq with rfl | rfl | rfl | rfl | rfl
    all_goals norm_num [prime_51376543]
  · exact PrimeChecksSmall.check_43670061551

theorem prime_9272813673901 : Nat.Prime 9272813673901 := by
  apply lucas_of_factor_list 9272813673901 2
    [2, 2, 3, 5, 5, 7, 7577, 582767]
  · norm_num
  · intro q hq
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hq
    rcases hq with rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl
    all_goals norm_num
  · exact PrimeChecksSmall.check_9272813673901

end CatBLST.AVX.MatrixNoK
