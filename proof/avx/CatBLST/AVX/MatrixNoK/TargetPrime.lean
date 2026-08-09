/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.TargetPrimeMedium
import CatBLST.AVX.MatrixNoK.TargetPrimeChecksLarge
import CatBLST.AVX.MatrixNoK.TargetPrimeLeaves
import Mathlib.Tactic

/-!
# BLS12-381 base-field primality certificate

`norm_num` proves primality by trial division and is not suitable for the
381-bit target modulus. This file completes a recursive Lucas certificate.
The imported check modules independently verify every modular-order witness;
the theorems below verify the factorizations and recursively establish that
all factors are prime.
-/

namespace CatBLST.AVX.MatrixNoK

theorem prime_92691255082156974996979 :
    Nat.Prime 92691255082156974996979 := by
  apply lucas_of_factor_list 92691255082156974996979 3
    [2, 3, 31, 467, 16447, 64881703735777]
  · norm_num
  · intro q hq
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hq
    rcases hq with rfl | rfl | rfl | rfl | rfl | rfl
    all_goals first
      | exact prime_64881703735777
      | exact target_prime_leaf_factor_prime (by decide)
  · exact PrimeChecksLarge.check_92691255082156974996979

theorem prime_1125266252156850182658904441386709967 :
    Nat.Prime 1125266252156850182658904441386709967 := by
  apply lucas_of_factor_list 1125266252156850182658904441386709967 5
    [2, 3373, 43670061551, 3819663927398918131021]
  · norm_num
  · intro q hq
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hq
    rcases hq with rfl | rfl | rfl | rfl
    all_goals first
      | exact prime_43670061551
      | exact prime_3819663927398918131021
      | exact target_prime_leaf_factor_prime (by decide)
  · exact
      PrimeChecksLarge.check_1125266252156850182658904441386709967

theorem prime_15778400344354997994418419698270088123916926905054652752758194827714659 :
    Nat.Prime
      15778400344354997994418419698270088123916926905054652752758194827714659 := by
  apply lucas_of_factor_list
    15778400344354997994418419698270088123916926905054652752758194827714659
    2 [2, 3, 53, 475709467, 92691255082156974996979,
      1125266252156850182658904441386709967]
  · norm_num
  · intro q hq
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hq
    rcases hq with rfl | rfl | rfl | rfl | rfl | rfl
    all_goals first
      | exact prime_475709467
      | exact prime_92691255082156974996979
      | exact prime_1125266252156850182658904441386709967
      | exact target_prime_leaf_factor_prime (by decide)
  · exact
      PrimeChecksLarge.check_15778400344354997994418419698270088123916926905054652752758194827714659

theorem target_modulus_prime : Nat.Prime Generated.target_modulus := by
  rw [show Generated.target_modulus =
    4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787 by
      rfl]
  apply lucas_of_factor_list
    4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787
    2 [2, 3, 3, 11, 23, 47, 10177, 859267, 52437899,
      2584487767265781317813,
      15778400344354997994418419698270088123916926905054652752758194827714659]
  · norm_num
  · intro q hq
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hq
    rcases hq with rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl |
      rfl | rfl
    all_goals first
      | exact prime_52437899
      | exact prime_2584487767265781317813
      | exact
          prime_15778400344354997994418419698270088123916926905054652752758194827714659
      | exact target_prime_leaf_factor_prime (by decide)
  · exact PrimeChecksLarge.check_target_modulus

instance targetModulusPrimeFact :
    Fact (Nat.Prime Generated.target_modulus) :=
  ⟨target_modulus_prime⟩

end CatBLST.AVX.MatrixNoK
