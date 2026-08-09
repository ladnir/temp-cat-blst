/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.TargetPrimeCheck
import Mathlib.Data.List.Prime
import Mathlib.NumberTheory.LucasPrimality

namespace CatBLST.AVX.MatrixNoK

theorem lucas_of_factor_list
    (n witness : Nat) (factors : List Nat)
    (product_eq : factors.prod = n - 1)
    (factors_prime : ∀ q ∈ factors, Nat.Prime q)
    (checks : lucasCheckValid ⟨n, witness, factors⟩ = true) :
    Nat.Prime n := by
  have fullOrder : (witness : ZMod n) ^ (n - 1) = 1 :=
    lucas_check_full_order checks
  apply lucas_primality n witness fullOrder
  intro q qPrime qDivides
  have qMember : q ∈ factors := by
    apply mem_list_primes_of_dvd_prod qPrime.prime
      (fun factor factorMem => (factors_prime factor factorMem).prime)
    rw [product_eq]
    exact qDivides
  exact lucas_check_factor_order checks qMember

end CatBLST.AVX.MatrixNoK
