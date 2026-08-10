/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.ScalarOrderPrimeChecks
import CatBLST.AVX.MatrixNoK.TargetPrimeLucas
import CatBLST.AVX.MatrixNoK.TargetPrimeSmall
import CatBLST.AVX.MatrixNoK.TargetPrimeLeaves
import Mathlib.Tactic

/-!
# BLS12-381 scalar-order primality

`norm_num` proves primality by trial division and deliberately does not scale to
the 255-bit scalar order. This file instead records a Lucas certificate. The
factorization of `r - 1` follows from the BLS parameter identity
`r - 1 = z² (z - 1) (z + 1)`; all of its distinct prime factors are small. The
integer `7` has order `r - 1` modulo `r`.

The modular certificate checks each distinct prime factor once. The separate
factorization list retains multiplicity so its product is exactly `r - 1`.
-/

namespace CatBLST.AVX.MatrixNoK

private theorem lucas_of_factorization_with_distinct_checks
    (n witness : Nat) (factorization checkedFactors : List Nat)
    (product_eq : factorization.prod = n - 1)
    (covered : ∀ q ∈ factorization, q ∈ checkedFactors)
    (checked_prime : ∀ q ∈ checkedFactors, Nat.Prime q)
    (checks : lucasCheckValid ⟨n, witness, checkedFactors⟩ = true) :
    Nat.Prime n := by
  have fullOrder : (witness : ZMod n) ^ (n - 1) = 1 :=
    lucas_check_full_order checks
  apply lucas_primality n witness fullOrder
  intro q qPrime qDivides
  have qInFactorization : q ∈ factorization := by
    apply mem_list_primes_of_dvd_prod qPrime.prime
      (fun factor factorMem =>
        (checked_prime factor (covered factor factorMem)).prime)
    rw [product_eq]
    exact qDivides
  exact lucas_check_factor_order checks (covered q qInFactorization)

private theorem scalar_order_prime_63690073 : Nat.Prime 63690073 := by
  apply lucas_of_factor_list 63690073 7 [2, 2, 2, 3, 2653753]
  · norm_num
  · intro q hq
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hq
    rcases hq with rfl | rfl | rfl | rfl | rfl
    all_goals norm_num
  · exact ScalarOrderPrimeChecks.check_63690073

private theorem scalar_order_prime_254760293 : Nat.Prime 254760293 := by
  apply lucas_of_factor_list 254760293 2 [2, 2, 63690073]
  · norm_num
  · intro q hq
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hq
    rcases hq with rfl | rfl | rfl
    all_goals first
      | exact scalar_order_prime_63690073
      | norm_num
  · exact ScalarOrderPrimeChecks.check_254760293

theorem bls_scalar_order_prime : blsScalarOrder.Prime := by
  apply lucas_of_factorization_with_distinct_checks blsScalarOrder 7
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
      2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
      3, 11, 19, 10177, 125527, 859267, 906349, 906349, 2508409,
      2529403, 52437899, 254760293, 254760293]
    [2, 3, 11, 19, 10177, 125527, 859267, 906349, 2508409,
      2529403, 52437899, 254760293]
  · norm_num [blsScalarOrder]
  · intro q hq
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hq ⊢
    aesop
  · intro q hq
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hq
    rcases hq with rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl | rfl |
      rfl | rfl | rfl
    all_goals first
      | exact prime_52437899
      | exact scalar_order_prime_254760293
      | exact target_prime_leaf_factor_prime (by decide)
      | norm_num
  · exact ScalarOrderPrimeChecks.check_scalar_order

instance blsScalarOrderPrimeFact : Fact blsScalarOrder.Prime :=
  ⟨bls_scalar_order_prime⟩

end CatBLST.AVX.MatrixNoK
