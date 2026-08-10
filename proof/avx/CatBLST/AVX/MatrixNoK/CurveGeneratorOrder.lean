/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.CurveGenerators
import CatBLST.AVX.MatrixNoK.ScalarOrderPrime
import Mathlib.GroupTheory.OrderOfElement

/-!
# Exact order of the selected BLS12-381 generators

The concrete G1 and G2 scalar-multiplication checks establish that the selected
generators are killed by `r`. The Lucas certificate for `r`, together with the
literal nonzero coordinates, upgrades divisibility to exact additive order.
-/

namespace CatBLST.AVX.MatrixNoK

theorem selected_g1_generator_point_ne_zero :
    selectedG1GeneratorPoint ≠ 0 := by
  intro equalZero
  have equalProjectiveZero :
      toSelectedG1Point selectedG1Generator selected_g1_generator_valid =
        toSelectedG1Point cppProjectiveZero selected_g1_zero_valid :=
    equalZero.trans to_selected_g1_point_zero.symm
  have equivalent := (to_selected_g1_point_eq_iff _ _
    selected_g1_generator_valid selected_g1_zero_valid).mp equalProjectiveZero
  simp [HomogeneousPoint.equivalent, cppProjectiveZero,
    selectedG1Generator] at equivalent

theorem selected_g2_generator_point_ne_zero :
    selectedG2GeneratorPoint ≠ 0 := by
  intro equalZero
  have equalProjectiveZero :
      toSelectedG2Point selectedG2Generator selected_g2_generator_valid =
        toSelectedG2Point cppProjectiveZero selected_g2_zero_valid :=
    equalZero.trans to_selected_g2_point_zero.symm
  have equivalent := (to_selected_g2_point_eq_iff _ _
    selected_g2_generator_valid selected_g2_zero_valid).mp equalProjectiveZero
  simp [HomogeneousPoint.equivalent, cppProjectiveZero,
    selectedG2Generator] at equivalent

theorem selected_g1_generator_point_exact_order :
    addOrderOf selectedG1GeneratorPoint = blsScalarOrder := by
  have divides : addOrderOf selectedG1GeneratorPoint ∣ blsScalarOrder :=
    addOrderOf_dvd_iff_nsmul_eq_zero.mpr
      selected_g1_generator_point_order_divides_r
  rcases (Nat.dvd_prime bls_scalar_order_prime).mp divides with
    orderOne | orderR
  · exact (selected_g1_generator_point_ne_zero
      (AddMonoid.addOrderOf_eq_one_iff.mp orderOne)).elim
  · exact orderR

theorem selected_g2_generator_point_exact_order :
    addOrderOf selectedG2GeneratorPoint = blsScalarOrder := by
  have divides : addOrderOf selectedG2GeneratorPoint ∣ blsScalarOrder :=
    addOrderOf_dvd_iff_nsmul_eq_zero.mpr
      selected_g2_generator_point_order_divides_r
  rcases (Nat.dvd_prime bls_scalar_order_prime).mp divides with
    orderOne | orderR
  · exact (selected_g2_generator_point_ne_zero
      (AddMonoid.addOrderOf_eq_one_iff.mp orderOne)).elim
  · exact orderR

theorem selected_g1_generator_span_natCard :
    Nat.card (AddSubgroup.zmultiples selectedG1GeneratorPoint) =
      blsScalarOrder := by
  rw [Nat.card_zmultiples, selected_g1_generator_point_exact_order]

theorem selected_g2_generator_span_natCard :
    Nat.card (AddSubgroup.zmultiples selectedG2GeneratorPoint) =
      blsScalarOrder := by
  rw [Nat.card_zmultiples, selected_g2_generator_point_exact_order]

end CatBLST.AVX.MatrixNoK
