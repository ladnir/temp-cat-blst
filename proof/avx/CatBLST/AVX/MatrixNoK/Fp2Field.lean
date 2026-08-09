/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Fp2Algebra
import CatBLST.AVX.MatrixNoK.TargetPrime
import CatBLST.AVX.MatrixNoK.MillerDouble
import Mathlib.NumberTheory.LegendreSymbol.Basic
import Mathlib.Tactic

/-!
# Field structure for the selected BLS12-381 quadratic extension

`Fp2Value` was already proved to be the ring `Fp[u]/(u² + 1)`. This file
closes the remaining field obligation for the selected base prime. Because the
BLS12-381 field modulus is `3 mod 4`, `-1` is not a square in `Fp`; therefore
the norm `a² + b²` is nonzero for every nonzero pair. The usual conjugate over
norm supplies an inverse without adding any executable implementation claim.
-/

namespace CatBLST.AVX.MatrixNoK

open Fp2Concrete

def selectedFp2Norm (value : SelectedMillerFp2) :
    ZMod Generated.target_modulus :=
  value.x ^ 2 + value.y ^ 2

theorem target_modulus_mod_four : Generated.target_modulus % 4 = 3 := by
  native_decide

theorem selected_fp2_norm_ne_zero (value : SelectedMillerFp2)
    (valueNonzero : value ≠ 0) : selectedFp2Norm value ≠ 0 := by
  intro normZero
  by_cases yZero : value.y = 0
  · have xSquareZero : value.x ^ 2 = 0 := by
      simpa [selectedFp2Norm, yZero] using normZero
    have xZero : value.x = 0 := eq_zero_of_pow_eq_zero xSquareZero
    apply valueNonzero
    ext
    · exact xZero
    · exact yZero
  · have squaresOpposite : value.x ^ 2 = -(value.y ^ 2) := by
      rw [selectedFp2Norm] at normZero
      linear_combination normZero
    exact (ZMod.mod_four_ne_three_of_sq_eq_neg_sq' yZero squaresOpposite)
      target_modulus_mod_four

def selectedFp2Inverse (value : SelectedMillerFp2) : SelectedMillerFp2 :=
  let normInverse := (selectedFp2Norm value)⁻¹
  ⟨value.x * normInverse, -value.y * normInverse⟩

theorem selected_fp2_mul_inverse (value : SelectedMillerFp2)
    (valueNonzero : value ≠ 0) :
    value * selectedFp2Inverse value = 1 := by
  have normNonzero := selected_fp2_norm_ne_zero value valueNonzero
  change Fp2Value.mul value (selectedFp2Inverse value) = Fp2Concrete.one
  ext
  · change value.x * (value.x * (selectedFp2Norm value)⁻¹) -
        value.y * (-value.y * (selectedFp2Norm value)⁻¹) = 1
    calc
      _ = selectedFp2Norm value * (selectedFp2Norm value)⁻¹ := by
        simp only [selectedFp2Norm]
        ring
      _ = 1 := mul_inv_cancel₀ normNonzero
  · change value.x * (-value.y * (selectedFp2Norm value)⁻¹) +
        value.y * (value.x * (selectedFp2Norm value)⁻¹) = 0
    ring

theorem selected_fp2_is_unit_or_zero (value : SelectedMillerFp2) :
    IsUnit value ∨ value = 0 := by
  by_cases valueZero : value = 0
  · exact Or.inr valueZero
  · left
    exact isUnit_iff_exists_inv.mpr
      ⟨selectedFp2Inverse value, selected_fp2_mul_inverse value valueZero⟩

instance selectedMillerFp2Nontrivial : Nontrivial SelectedMillerFp2 := by
  constructor
  refine ⟨0, 1, ?_⟩
  intro equal
  have xEqual := congrArg Fp2Value.x equal
  change (0 : ZMod Generated.target_modulus) = 1 at xEqual
  exact zero_ne_one xEqual

noncomputable instance selectedMillerFp2Field : Field SelectedMillerFp2 :=
  Field.ofIsUnitOrEqZero selected_fp2_is_unit_or_zero

end CatBLST.AVX.MatrixNoK
