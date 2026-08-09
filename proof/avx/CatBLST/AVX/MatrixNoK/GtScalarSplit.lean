/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Fp12FinalExponentiationPhysical
import Mathlib.Tactic

set_option maxHeartbeats 800000

/-!
# GT unsigned radix-z scalar split

This models the production CT and AVX-VT splitter after scalar preparation:
one division by `z^2` separates two 128-bit halves and one division by `z` in
each half yields four 64-bit digits.  The source stores remainder before
quotient in each pair, exactly as the definition below does.
-/

namespace CatBLST.AVX.MatrixNoK

def gtRadix : Nat := 0xd201000000010000

abbrev GtExponent := ZMod blsScalarOrder

def gtCanonicalScalar (scalar : Nat) : Nat := scalar % blsScalarOrder

def gtRadixZSplit (scalar : Nat) : Fin 4 -> Nat :=
  let low := scalar % (gtRadix ^ 2)
  let high := scalar / (gtRadix ^ 2)
  ![low % gtRadix, low / gtRadix,
    high % gtRadix, high / gtRadix]

theorem gt_radix_positive : 0 < gtRadix := by
  norm_num [gtRadix]

theorem gt_scalar_order_lt_radix_fourth :
    blsScalarOrder < gtRadix ^ 4 := by
  native_decide

theorem gt_radix_z_split_recompose (scalar : Nat) :
    scalar =
      gtRadixZSplit scalar 0 +
        gtRadix * gtRadixZSplit scalar 1 +
        gtRadix ^ 2 * gtRadixZSplit scalar 2 +
        gtRadix ^ 3 * gtRadixZSplit scalar 3 := by
  let low := scalar % (gtRadix ^ 2)
  let high := scalar / (gtRadix ^ 2)
  have splitScalar :
      low + gtRadix ^ 2 * high = scalar := by
    exact Nat.mod_add_div scalar (gtRadix ^ 2)
  have splitLow :
      low % gtRadix + gtRadix * (low / gtRadix) = low := by
    exact Nat.mod_add_div low gtRadix
  have splitHigh :
      high % gtRadix + gtRadix * (high / gtRadix) = high := by
    exact Nat.mod_add_div high gtRadix
  simp only [gtRadixZSplit]
  change scalar =
    low % gtRadix + gtRadix * (low / gtRadix) +
      gtRadix ^ 2 * (high % gtRadix) +
      gtRadix ^ 3 * (high / gtRadix)
  rw [splitLow]
  calc
    scalar = low + gtRadix ^ 2 * high := splitScalar.symm
    _ = low + gtRadix ^ 2 *
        (high % gtRadix + gtRadix * (high / gtRadix)) := by
      rw [splitHigh]
    _ = low + gtRadix ^ 2 * (high % gtRadix) +
        gtRadix ^ 3 * (high / gtRadix) := by ring

theorem gt_radix_z_split_digit_bounds
    (scalar : Nat) (scalarBound : scalar < gtRadix ^ 4) :
    forall index : Fin 4, gtRadixZSplit scalar index < gtRadix := by
  intro index
  fin_cases index
  all_goals simp only [gtRadixZSplit]
  all_goals first
    | exact Nat.mod_lt _ gt_radix_positive
    | apply (Nat.div_lt_iff_lt_mul gt_radix_positive).2
      have lowBound := Nat.mod_lt scalar (pow_pos gt_radix_positive 2)
      simpa [pow_two] using lowBound
    | change (scalar / (gtRadix ^ 2)) / gtRadix < gtRadix
      rw [Nat.div_div_eq_div_mul]
      have radixCubedPositive : 0 < gtRadix ^ 2 * gtRadix :=
        Nat.mul_pos (pow_pos gt_radix_positive 2) gt_radix_positive
      apply (Nat.div_lt_iff_lt_mul radixCubedPositive).2
      simpa [pow_succ, Nat.mul_assoc] using scalarBound

theorem gt_canonical_scalar_bound (scalar : Nat) :
    gtCanonicalScalar scalar < blsScalarOrder := by
  exact Nat.mod_lt scalar (by native_decide : 0 < blsScalarOrder)

theorem gt_canonical_split_digit_bounds (scalar : Nat) :
    forall index : Fin 4,
      gtRadixZSplit (gtCanonicalScalar scalar) index < gtRadix := by
  apply gt_radix_z_split_digit_bounds
  exact (gt_canonical_scalar_bound scalar).trans
    gt_scalar_order_lt_radix_fourth

/-- Frobenius powers used before conjugating odd table bases. -/
def gtFrobeniusBasisExponent : Fin 4 -> GtExponent :=
  ![1, Generated.target_modulus,
    Generated.target_modulus ^ 2, Generated.target_modulus ^ 3]

/-- Exact CT/unsigned-VT source basis after conjugating odd Frobenius powers. -/
def gtUnsignedSourceBasisExponent : Fin 4 -> GtExponent :=
  ![gtFrobeniusBasisExponent 0, -gtFrobeniusBasisExponent 1,
    gtFrobeniusBasisExponent 2, -gtFrobeniusBasisExponent 3]

def gtRadixBasisExponent : Fin 4 -> GtExponent :=
  ![1, gtRadix, gtRadix ^ 2, gtRadix ^ 3]

theorem gt_frobenius_eigenvalue_is_negative_z :
    (Generated.target_modulus : GtExponent) = -(gtRadix : GtExponent) := by
  native_decide

theorem gt_unsigned_source_basis_is_radix_basis :
    gtUnsignedSourceBasisExponent = gtRadixBasisExponent := by
  funext index
  fin_cases index <;>
    native_decide

def gtUnsignedSplitExponent (scalar : Nat) : GtExponent :=
  Finset.univ.sum (fun index : Fin 4 =>
    (gtRadixZSplit (gtCanonicalScalar scalar) index : GtExponent) *
      gtUnsignedSourceBasisExponent index)

theorem gt_unsigned_split_exponent_correct (scalar : Nat) :
    gtUnsignedSplitExponent scalar = (gtCanonicalScalar scalar : GtExponent) := by
  rw [gtUnsignedSplitExponent, gt_unsigned_source_basis_is_radix_basis]
  have recomposed := gt_radix_z_split_recompose (gtCanonicalScalar scalar)
  have castRecomposed := congrArg
    (fun value : Nat => (value : GtExponent)) recomposed
  push_cast at castRecomposed
  simp [gtRadixBasisExponent, Fin.sum_univ_succ]
  convert castRecomposed.symm using 1
  all_goals ring

end CatBLST.AVX.MatrixNoK
