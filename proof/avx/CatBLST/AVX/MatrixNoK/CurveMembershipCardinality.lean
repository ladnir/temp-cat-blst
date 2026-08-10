/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.CurveMembership
import Mathlib.FieldTheory.Finiteness

/-!
# Cardinality closure for the selected G1 relation kernel

This file closes the remaining G1 membership obligation without assuming the
published BLS12-381 curve order.  An affine Weierstrass curve over a finite
field has at most two nonzero points for each x-coordinate, hence at most
`2 * #F + 1` points including infinity.  The accepted relation kernel is a
vector space over `ZMod r` containing a point of exact order `r`; if it were
larger than `r`, it would contain at least `r^2` points.  For BLS12-381,
`r^2 > 2p + 1`, which is impossible.

The analogous coarse bound over the G2 base field `Fp2` is not strong enough,
so this module intentionally proves only the G1 result.
-/

namespace CatBLST.AVX.MatrixNoK

open WeierstrassCurve
open scoped WeierstrassCurve.Projective

namespace AffineCardinality

variable {F : Type*} [Field F]
variable {W : WeierstrassCurve F} [W.IsElliptic]

private abbrev NonzeroPoint := {point : W.toAffine.Point // point ≠ 0}

private noncomputable instance affinePointFinite [Finite F] :
    Finite W.toAffine.Point := by
  let coordinates : W.toAffine.Point → Option (F × F)
    | .zero => none
    | .some x y _ => some (x, y)
  apply Finite.of_injective coordinates
  intro left right equal
  rcases left with _ | ⟨leftX, leftY, leftValid⟩
  · rcases right with _ | ⟨rightX, rightY, rightValid⟩
    · rfl
    · simp [coordinates] at equal
  · rcases right with _ | ⟨rightX, rightY, rightValid⟩
    · simp [coordinates] at equal
    · simp only [coordinates, Option.some.injEq, Prod.mk.injEq] at equal
      rcases equal with ⟨rfl, rfl⟩
      rfl

private def xCoordinate : NonzeroPoint (W := W) → F
  | ⟨.zero, notZero⟩ => (notZero rfl).elim
  | ⟨.some x _ _, _⟩ => x

omit [W.IsElliptic] in
private theorem same_x_eq_or_eq_neg
    (left right : NonzeroPoint (W := W))
    (sameX : xCoordinate left = xCoordinate right) :
    (left : W.toAffine.Point) = right ∨
      (left : W.toAffine.Point) = -(right : W.toAffine.Point) := by
  rcases left with ⟨(_ | ⟨leftX, leftY, leftValid⟩), leftNotZero⟩
  · exact (leftNotZero rfl).elim
  rcases right with ⟨(_ | ⟨rightX, rightY, rightValid⟩), rightNotZero⟩
  · exact (rightNotZero rfl).elim
  exact (WeierstrassCurve.Affine.Point.X_eq_iff).mp sameX

private theorem x_fiber_card_le_two [Finite F]
    (W : WeierstrassCurve F) [W.IsElliptic] (x : F) :
    Nat.card {point : NonzeroPoint (W := W) // xCoordinate point = x} ≤ 2 := by
  classical
  let fiber := {point : NonzeroPoint (W := W) // xCoordinate point = x}
  letI := Fintype.ofFinite fiber
  rw [Nat.card_eq_fintype_card]
  cases isEmpty_or_nonempty fiber with
  | inl empty =>
      letI : IsEmpty fiber := empty
      simp
  | inr inhabited =>
      letI : Nonempty fiber := inhabited
      let base : fiber := Classical.choice inhabited
      let code : fiber → Bool := fun point => decide (point.1 = base.1)
      apply Fintype.card_le_of_injective code
      intro left right equalCode
      by_cases leftBase : left.1 = base.1
      · have rightBase : right.1 = base.1 := by
          by_contra rightNotBase
          simp [code, leftBase, rightNotBase] at equalCode
        exact Subtype.ext (leftBase.trans rightBase.symm)
      · have rightNotBase : right.1 ≠ base.1 := by
          intro rightBase
          simp [code, leftBase, rightBase] at equalCode
        have leftSameX : xCoordinate left.1 = xCoordinate base.1 :=
          left.property.trans base.property.symm
        have rightSameX : xCoordinate right.1 = xCoordinate base.1 :=
          right.property.trans base.property.symm
        have leftNeg : (left.1.1 : W.toAffine.Point) =
            -(base.1.1 : W.toAffine.Point) := by
          rcases same_x_eq_or_eq_neg left.1 base.1 leftSameX with equal | neg
          · exact (leftBase (Subtype.ext equal)).elim
          · exact neg
        have rightNeg : (right.1.1 : W.toAffine.Point) =
            -(base.1.1 : W.toAffine.Point) := by
          rcases same_x_eq_or_eq_neg right.1 base.1 rightSameX with equal | neg
          · exact (rightNotBase (Subtype.ext equal)).elim
          · exact neg
        exact Subtype.ext (Subtype.ext (leftNeg.trans rightNeg.symm))

theorem affine_point_card_le [Finite F] :
    Nat.card W.toAffine.Point ≤ 2 * Nat.card F + 1 := by
  classical
  letI := Fintype.ofFinite F
  letI := Fintype.ofFinite W.toAffine.Point
  letI := Fintype.ofFinite (NonzeroPoint (W := W))
  have nonzeroBound :
      Fintype.card (NonzeroPoint (W := W)) ≤ 2 * Fintype.card F := by
    change Finset.univ.card ≤ 2 * Fintype.card F
    rw [Finset.card_eq_sum_card_fiberwise
      (s := Finset.univ) (t := Finset.univ)
      (f := xCoordinate (W := W)) (by simp)]
    calc
      _ ≤ ∑ _x ∈ (Finset.univ : Finset F), 2 := by
        apply Finset.sum_le_sum
        intro x _
        have bound := x_fiber_card_le_two W x
        rw [Nat.card_eq_fintype_card, Fintype.card_subtype] at bound
        simpa using bound
      _ = 2 * Fintype.card F := by simp [Nat.mul_comm]
  have pointCard :
      Fintype.card W.toAffine.Point =
        Fintype.card (NonzeroPoint (W := W)) + 1 := by
    rw [← Fintype.card_option]
    exact Fintype.card_congr (Equiv.optionSubtypeNe (0 : W.toAffine.Point)).symm
  simp only [Nat.card_eq_fintype_card]
  rw [pointCard]
  omega

end AffineCardinality

theorem selected_g1_projective_point_card_le :
    Nat.card selectedG1Weierstrass.toProjective.Point ≤
      2 * Generated.target_modulus + 1 := by
  letI : selectedG1Weierstrass.IsElliptic := selectedG1Elliptic
  rw [Nat.card_congr
    (WeierstrassCurve.Projective.Point.toAffineAddEquiv
      selectedG1Weierstrass).toEquiv]
  calc
    Nat.card selectedG1Weierstrass.toAffine.Point ≤
        2 * Nat.card SelectedCurveFp + 1 :=
      AffineCardinality.affine_point_card_le
    _ = 2 * Generated.target_modulus + 1 := by
      rw [Nat.card_zmod]

theorem selected_g1_accepted_card_le_projective :
    Nat.card SelectedG1AcceptedSubgroup ≤
      Nat.card selectedG1Weierstrass.toProjective.Point := by
  exact Nat.card_le_card_of_injective
    (fun point : SelectedG1AcceptedSubgroup =>
      (point : selectedG1Weierstrass.toProjective.Point)) Subtype.coe_injective

theorem scalar_order_squared_gt_g1_coarse_bound :
    2 * Generated.target_modulus + 1 < blsScalarOrder ^ 2 := by
  norm_num [Generated.target_modulus, blsScalarOrder]

theorem selected_g1_accepted_card_le_scalar_order :
    Nat.card SelectedG1AcceptedSubgroup ≤ blsScalarOrder := by
  by_contra notBound
  have strict : blsScalarOrder < Nat.card SelectedG1AcceptedSubgroup :=
    Nat.lt_of_not_ge notBound
  have cardPower : Nat.card SelectedG1AcceptedSubgroup =
      blsScalarOrder ^ Module.finrank CurveExponent SelectedG1AcceptedSubgroup := by
    rw [Module.natCard_eq_pow_finrank
      (K := CurveExponent) (V := SelectedG1AcceptedSubgroup), Nat.card_zmod]
  have rankAtLeastTwo :
      2 ≤ Module.finrank CurveExponent SelectedG1AcceptedSubgroup := by
    by_contra rankSmall
    have rankCases : Module.finrank CurveExponent SelectedG1AcceptedSubgroup = 0 ∨
        Module.finrank CurveExponent SelectedG1AcceptedSubgroup = 1 := by omega
    rcases rankCases with rankZero | rankOne
    · rw [cardPower, rankZero, pow_zero] at strict
      have scalarOrderPositive : 0 < blsScalarOrder := by
        norm_num [blsScalarOrder]
      omega
    · rw [cardPower, rankOne, pow_one] at strict
      exact (Nat.lt_irrefl blsScalarOrder strict).elim
  have squareLeCard : blsScalarOrder ^ 2 ≤
      Nat.card SelectedG1AcceptedSubgroup := by
    rw [cardPower]
    exact pow_le_pow_right' (by norm_num [blsScalarOrder]) rankAtLeastTwo
  have acceptedLeCoarse : Nat.card SelectedG1AcceptedSubgroup ≤
      2 * Generated.target_modulus + 1 :=
    selected_g1_accepted_card_le_projective.trans
      selected_g1_projective_point_card_le
  have coarseLtSquare := scalar_order_squared_gt_g1_coarse_bound
  omega

theorem selected_g1_cyclic_eq_accepted :
    AddSubgroup.zmultiples selectedG1GeneratorPoint =
      SelectedG1AcceptedSubgroup :=
  selected_g1_cyclic_eq_accepted_of_card_le
    selected_g1_accepted_card_le_scalar_order

end CatBLST.AVX.MatrixNoK
