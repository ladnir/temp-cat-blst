/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.CurveDispatch
import CatBLST.AVX.MatrixNoK.CurveGeneratorOrder
import Mathlib.GroupTheory.OrderOfElement
import Mathlib.Tactic

/-!
# Concrete subgroup-relation completeness

The inherited subgroup tests use short endomorphism relations.  The previous
reduction proved completeness from the full ambient curve cardinality.  This
file takes a narrower route: after removing the prime-order component, the G1
relation is fixed by `sigma²`, whose selected fixed locus is only infinity and
3-torsion, while the G2 relation satisfies `psi²(P) = -P`, whose selected fixed
locus is likewise only infinity and 3-torsion.  The extra scalar factor
`z² + 1` is coprime to 3.
-/

namespace CatBLST.AVX.MatrixNoK

open WeierstrassCurve
open scoped WeierstrassCurve.Projective

def curveRelationExtraFactor : Nat := curveRadix ^ 2 + 1

theorem curve_relation_extra_factor_coprime_three :
    Nat.Coprime curveRelationExtraFactor 3 := by
  native_decide

theorem curve_relation_scalar_factorization :
    curveRadix ^ 6 + 1 =
      blsScalarOrder * curveRelationExtraFactor := by
  exact curve_radix_six_factorization

theorem point_eq_zero_of_coprime_annihilators
    {G : Type*} [AddCommGroup G] (point : G) (left right : Nat)
    (coprime : Nat.Coprime left right)
    (leftKills : left • point = 0) (rightKills : right • point = 0) :
    point = 0 := by
  have leftDivides : addOrderOf point ∣ left :=
    addOrderOf_dvd_iff_nsmul_eq_zero.mpr leftKills
  have rightDivides : addOrderOf point ∣ right :=
    addOrderOf_dvd_iff_nsmul_eq_zero.mpr rightKills
  have dividesOne : addOrderOf point ∣ 1 := by
    simpa [coprime.gcd_eq_one] using Nat.dvd_gcd leftDivides rightDivides
  have orderOne : addOrderOf point = 1 := Nat.eq_one_of_dvd_one dividesOne
  exact AddMonoid.addOrderOf_eq_one_iff.mp orderOne

theorem g1_relation_implies_prime_order_of_small_fixed_kernel
    {G : Type*} [AddCommGroup G]
    (sigma : G →+ G) (point : G)
    (sigmaCube : ∀ value, sigma (sigma (sigma value)) = value)
    (fixedKernel : ∀ value, sigma (sigma value) = value → 3 • value = 0)
    (relation : sigma (sigma point) =
      -((curveRadix ^ 2) • point)) :
    blsScalarOrder • point = 0 := by
  let residual := blsScalarOrder • point
  have productKills :
      (blsScalarOrder * curveRelationExtraFactor) • point = 0 := by
    rw [← curve_relation_scalar_factorization]
    exact g1_relation_annihilates_radix_six_plus_one
      sigma point sigmaCube relation
  have extraKills : curveRelationExtraFactor • residual = 0 := by
    change curveRelationExtraFactor • (blsScalarOrder • point) = 0
    rw [← mul_nsmul, Nat.mul_comm]
    exact productKills
  have residualRelation :
      sigma (sigma residual) = -((curveRadix ^ 2) • residual) := by
    simp only [residual, map_nsmul]
    rw [relation]
    simp only [smul_neg, ← mul_nsmul]
    rw [Nat.mul_comm blsScalarOrder (curveRadix ^ 2)]
  have radixNeg : -((curveRadix ^ 2) • residual) = residual := by
    have sumZero : (curveRadix ^ 2) • residual + residual = 0 := by
      simpa [curveRelationExtraFactor, add_nsmul, one_nsmul] using extraKills
    calc
      -((curveRadix ^ 2) • residual) =
          -((curveRadix ^ 2) • residual) + 0 := by simp
      _ = -((curveRadix ^ 2) • residual) +
          ((curveRadix ^ 2) • residual + residual) := by rw [sumZero]
      _ = residual := by abel
  have fixed : sigma (sigma residual) = residual :=
    residualRelation.trans radixNeg
  exact point_eq_zero_of_coprime_annihilators residual
    curveRelationExtraFactor 3 curve_relation_extra_factor_coprime_three
    extraKills (fixedKernel residual fixed)

theorem g2_relation_implies_prime_order_of_small_fixed_kernel
    {G : Type*} [AddCommGroup G]
    (psi : G →+ G) (point : G)
    (psiSix : ∀ value,
      psi (psi (psi (psi (psi (psi value))))) = -value)
    (fixedKernel : ∀ value, psi (psi value) = -value → 3 • value = 0)
    (relation : psi point = -(curveRadix • point)) :
    blsScalarOrder • point = 0 := by
  let residual := blsScalarOrder • point
  have productKills :
      (blsScalarOrder * curveRelationExtraFactor) • point = 0 := by
    rw [← curve_relation_scalar_factorization]
    exact g2_relation_annihilates_radix_six_plus_one
      psi point psiSix relation
  have extraKills : curveRelationExtraFactor • residual = 0 := by
    change curveRelationExtraFactor • (blsScalarOrder • point) = 0
    rw [← mul_nsmul, Nat.mul_comm]
    exact productKills
  have residualRelation : psi residual = -(curveRadix • residual) := by
    simp only [residual, map_nsmul]
    rw [relation]
    simp only [smul_neg, ← mul_nsmul]
    rw [Nat.mul_comm blsScalarOrder curveRadix]
  have squaredRelation :
      psi (psi residual) = (curveRadix ^ 2) • residual := by
    rw [residualRelation]
    simp only [map_neg, map_nsmul, residualRelation, smul_neg, neg_neg]
    simp [pow_two, mul_nsmul]
  have radixNeg : (curveRadix ^ 2) • residual = -residual := by
    have sumZero : (curveRadix ^ 2) • residual + residual = 0 := by
      simpa [curveRelationExtraFactor, add_nsmul, one_nsmul] using extraKills
    calc
      (curveRadix ^ 2) • residual =
          (curveRadix ^ 2) • residual + 0 := by simp
      _ = (curveRadix ^ 2) • residual + (residual + -residual) := by simp
      _ = ((curveRadix ^ 2) • residual + residual) + -residual := by abel
      _ = -residual := by rw [sumZero, zero_add]
  exact point_eq_zero_of_coprime_annihilators residual
    curveRelationExtraFactor 3 curve_relation_extra_factor_coprime_three
    extraKills (fixedKernel residual (squaredRelation.trans radixNeg))

theorem selected_g1_sigma_squared_fixed_x_zero_or_infinity
    (point : HomogeneousPoint SelectedCurveFp)
    (fixed : (cppG1Sigma (cppG1Sigma point)).equivalent point) :
    point.z = 0 ∨ point.x = 0 := by
  rcases fixed with atInfinity | finite
  · exact Or.inl atInfinity.2
  · by_cases zZero : point.z = 0
    · exact Or.inl zZero
    · right
      have xScaled :
          scalarG1SigmaBeta * (scalarG1SigmaBeta * point.x) = point.x := by
        apply mul_right_cancel₀ zZero
        exact finite.1
      have productZero :
          (scalarG1SigmaBeta ^ 2 - 1) * point.x = 0 := by
        linear_combination xScaled
      exact (mul_eq_zero.mp productZero).resolve_left (by native_decide)

theorem selected_g1_x_zero_three_torsion
    (point : HomogeneousPoint SelectedCurveFp)
    (valid : point.projectivelyValid 4) (xZero : point.x = 0) :
    3 • toSelectedG1Point point valid = 0 := by
  by_cases zZero : point.z = 0
  · have pointZero : toSelectedG1Point point valid = 0 := by
      rw [← to_selected_g1_point_zero]
      apply (to_selected_g1_point_eq_iff _ _ valid selected_g1_zero_valid).2
      exact Or.inl ⟨zZero, rfl⟩
    simp [pointZero]
  · have ySquared : point.y ^ 2 = 4 * point.z ^ 2 := by
      apply mul_right_cancel₀ zZero
      calc
        point.y ^ 2 * point.z = point.x ^ 3 + 4 * point.z ^ 3 := valid.1
        _ = (4 * point.z ^ 2) * point.z := by rw [xZero]; ring
    have yFourth : point.y ^ 4 = 16 * point.z ^ 4 := by
      calc
        point.y ^ 4 = (point.y ^ 2) ^ 2 := by ring
        _ = (4 * point.z ^ 2) ^ 2 := by rw [ySquared]
        _ = 16 * point.z ^ 4 := by ring
    have doubledValid := cpp_complete_add_projectively_valid 4
      selected_g1_no_affine_two_torsion selected_g1_two_ne_zero
      point point valid valid
    have negatedValid := negate_projectively_valid_over 4 point valid
    have doubledEquivalentNegated :
        (cppCompleteAdd 4 point point).equivalent point.negate := by
      rw [cpp_complete_add_coordinates]
      right
      constructor
      · simp only [HomogeneousPoint.negate]
        rw [xZero]
        ring
      · simp only [HomogeneousPoint.negate]
        rw [xZero]
        ring_nf
        rw [yFourth, ySquared]
        ring
    have doubledEqualsNegated :
        toSelectedG1Point (cppCompleteAdd 4 point point) doubledValid =
          toSelectedG1Point point.negate negatedValid :=
      (to_selected_g1_point_eq_iff _ _ doubledValid negatedValid).2
        doubledEquivalentNegated
    have twoEqualsNeg :
        toSelectedG1Point point valid + toSelectedG1Point point valid =
          -toSelectedG1Point point valid := by
      rw [← selected_g1_complete_add_eq_add point point valid valid,
        doubledEqualsNegated, to_selected_g1_point_negate]
    calc
      3 • toSelectedG1Point point valid =
          (toSelectedG1Point point valid + toSelectedG1Point point valid) +
            toSelectedG1Point point valid := by
              rw [three_nsmul]
              abel
      _ = -toSelectedG1Point point valid +
            toSelectedG1Point point valid := by rw [twoEqualsNeg]
      _ = 0 := neg_add_cancel _

theorem selected_g1_sigma_squared_fixed_three_torsion
    (point : HomogeneousPoint SelectedCurveFp)
    (valid : point.projectivelyValid 4)
    (fixed : (cppG1Sigma (cppG1Sigma point)).equivalent point) :
    3 • toSelectedG1Point point valid = 0 := by
  rcases selected_g1_sigma_squared_fixed_x_zero_or_infinity point fixed with
    atInfinity | xZero
  · have pointZero : toSelectedG1Point point valid = 0 := by
      rw [← to_selected_g1_point_zero]
      apply (to_selected_g1_point_eq_iff _ _ valid selected_g1_zero_valid).2
      exact Or.inl ⟨atInfinity, rfl⟩
    simp [pointZero]
  · exact selected_g1_x_zero_three_torsion point valid xZero

theorem g2_psi_squared_x_coefficient_ne_one :
    g2PsiSquaredXCoefficient ≠ 1 := by
  native_decide

theorem selected_g2_psi_squared_neg_fixed_x_zero_or_infinity
    (point : HomogeneousPoint SelectedMillerFp2)
    (_valid : point.projectivelyValid selectedG2CurveB)
    (fixed : (cppG2Psi (cppG2Psi point)).equivalent point.negate) :
    point.z = 0 ∨ point.x = 0 := by
  rw [cpp_g2_psi_squared] at fixed
  rcases fixed with atInfinity | finite
  · exact Or.inl atInfinity.1
  · by_cases zZero : point.z = 0
    · exact Or.inl zZero
    · right
      have xScaled : g2PsiSquaredXCoefficient * point.x = point.x := by
        apply mul_right_cancel₀ zZero
        exact finite.1
      have xProduct : (g2PsiSquaredXCoefficient - 1) * point.x = 0 := by
        linear_combination xScaled
      exact (mul_eq_zero.mp xProduct).resolve_left (sub_ne_zero.mpr
        g2_psi_squared_x_coefficient_ne_one)

theorem selected_g2_x_zero_three_torsion
    (point : HomogeneousPoint SelectedMillerFp2)
    (valid : point.projectivelyValid selectedG2CurveB)
    (xZero : point.x = 0) :
    3 • toSelectedG2Point point valid = 0 := by
  by_cases zZero : point.z = 0
  · have pointZero : toSelectedG2Point point valid = 0 := by
      rw [← to_selected_g2_point_zero]
      apply (to_selected_g2_point_eq_iff _ _ valid selected_g2_zero_valid).2
      exact Or.inl ⟨zZero, rfl⟩
    simp [pointZero]
  · have ySquared : point.y ^ 2 = selectedG2CurveB * point.z ^ 2 := by
      apply mul_right_cancel₀ zZero
      calc
        point.y ^ 2 * point.z =
            point.x ^ 3 + selectedG2CurveB * point.z ^ 3 := valid.1
        _ = (selectedG2CurveB * point.z ^ 2) * point.z := by
          rw [xZero]
          ring
    have yFourth :
        point.y ^ 4 = selectedG2CurveB ^ 2 * point.z ^ 4 := by
      calc
        point.y ^ 4 = (point.y ^ 2) ^ 2 := by ring
        _ = (selectedG2CurveB * point.z ^ 2) ^ 2 := by rw [ySquared]
        _ = selectedG2CurveB ^ 2 * point.z ^ 4 := by ring
    have doubledValid := cpp_complete_add_projectively_valid selectedG2CurveB
      selected_g2_no_affine_two_torsion selected_g2_two_ne_zero
      point point valid valid
    have negatedValid :=
      negate_projectively_valid_over selectedG2CurveB point valid
    have doubledEquivalentNegated :
        (cppCompleteAdd selectedG2CurveB point point).equivalent
          point.negate := by
      rw [cpp_complete_add_coordinates]
      right
      constructor
      · simp only [HomogeneousPoint.negate]
        rw [xZero]
        ring
      · simp only [HomogeneousPoint.negate]
        rw [xZero]
        ring_nf
        rw [yFourth, ySquared]
        ring
    have doubledEqualsNegated :
        toSelectedG2Point
            (cppCompleteAdd selectedG2CurveB point point) doubledValid =
          toSelectedG2Point point.negate negatedValid :=
      (to_selected_g2_point_eq_iff _ _ doubledValid negatedValid).2
        doubledEquivalentNegated
    have twoEqualsNeg :
        toSelectedG2Point point valid + toSelectedG2Point point valid =
          -toSelectedG2Point point valid := by
      rw [← selected_g2_complete_add_eq_add point point valid valid,
        doubledEqualsNegated, to_selected_g2_point_negate]
    calc
      3 • toSelectedG2Point point valid =
          (toSelectedG2Point point valid + toSelectedG2Point point valid) +
            toSelectedG2Point point valid := by
              rw [three_nsmul]
              abel
      _ = -toSelectedG2Point point valid +
            toSelectedG2Point point valid := by rw [twoEqualsNeg]
      _ = 0 := neg_add_cancel _

theorem selected_g2_psi_squared_neg_fixed_three_torsion
    (point : HomogeneousPoint SelectedMillerFp2)
    (valid : point.projectivelyValid selectedG2CurveB)
    (fixed : (cppG2Psi (cppG2Psi point)).equivalent point.negate) :
    3 • toSelectedG2Point point valid = 0 := by
  rcases selected_g2_psi_squared_neg_fixed_x_zero_or_infinity
      point valid fixed with atInfinity | xZero
  · have pointZero : toSelectedG2Point point valid = 0 := by
      rw [← to_selected_g2_point_zero]
      apply (to_selected_g2_point_eq_iff _ _ valid selected_g2_zero_valid).2
      exact Or.inl ⟨atInfinity, rfl⟩
    simp [pointZero]
  · exact selected_g2_x_zero_three_torsion point valid xZero

/-! ## Quotient-safe point endomorphisms

The source predicates operate on homogeneous representatives, while the
subgroup contract is stated on Mathlib's quotient of nonsingular projective
points.  The definitions below are deliberately literal coordinate maps.  In
particular, they introduce no normalized affine representative and therefore
match the C++ source without an inversion or a hidden exceptional case.
-/

def homogeneousOfProjectiveCoordinates (coordinates : Fin 3 → K) :
    HomogeneousPoint K :=
  ⟨coordinates 0, coordinates 1, coordinates 2⟩

@[simp] theorem projective_coordinates_homogeneous_of_projective_coordinates
    (coordinates : Fin 3 → K) :
    projectiveCoordinates (homogeneousOfProjectiveCoordinates coordinates) =
      coordinates := by
  funext coordinate
  fin_cases coordinate <;> rfl

theorem homogeneous_of_projective_coordinates_valid [Field K]
    (b : K) [(shortWeierstrass b).IsElliptic]
    (coordinates : Fin 3 → K)
    (nonsingular :
      (shortWeierstrass b).toProjective.Nonsingular coordinates) :
    (homogeneousOfProjectiveCoordinates coordinates).projectivelyValid b := by
  constructor
  · apply (projective_equation_iff b _).mp
    rw [projective_coordinates_homogeneous_of_projective_coordinates]
    exact nonsingular.1
  · by_contra allZero
    push Not at allZero
    have derivatives :=
      ((shortWeierstrass b).toProjective.nonsingular_iff coordinates).mp
        nonsingular
    simp only [shortWeierstrass] at derivatives
    rcases allZero with ⟨xZero, yZero, zZero⟩
    simp [homogeneousOfProjectiveCoordinates] at xZero yZero zZero
    simp [xZero, yZero, zZero] at derivatives

theorem homogeneous_equivalent_self [CommRing K]
    (point : HomogeneousPoint K) : point.equivalent point := by
  by_cases zZero : point.z = 0
  · exact Or.inl ⟨zZero, zZero⟩
  · exact Or.inr ⟨rfl, rfl⟩

theorem homogeneous_equivalent_symm [CommRing K]
    {left right : HomogeneousPoint K} (equivalent : left.equivalent right) :
    right.equivalent left := by
  rcases equivalent with atInfinity | finite
  · exact Or.inl ⟨atInfinity.2, atInfinity.1⟩
  · exact Or.inr ⟨finite.1.symm, finite.2.symm⟩

def g1SigmaCoordinates (coordinates : Fin 3 → SelectedCurveFp) :
    Fin 3 → SelectedCurveFp :=
  ![scalarG1SigmaBeta * coordinates 0, coordinates 1, coordinates 2]

theorem g1_sigma_coordinates_equiv
    {left right : Fin 3 → SelectedCurveFp} (equivalent : left ≈ right) :
    g1SigmaCoordinates left ≈ g1SigmaCoordinates right := by
  rcases equivalent with ⟨unit, rfl⟩
  refine ⟨unit, ?_⟩
  funext coordinate
  fin_cases coordinate
  all_goals simp [g1SigmaCoordinates, Units.smul_def]
  all_goals ring

def selectedG1SigmaPointClass :
    WeierstrassCurve.Projective.PointClass SelectedCurveFp →
      WeierstrassCurve.Projective.PointClass SelectedCurveFp :=
  Quotient.map g1SigmaCoordinates (fun _ _ => g1_sigma_coordinates_equiv)

theorem cpp_g1_sigma_projectively_valid
    (point : HomogeneousPoint SelectedCurveFp)
    (valid : point.projectivelyValid 4) :
    (cppG1Sigma point).projectivelyValid 4 := by
  constructor
  · simp only [HomogeneousPoint.onCurve, cppG1Sigma] at valid ⊢
    rw [mul_pow, scalar_g1_sigma_beta_cube, one_mul]
    exact valid.1
  · rcases valid.2 with xNonzero | yNonzero | zNonzero
    · left
      exact mul_ne_zero (by native_decide) xNonzero
    · exact Or.inr (Or.inl yNonzero)
    · exact Or.inr (Or.inr zNonzero)

theorem selected_g1_sigma_point_class_nonsingular
    (point : WeierstrassCurve.Projective.PointClass SelectedCurveFp)
    (nonsingular :
      selectedG1Weierstrass.toProjective.NonsingularLift point) :
    selectedG1Weierstrass.toProjective.NonsingularLift
      (selectedG1SigmaPointClass point) := by
  revert nonsingular
  refine Quotient.inductionOn point ?_
  intro coordinates nonsingular
  have sourceNonsingular :
      selectedG1Weierstrass.toProjective.Nonsingular coordinates :=
    nonsingular
  have sourceValid := homogeneous_of_projective_coordinates_valid 4
    coordinates sourceNonsingular
  have targetValid := cpp_g1_sigma_projectively_valid
    (homogeneousOfProjectiveCoordinates coordinates) sourceValid
  exact (WeierstrassCurve.Projective.nonsingularLift_iff _).mpr (by
    simpa [selectedG1SigmaPointClass, g1SigmaCoordinates, cppG1Sigma,
      homogeneousOfProjectiveCoordinates, projectiveCoordinates] using
        selected_g1_valid_nonsingular _ targetValid)

noncomputable def selectedG1SigmaPoint :
    selectedG1Weierstrass.toProjective.Point →
      selectedG1Weierstrass.toProjective.Point :=
  fun point => ⟨selected_g1_sigma_point_class_nonsingular point.point
    point.nonsingular⟩

theorem selected_g1_sigma_point_to_selected
    (point : HomogeneousPoint SelectedCurveFp)
    (valid : point.projectivelyValid 4) :
    selectedG1SigmaPoint (toSelectedG1Point point valid) =
      toSelectedG1Point (cppG1Sigma point)
        (cpp_g1_sigma_projectively_valid point valid) := by
  apply WeierstrassCurve.Projective.Point.ext
  rfl

theorem selected_g1_point_has_homogeneous_representative
    (point : selectedG1Weierstrass.toProjective.Point) :
    ∃ (coordinates : HomogeneousPoint SelectedCurveFp)
      (valid : coordinates.projectivelyValid 4),
      point = toSelectedG1Point coordinates valid := by
  classical
  obtain ⟨affine, rfl⟩ :=
    (WeierstrassCurve.Projective.Point.toAffineAddEquiv
      selectedG1Weierstrass.toProjective).symm.surjective point
  cases affine with
  | zero =>
      exact ⟨cppProjectiveZero, selected_g1_zero_valid, by
        change (0 : selectedG1Weierstrass.toProjective.Point) =
          toSelectedG1Point cppProjectiveZero selected_g1_zero_valid
        exact to_selected_g1_point_zero.symm⟩
  | some x y nonsingular =>
      let coordinates : HomogeneousPoint SelectedCurveFp := ⟨x, y, 1⟩
      have projectiveNonsingular :
          selectedG1Weierstrass.toProjective.Nonsingular
            (projectiveCoordinates coordinates) := by
        simpa [coordinates, projectiveCoordinates] using
          (WeierstrassCurve.Projective.nonsingular_some x y).mpr nonsingular
      have valid := homogeneous_of_projective_coordinates_valid 4
        (projectiveCoordinates coordinates) projectiveNonsingular
      refine ⟨coordinates, valid, ?_⟩
      apply WeierstrassCurve.Projective.Point.ext
      rfl

theorem selected_g1_sigma_point_zero :
    selectedG1SigmaPoint 0 = 0 := by
  rw [← to_selected_g1_point_zero,
    selected_g1_sigma_point_to_selected]
  apply (to_selected_g1_point_eq_iff _ _ _ _).2
  left
  constructor <;> rfl

theorem selected_g1_sigma_point_add
    (left right : selectedG1Weierstrass.toProjective.Point) :
    selectedG1SigmaPoint (left + right) =
      selectedG1SigmaPoint left + selectedG1SigmaPoint right := by
  rcases selected_g1_point_has_homogeneous_representative left with
    ⟨leftCoordinates, leftValid, rfl⟩
  rcases selected_g1_point_has_homogeneous_representative right with
    ⟨rightCoordinates, rightValid, rfl⟩
  let sumValid := cpp_complete_add_projectively_valid 4
    selected_g1_no_affine_two_torsion selected_g1_two_ne_zero
    leftCoordinates rightCoordinates leftValid rightValid
  let sigmaLeftValid := cpp_g1_sigma_projectively_valid
    leftCoordinates leftValid
  let sigmaRightValid := cpp_g1_sigma_projectively_valid
    rightCoordinates rightValid
  let sigmaSumValid := cpp_g1_sigma_projectively_valid
    (cppCompleteAdd 4 leftCoordinates rightCoordinates) sumValid
  let targetSumValid := cpp_complete_add_projectively_valid 4
    selected_g1_no_affine_two_torsion selected_g1_two_ne_zero
    (cppG1Sigma leftCoordinates) (cppG1Sigma rightCoordinates)
    sigmaLeftValid sigmaRightValid
  calc
    selectedG1SigmaPoint
        (toSelectedG1Point leftCoordinates leftValid +
          toSelectedG1Point rightCoordinates rightValid) =
      selectedG1SigmaPoint
        (toSelectedG1Point
          (cppCompleteAdd 4 leftCoordinates rightCoordinates) sumValid) := by
            rw [selected_g1_complete_add_eq_add]
    _ = toSelectedG1Point
        (cppG1Sigma (cppCompleteAdd 4 leftCoordinates rightCoordinates))
          sigmaSumValid := selected_g1_sigma_point_to_selected _ _
    _ = toSelectedG1Point
        (cppCompleteAdd 4 (cppG1Sigma leftCoordinates)
          (cppG1Sigma rightCoordinates)) targetSumValid := by
            apply (to_selected_g1_point_eq_iff _ _ _ _).2
            rw [cpp_g1_sigma_complete_add]
            exact homogeneous_equivalent_self
              (point := cppCompleteAdd 4 (cppG1Sigma leftCoordinates)
                (cppG1Sigma rightCoordinates))
    _ = toSelectedG1Point (cppG1Sigma leftCoordinates) sigmaLeftValid +
        toSelectedG1Point (cppG1Sigma rightCoordinates) sigmaRightValid :=
          selected_g1_complete_add_eq_add _ _ _ _
    _ = selectedG1SigmaPoint (toSelectedG1Point leftCoordinates leftValid) +
        selectedG1SigmaPoint
          (toSelectedG1Point rightCoordinates rightValid) := by
            rw [selected_g1_sigma_point_to_selected,
              selected_g1_sigma_point_to_selected]

noncomputable def selectedG1SigmaHom :
    selectedG1Weierstrass.toProjective.Point →+
      selectedG1Weierstrass.toProjective.Point where
  toFun := selectedG1SigmaPoint
  map_zero' := selected_g1_sigma_point_zero
  map_add' := selected_g1_sigma_point_add

theorem selected_g1_sigma_hom_cube
    (point : selectedG1Weierstrass.toProjective.Point) :
    selectedG1SigmaHom
        (selectedG1SigmaHom (selectedG1SigmaHom point)) = point := by
  rcases selected_g1_point_has_homogeneous_representative point with
    ⟨coordinates, valid, rfl⟩
  change selectedG1SigmaPoint
      (selectedG1SigmaPoint
        (selectedG1SigmaPoint (toSelectedG1Point coordinates valid))) =
    toSelectedG1Point coordinates valid
  rw [selected_g1_sigma_point_to_selected,
    selected_g1_sigma_point_to_selected,
    selected_g1_sigma_point_to_selected]
  apply (to_selected_g1_point_eq_iff _ _ _ _).2
  rw [cpp_g1_sigma_cube]
  exact homogeneous_equivalent_self (point := coordinates)

theorem selected_g1_sigma_hom_squared_fixed_three_torsion
    (point : selectedG1Weierstrass.toProjective.Point)
    (fixed : selectedG1SigmaHom (selectedG1SigmaHom point) = point) :
    3 • point = 0 := by
  rcases selected_g1_point_has_homogeneous_representative point with
    ⟨coordinates, valid, rfl⟩
  change selectedG1SigmaPoint
      (selectedG1SigmaPoint (toSelectedG1Point coordinates valid)) =
    toSelectedG1Point coordinates valid at fixed
  rw [selected_g1_sigma_point_to_selected,
    selected_g1_sigma_point_to_selected] at fixed
  have fixedCoordinates :
      (cppG1Sigma (cppG1Sigma coordinates)).equivalent coordinates :=
    (to_selected_g1_point_eq_iff _ _ _ _).mp fixed
  exact selected_g1_sigma_squared_fixed_three_torsion
    coordinates valid fixedCoordinates

theorem selected_g1_source_relation_implies_prime_order
    (point : HomogeneousPoint SelectedCurveFp)
    (valid : point.projectivelyValid 4)
    (relation :
      (cppG1Sigma (cppG1Sigma point)).equivalent
        (cppBinaryMul 4 point (curveRadix ^ 2)).negate) :
    blsScalarOrder • toSelectedG1Point point valid = 0 := by
  apply g1_relation_implies_prime_order_of_small_fixed_kernel
    selectedG1SigmaHom (toSelectedG1Point point valid)
    selected_g1_sigma_hom_cube
    selected_g1_sigma_hom_squared_fixed_three_torsion
  change selectedG1SigmaPoint
      (selectedG1SigmaPoint (toSelectedG1Point point valid)) =
    -((curveRadix ^ 2) • toSelectedG1Point point valid)
  rw [selected_g1_sigma_point_to_selected,
    selected_g1_sigma_point_to_selected]
  rw [← to_selected_g1_point_cpp_binary_mul point valid (curveRadix ^ 2),
    ← to_selected_g1_point_negate]
  exact (to_selected_g1_point_eq_iff _ _ _ _).mpr relation

def g2PsiCoordinates (coordinates : Fin 3 → SelectedMillerFp2) :
    Fin 3 → SelectedMillerFp2 :=
  projectiveCoordinates
    (cppG2Psi (homogeneousOfProjectiveCoordinates coordinates))

theorem g2_psi_coordinates_equiv
    {left right : Fin 3 → SelectedMillerFp2} (equivalent : left ≈ right) :
    g2PsiCoordinates left ≈ g2PsiCoordinates right := by
  rcases equivalent with ⟨unit, rfl⟩
  refine ⟨Units.map fp2ConjugateHom unit, ?_⟩
  funext coordinate
  fin_cases coordinate <;>
    simp [g2PsiCoordinates, cppG2Psi, homogeneousOfProjectiveCoordinates,
      projectiveCoordinates, Units.smul_def, fp2ConjugateHom,
      fp2_conjugate_mul] <;> ring

def selectedG2PsiPointClass :
    WeierstrassCurve.Projective.PointClass SelectedMillerFp2 →
      WeierstrassCurve.Projective.PointClass SelectedMillerFp2 :=
  Quotient.map g2PsiCoordinates (fun _ _ => g2_psi_coordinates_equiv)

theorem fp2_conjugate_ne_zero {value : SelectedMillerFp2}
    (nonzero : value ≠ 0) : fp2Conjugate value ≠ 0 := by
  intro conjugateZero
  apply nonzero
  apply Fp2Value.ext
  · simpa [fp2Conjugate] using congrArg Fp2Value.x conjugateZero
  · have imaginary := congrArg Fp2Value.y conjugateZero
    change -value.y = 0 at imaginary
    exact neg_eq_zero.mp imaginary

theorem cpp_g2_psi_projectively_valid
    (point : HomogeneousPoint SelectedMillerFp2)
    (valid : point.projectivelyValid selectedG2CurveB) :
    (cppG2Psi point).projectivelyValid selectedG2CurveB := by
  constructor
  · exact cpp_g2_psi_preserves_curve point valid.1
  · rcases valid.2 with xNonzero | yNonzero | zNonzero
    · left
      exact mul_ne_zero (by native_decide) (fp2_conjugate_ne_zero xNonzero)
    · right
      left
      exact mul_ne_zero g2_psi_y_coefficient_nonzero
        (fp2_conjugate_ne_zero yNonzero)
    · exact Or.inr (Or.inr (fp2_conjugate_ne_zero zNonzero))

theorem selected_g2_psi_point_class_nonsingular
    (point : WeierstrassCurve.Projective.PointClass SelectedMillerFp2)
    (nonsingular :
      selectedG2Weierstrass.toProjective.NonsingularLift point) :
    selectedG2Weierstrass.toProjective.NonsingularLift
      (selectedG2PsiPointClass point) := by
  revert nonsingular
  refine Quotient.inductionOn point ?_
  intro coordinates nonsingular
  have sourceNonsingular :
      selectedG2Weierstrass.toProjective.Nonsingular coordinates :=
    nonsingular
  have sourceValid := homogeneous_of_projective_coordinates_valid
    selectedG2CurveB coordinates sourceNonsingular
  have targetValid := cpp_g2_psi_projectively_valid
    (homogeneousOfProjectiveCoordinates coordinates) sourceValid
  exact (WeierstrassCurve.Projective.nonsingularLift_iff _).mpr (by
    simpa [selectedG2PsiPointClass, g2PsiCoordinates] using
      selected_g2_valid_nonsingular _ targetValid)

noncomputable def selectedG2PsiPoint :
    selectedG2Weierstrass.toProjective.Point →
      selectedG2Weierstrass.toProjective.Point :=
  fun point => ⟨selected_g2_psi_point_class_nonsingular point.point
    point.nonsingular⟩

theorem selected_g2_psi_point_to_selected
    (point : HomogeneousPoint SelectedMillerFp2)
    (valid : point.projectivelyValid selectedG2CurveB) :
    selectedG2PsiPoint (toSelectedG2Point point valid) =
      toSelectedG2Point (cppG2Psi point)
        (cpp_g2_psi_projectively_valid point valid) := by
  apply WeierstrassCurve.Projective.Point.ext
  rfl

theorem selected_g2_point_has_homogeneous_representative
    (point : selectedG2Weierstrass.toProjective.Point) :
    ∃ (coordinates : HomogeneousPoint SelectedMillerFp2)
      (valid : coordinates.projectivelyValid selectedG2CurveB),
      point = toSelectedG2Point coordinates valid := by
  classical
  obtain ⟨affine, rfl⟩ :=
    (WeierstrassCurve.Projective.Point.toAffineAddEquiv
      selectedG2Weierstrass.toProjective).symm.surjective point
  cases affine with
  | zero =>
      exact ⟨cppProjectiveZero, selected_g2_zero_valid, by
        change (0 : selectedG2Weierstrass.toProjective.Point) =
          toSelectedG2Point cppProjectiveZero selected_g2_zero_valid
        exact to_selected_g2_point_zero.symm⟩
  | some x y nonsingular =>
      let coordinates : HomogeneousPoint SelectedMillerFp2 := ⟨x, y, 1⟩
      have projectiveNonsingular :
          selectedG2Weierstrass.toProjective.Nonsingular
            (projectiveCoordinates coordinates) := by
        simpa [coordinates, projectiveCoordinates] using
          (WeierstrassCurve.Projective.nonsingular_some x y).mpr nonsingular
      have valid := homogeneous_of_projective_coordinates_valid
        selectedG2CurveB (projectiveCoordinates coordinates)
          projectiveNonsingular
      refine ⟨coordinates, valid, ?_⟩
      apply WeierstrassCurve.Projective.Point.ext
      rfl

private instance selectedCurveHomogeneousFinite :
    Finite (HomogeneousPoint SelectedCurveFp) :=
  Finite.of_injective (fun point => (point.x, point.y, point.z)) (by
    rintro ⟨leftX, leftY, leftZ⟩ ⟨rightX, rightY, rightZ⟩ equal
    simp only [Prod.mk.injEq] at equal
    rcases equal with ⟨rfl, rfl, rfl⟩
    rfl)

private instance selectedMillerFp2Finite : Finite SelectedMillerFp2 :=
  Finite.of_injective (fun value => (value.x, value.y)) (by
    rintro ⟨leftX, leftY⟩ ⟨rightX, rightY⟩ equal
    simp only [Prod.mk.injEq] at equal
    rcases equal with ⟨rfl, rfl⟩
    rfl)

private instance selectedMillerHomogeneousFinite :
    Finite (HomogeneousPoint SelectedMillerFp2) :=
  Finite.of_injective (fun point => (point.x, point.y, point.z)) (by
    rintro ⟨leftX, leftY, leftZ⟩ ⟨rightX, rightY, rightZ⟩ equal
    simp only [Prod.mk.injEq] at equal
    rcases equal with ⟨rfl, rfl, rfl⟩
    rfl)

private instance selectedG1ProjectivePointFinite :
    Finite selectedG1Weierstrass.toProjective.Point := by
  let toPoint :
      {coordinates : HomogeneousPoint SelectedCurveFp //
        coordinates.projectivelyValid 4} →
        selectedG1Weierstrass.toProjective.Point :=
    fun coordinates => toSelectedG1Point coordinates coordinates.property
  apply Finite.of_surjective toPoint
  intro point
  rcases selected_g1_point_has_homogeneous_representative point with
    ⟨coordinates, valid, equal⟩
  exact ⟨⟨coordinates, valid⟩, equal.symm⟩

private instance selectedG2ProjectivePointFinite :
    Finite selectedG2Weierstrass.toProjective.Point := by
  let toPoint :
      {coordinates : HomogeneousPoint SelectedMillerFp2 //
        coordinates.projectivelyValid selectedG2CurveB} →
        selectedG2Weierstrass.toProjective.Point :=
    fun coordinates => toSelectedG2Point coordinates coordinates.property
  apply Finite.of_surjective toPoint
  intro point
  rcases selected_g2_point_has_homogeneous_representative point with
    ⟨coordinates, valid, equal⟩
  exact ⟨⟨coordinates, valid⟩, equal.symm⟩

theorem selected_g2_psi_point_zero :
    selectedG2PsiPoint 0 = 0 := by
  rw [← to_selected_g2_point_zero,
    selected_g2_psi_point_to_selected]
  apply (to_selected_g2_point_eq_iff _ _ _ _).2
  left
  constructor <;> rfl

theorem selected_g2_psi_point_add
    (left right : selectedG2Weierstrass.toProjective.Point) :
    selectedG2PsiPoint (left + right) =
      selectedG2PsiPoint left + selectedG2PsiPoint right := by
  rcases selected_g2_point_has_homogeneous_representative left with
    ⟨leftCoordinates, leftValid, rfl⟩
  rcases selected_g2_point_has_homogeneous_representative right with
    ⟨rightCoordinates, rightValid, rfl⟩
  let sumValid := cpp_complete_add_projectively_valid selectedG2CurveB
    selected_g2_no_affine_two_torsion selected_g2_two_ne_zero
    leftCoordinates rightCoordinates leftValid rightValid
  let psiLeftValid := cpp_g2_psi_projectively_valid
    leftCoordinates leftValid
  let psiRightValid := cpp_g2_psi_projectively_valid
    rightCoordinates rightValid
  let psiSumValid := cpp_g2_psi_projectively_valid
    (cppCompleteAdd selectedG2CurveB leftCoordinates rightCoordinates) sumValid
  let targetSumValid := cpp_complete_add_projectively_valid selectedG2CurveB
    selected_g2_no_affine_two_torsion selected_g2_two_ne_zero
    (cppG2Psi leftCoordinates) (cppG2Psi rightCoordinates)
    psiLeftValid psiRightValid
  calc
    selectedG2PsiPoint
        (toSelectedG2Point leftCoordinates leftValid +
          toSelectedG2Point rightCoordinates rightValid) =
      selectedG2PsiPoint
        (toSelectedG2Point
          (cppCompleteAdd selectedG2CurveB leftCoordinates rightCoordinates)
            sumValid) := by rw [selected_g2_complete_add_eq_add]
    _ = toSelectedG2Point
        (cppG2Psi
          (cppCompleteAdd selectedG2CurveB leftCoordinates rightCoordinates))
          psiSumValid := selected_g2_psi_point_to_selected _ _
    _ = toSelectedG2Point
        (cppCompleteAdd selectedG2CurveB (cppG2Psi leftCoordinates)
          (cppG2Psi rightCoordinates)) targetSumValid := by
            apply (to_selected_g2_point_eq_iff _ _ _ _).2
            exact homogeneous_equivalent_symm
              (cpp_g2_psi_complete_add leftCoordinates rightCoordinates)
    _ = toSelectedG2Point (cppG2Psi leftCoordinates) psiLeftValid +
        toSelectedG2Point (cppG2Psi rightCoordinates) psiRightValid :=
          selected_g2_complete_add_eq_add _ _ _ _
    _ = selectedG2PsiPoint (toSelectedG2Point leftCoordinates leftValid) +
        selectedG2PsiPoint
          (toSelectedG2Point rightCoordinates rightValid) := by
            rw [selected_g2_psi_point_to_selected,
              selected_g2_psi_point_to_selected]

noncomputable def selectedG2PsiHom :
    selectedG2Weierstrass.toProjective.Point →+
      selectedG2Weierstrass.toProjective.Point where
  toFun := selectedG2PsiPoint
  map_zero' := selected_g2_psi_point_zero
  map_add' := selected_g2_psi_point_add

theorem selected_g2_psi_hom_six
    (point : selectedG2Weierstrass.toProjective.Point) :
    selectedG2PsiHom (selectedG2PsiHom (selectedG2PsiHom
      (selectedG2PsiHom (selectedG2PsiHom (selectedG2PsiHom point))))) =
        -point := by
  rcases selected_g2_point_has_homogeneous_representative point with
    ⟨coordinates, valid, rfl⟩
  change selectedG2PsiPoint (selectedG2PsiPoint (selectedG2PsiPoint
      (selectedG2PsiPoint (selectedG2PsiPoint
        (selectedG2PsiPoint (toSelectedG2Point coordinates valid)))))) =
    -toSelectedG2Point coordinates valid
  rw [selected_g2_psi_point_to_selected,
    selected_g2_psi_point_to_selected,
    selected_g2_psi_point_to_selected,
    selected_g2_psi_point_to_selected,
    selected_g2_psi_point_to_selected,
    selected_g2_psi_point_to_selected,
    ← to_selected_g2_point_negate]
  apply (to_selected_g2_point_eq_iff _ _ _ _).2
  rw [cpp_g2_psi_six]
  exact homogeneous_equivalent_self (point := coordinates.negate)

theorem selected_g2_psi_hom_squared_neg_fixed_three_torsion
    (point : selectedG2Weierstrass.toProjective.Point)
    (fixed : selectedG2PsiHom (selectedG2PsiHom point) = -point) :
    3 • point = 0 := by
  rcases selected_g2_point_has_homogeneous_representative point with
    ⟨coordinates, valid, rfl⟩
  change selectedG2PsiPoint
      (selectedG2PsiPoint (toSelectedG2Point coordinates valid)) =
    -toSelectedG2Point coordinates valid at fixed
  rw [selected_g2_psi_point_to_selected,
    selected_g2_psi_point_to_selected,
    ← to_selected_g2_point_negate] at fixed
  have fixedCoordinates :
      (cppG2Psi (cppG2Psi coordinates)).equivalent coordinates.negate :=
    (to_selected_g2_point_eq_iff _ _ _ _).mp fixed
  exact selected_g2_psi_squared_neg_fixed_three_torsion
    coordinates valid fixedCoordinates

theorem selected_g2_source_relation_implies_prime_order
    (point : HomogeneousPoint SelectedMillerFp2)
    (valid : point.projectivelyValid selectedG2CurveB)
    (relation :
      (cppG2Psi point).equivalent
        (cppBinaryMul selectedG2CurveB point curveRadix).negate) :
    blsScalarOrder • toSelectedG2Point point valid = 0 := by
  apply g2_relation_implies_prime_order_of_small_fixed_kernel
    selectedG2PsiHom (toSelectedG2Point point valid)
    selected_g2_psi_hom_six
    selected_g2_psi_hom_squared_neg_fixed_three_torsion
  change selectedG2PsiPoint (toSelectedG2Point point valid) =
    -(curveRadix • toSelectedG2Point point valid)
  rw [selected_g2_psi_point_to_selected,
    ← to_selected_g2_point_cpp_binary_mul point valid curveRadix,
    ← to_selected_g2_point_negate]
  exact (to_selected_g2_point_eq_iff _ _ _ _).mpr relation

/-! ## Conventional cyclic subgroup inclusion

The reverse inclusion needed by the public contract is constructive: the
standard generators satisfy the literal source relations, and additivity makes
the same relation hold for every integral multiple.  Combined with the two
theorems above, this gives the formally useful sandwich

`standard cyclic subgroup ⊆ accepted relation ⊆ r-torsion`.

Identifying the right-hand `r`-torsion set with the standard cyclic subgroup is
the sole remaining global elliptic-curve fact.
-/

theorem g1_integer_span_implies_relation
    {G : Type*} [AddCommGroup G]
    (sigma : G →+ G) (generator point : G)
    (generatorRelation : sigma (sigma generator) =
      -((curveRadix ^ 2) • generator))
    (span : ∃ scalar : Int, point = scalar • generator) :
    sigma (sigma point) = -((curveRadix ^ 2) • point) := by
  rcases span with ⟨scalar, rfl⟩
  cases scalar with
  | ofNat scalar =>
      have positive := g1_generator_span_implies_relation sigma generator
        (scalar • generator) generatorRelation ⟨scalar, rfl⟩
      simpa only [Int.ofNat_eq_natCast, natCast_zsmul] using positive
  | negSucc scalar =>
      have positive := g1_generator_span_implies_relation sigma generator
        ((scalar + 1) • generator) generatorRelation ⟨scalar + 1, rfl⟩
      simpa only [negSucc_zsmul, map_neg, smul_neg, neg_neg] using
        congrArg Neg.neg positive

theorem g2_integer_span_implies_relation
    {G : Type*} [AddCommGroup G]
    (psi : G →+ G) (generator point : G)
    (generatorRelation : psi generator = -(curveRadix • generator))
    (span : ∃ scalar : Int, point = scalar • generator) :
    psi point = -(curveRadix • point) := by
  rcases span with ⟨scalar, rfl⟩
  cases scalar with
  | ofNat scalar =>
      have positive := g2_generator_span_implies_relation psi generator
        (scalar • generator) generatorRelation ⟨scalar, rfl⟩
      simpa only [Int.ofNat_eq_natCast, natCast_zsmul] using positive
  | negSucc scalar =>
      have positive := g2_generator_span_implies_relation psi generator
        ((scalar + 1) • generator) generatorRelation ⟨scalar + 1, rfl⟩
      simpa only [negSucc_zsmul, map_neg, smul_neg, neg_neg] using
        congrArg Neg.neg positive

theorem selected_g1_generator_source_relation :
    (cppG1Sigma (cppG1Sigma selectedG1Generator)).equivalent
      (cppBinaryMul 4 selectedG1Generator (curveRadix ^ 2)).negate := by
  have avxBasis := selected_g1_generator_avx_basis
  rw [cpp_g1_avx_endomorphism_eq_neg_sigma_squared] at avxBasis
  have negated := homogeneous_equivalent_negate _ _ avxBasis
  simpa [homogeneous_negate_negate, curveRadixSquared] using negated

theorem selected_g2_generator_source_relation :
    (cppG2Psi selectedG2Generator).equivalent
      (cppBinaryMul selectedG2CurveB selectedG2Generator curveRadix).negate := by
  have avxBasis := selected_g2_generator_avx_basis
  rw [cpp_g2_avx_minus_psi_eq_neg_psi] at avxBasis
  have negated := homogeneous_equivalent_negate _ _ avxBasis
  simpa [homogeneous_negate_negate] using negated

theorem selected_g1_generator_relation_group :
    selectedG1SigmaHom
        (selectedG1SigmaHom selectedG1GeneratorPoint) =
      -((curveRadix ^ 2) • selectedG1GeneratorPoint) := by
  unfold selectedG1GeneratorPoint
  change selectedG1SigmaPoint
      (selectedG1SigmaPoint
        (toSelectedG1Point selectedG1Generator selected_g1_generator_valid)) =
    -((curveRadix ^ 2) •
      toSelectedG1Point selectedG1Generator selected_g1_generator_valid)
  rw [selected_g1_sigma_point_to_selected,
    selected_g1_sigma_point_to_selected,
    ← to_selected_g1_point_cpp_binary_mul,
    ← to_selected_g1_point_negate]
  exact (to_selected_g1_point_eq_iff _ _ _ _).mpr
    selected_g1_generator_source_relation

theorem selected_g2_generator_relation_group :
    selectedG2PsiHom selectedG2GeneratorPoint =
      -(curveRadix • selectedG2GeneratorPoint) := by
  unfold selectedG2GeneratorPoint
  change selectedG2PsiPoint
      (toSelectedG2Point selectedG2Generator selected_g2_generator_valid) =
    -(curveRadix •
      toSelectedG2Point selectedG2Generator selected_g2_generator_valid)
  rw [selected_g2_psi_point_to_selected,
    ← to_selected_g2_point_cpp_binary_mul,
    ← to_selected_g2_point_negate]
  exact (to_selected_g2_point_eq_iff _ _ _ _).mpr
    selected_g2_generator_source_relation

theorem selected_g1_cyclic_membership_implies_source_relation
    (point : HomogeneousPoint SelectedCurveFp)
    (valid : point.projectivelyValid 4)
    (membership : toSelectedG1Point point valid ∈
      AddSubgroup.zmultiples selectedG1GeneratorPoint) :
    (cppG1Sigma (cppG1Sigma point)).equivalent
      (cppBinaryMul 4 point (curveRadix ^ 2)).negate := by
  have groupRelation := g1_integer_span_implies_relation
    selectedG1SigmaHom selectedG1GeneratorPoint
      (toSelectedG1Point point valid)
    selected_g1_generator_relation_group (by
      rcases AddSubgroup.mem_zmultiples_iff.mp membership with ⟨scalar, span⟩
      exact ⟨scalar, span.symm⟩)
  change selectedG1SigmaPoint
      (selectedG1SigmaPoint (toSelectedG1Point point valid)) =
    -((curveRadix ^ 2) • toSelectedG1Point point valid) at groupRelation
  rw [selected_g1_sigma_point_to_selected,
    selected_g1_sigma_point_to_selected,
    ← to_selected_g1_point_cpp_binary_mul,
    ← to_selected_g1_point_negate] at groupRelation
  exact (to_selected_g1_point_eq_iff _ _ _ _).mp groupRelation

theorem selected_g2_cyclic_membership_implies_source_relation
    (point : HomogeneousPoint SelectedMillerFp2)
    (valid : point.projectivelyValid selectedG2CurveB)
    (membership : toSelectedG2Point point valid ∈
      AddSubgroup.zmultiples selectedG2GeneratorPoint) :
    (cppG2Psi point).equivalent
      (cppBinaryMul selectedG2CurveB point curveRadix).negate := by
  have groupRelation := g2_integer_span_implies_relation
    selectedG2PsiHom selectedG2GeneratorPoint
      (toSelectedG2Point point valid)
    selected_g2_generator_relation_group (by
      rcases AddSubgroup.mem_zmultiples_iff.mp membership with ⟨scalar, span⟩
      exact ⟨scalar, span.symm⟩)
  change selectedG2PsiPoint (toSelectedG2Point point valid) =
    -(curveRadix • toSelectedG2Point point valid) at groupRelation
  rw [selected_g2_psi_point_to_selected,
    ← to_selected_g2_point_cpp_binary_mul,
    ← to_selected_g2_point_negate] at groupRelation
  exact (to_selected_g2_point_eq_iff _ _ _ _).mp groupRelation

/-! ## Validator-accepted subgroup carriers

These kernels are the exact additive domains recognized by the source
relations after the public on-curve check.  They are intentionally distinct
from the conventional generator spans until the final global eigenspace-size
theorem is available.
-/

noncomputable def selectedG1RelationHom :
    selectedG1Weierstrass.toProjective.Point →+
      selectedG1Weierstrass.toProjective.Point :=
  selectedG1SigmaHom.comp selectedG1SigmaHom +
    nsmulAddMonoidHom (curveRadix ^ 2)

noncomputable def selectedG2RelationHom :
    selectedG2Weierstrass.toProjective.Point →+
      selectedG2Weierstrass.toProjective.Point :=
  selectedG2PsiHom + nsmulAddMonoidHom curveRadix

noncomputable abbrev SelectedG1AcceptedSubgroup :=
  selectedG1RelationHom.ker

noncomputable abbrev SelectedG2AcceptedSubgroup :=
  selectedG2RelationHom.ker

theorem selected_g1_mem_accepted_iff_relation
    (point : selectedG1Weierstrass.toProjective.Point) :
    point ∈ SelectedG1AcceptedSubgroup ↔
      selectedG1SigmaHom (selectedG1SigmaHom point) =
        -((curveRadix ^ 2) • point) := by
  rw [AddMonoidHom.mem_ker]
  change selectedG1SigmaHom (selectedG1SigmaHom point) +
      (curveRadix ^ 2) • point = 0 ↔ _
  exact add_eq_zero_iff_eq_neg

theorem selected_g2_mem_accepted_iff_relation
    (point : selectedG2Weierstrass.toProjective.Point) :
    point ∈ SelectedG2AcceptedSubgroup ↔
      selectedG2PsiHom point = -(curveRadix • point) := by
  rw [AddMonoidHom.mem_ker]
  change selectedG2PsiHom point + curveRadix • point = 0 ↔ _
  exact add_eq_zero_iff_eq_neg

theorem selected_g1_source_relation_iff_accepted
    (point : HomogeneousPoint SelectedCurveFp)
    (valid : point.projectivelyValid 4) :
    (cppG1Sigma (cppG1Sigma point)).equivalent
        (cppBinaryMul 4 point (curveRadix ^ 2)).negate ↔
      toSelectedG1Point point valid ∈ SelectedG1AcceptedSubgroup := by
  rw [selected_g1_mem_accepted_iff_relation]
  change _ ↔ selectedG1SigmaPoint
      (selectedG1SigmaPoint (toSelectedG1Point point valid)) =
    -((curveRadix ^ 2) • toSelectedG1Point point valid)
  rw [selected_g1_sigma_point_to_selected,
    selected_g1_sigma_point_to_selected,
    ← to_selected_g1_point_cpp_binary_mul,
    ← to_selected_g1_point_negate,
    to_selected_g1_point_eq_iff]

theorem selected_g2_source_relation_iff_accepted
    (point : HomogeneousPoint SelectedMillerFp2)
    (valid : point.projectivelyValid selectedG2CurveB) :
    (cppG2Psi point).equivalent
        (cppBinaryMul selectedG2CurveB point curveRadix).negate ↔
      toSelectedG2Point point valid ∈ SelectedG2AcceptedSubgroup := by
  rw [selected_g2_mem_accepted_iff_relation]
  change _ ↔ selectedG2PsiPoint (toSelectedG2Point point valid) =
    -(curveRadix • toSelectedG2Point point valid)
  rw [selected_g2_psi_point_to_selected,
    ← to_selected_g2_point_cpp_binary_mul,
    ← to_selected_g2_point_negate,
    to_selected_g2_point_eq_iff]

theorem selected_g1_cyclic_le_accepted :
    AddSubgroup.zmultiples selectedG1GeneratorPoint ≤
      SelectedG1AcceptedSubgroup := by
  intro point membership
  rw [selected_g1_mem_accepted_iff_relation]
  exact g1_integer_span_implies_relation selectedG1SigmaHom
    selectedG1GeneratorPoint point selected_g1_generator_relation_group (by
      rcases AddSubgroup.mem_zmultiples_iff.mp membership with ⟨scalar, span⟩
      exact ⟨scalar, span.symm⟩)

theorem selected_g2_cyclic_le_accepted :
    AddSubgroup.zmultiples selectedG2GeneratorPoint ≤
      SelectedG2AcceptedSubgroup := by
  intro point membership
  rw [selected_g2_mem_accepted_iff_relation]
  exact g2_integer_span_implies_relation selectedG2PsiHom
    selectedG2GeneratorPoint point selected_g2_generator_relation_group (by
      rcases AddSubgroup.mem_zmultiples_iff.mp membership with ⟨scalar, span⟩
      exact ⟨scalar, span.symm⟩)

/-!
The remaining global identification obligation is now exactly a cardinality
upper bound. The reverse inequality is already forced by the concrete
generator's exact order and the inclusions above.
-/

theorem selected_g1_scalar_order_le_accepted_card :
    blsScalarOrder ≤ Nat.card SelectedG1AcceptedSubgroup := by
  rw [← selected_g1_generator_span_natCard]
  apply Nat.card_le_card_of_injective
    (fun point : AddSubgroup.zmultiples selectedG1GeneratorPoint =>
      (⟨point, selected_g1_cyclic_le_accepted point.property⟩ :
        SelectedG1AcceptedSubgroup))
  intro left right equal
  apply Subtype.ext
  exact congrArg
    (fun point : SelectedG1AcceptedSubgroup =>
      (point : selectedG1Weierstrass.toProjective.Point)) equal

theorem selected_g2_scalar_order_le_accepted_card :
    blsScalarOrder ≤ Nat.card SelectedG2AcceptedSubgroup := by
  rw [← selected_g2_generator_span_natCard]
  apply Nat.card_le_card_of_injective
    (fun point : AddSubgroup.zmultiples selectedG2GeneratorPoint =>
      (⟨point, selected_g2_cyclic_le_accepted point.property⟩ :
        SelectedG2AcceptedSubgroup))
  intro left right equal
  apply Subtype.ext
  exact congrArg
    (fun point : SelectedG2AcceptedSubgroup =>
      (point : selectedG2Weierstrass.toProjective.Point)) equal

theorem selected_g1_cyclic_eq_accepted_of_card_le
    (cardBound : Nat.card SelectedG1AcceptedSubgroup ≤ blsScalarOrder) :
    AddSubgroup.zmultiples selectedG1GeneratorPoint =
      SelectedG1AcceptedSubgroup := by
  apply AddSubgroup.eq_of_le_of_card_ge selected_g1_cyclic_le_accepted
  rw [selected_g1_generator_span_natCard]
  exact cardBound

theorem selected_g2_cyclic_eq_accepted_of_card_le
    (cardBound : Nat.card SelectedG2AcceptedSubgroup ≤ blsScalarOrder) :
    AddSubgroup.zmultiples selectedG2GeneratorPoint =
      SelectedG2AcceptedSubgroup := by
  apply AddSubgroup.eq_of_le_of_card_ge selected_g2_cyclic_le_accepted
  rw [selected_g2_generator_span_natCard]
  exact cardBound

theorem selected_g1_cyclic_eq_accepted_iff_card_le :
    AddSubgroup.zmultiples selectedG1GeneratorPoint =
        SelectedG1AcceptedSubgroup ↔
      Nat.card SelectedG1AcceptedSubgroup ≤ blsScalarOrder := by
  constructor
  · intro equalSubgroups
    rw [← equalSubgroups, selected_g1_generator_span_natCard]
  · exact selected_g1_cyclic_eq_accepted_of_card_le

theorem selected_g2_cyclic_eq_accepted_iff_card_le :
    AddSubgroup.zmultiples selectedG2GeneratorPoint =
        SelectedG2AcceptedSubgroup ↔
      Nat.card SelectedG2AcceptedSubgroup ≤ blsScalarOrder := by
  constructor
  · intro equalSubgroups
    rw [← equalSubgroups, selected_g2_generator_span_natCard]
  · exact selected_g2_cyclic_eq_accepted_of_card_le

noncomputable instance selectedG1AcceptedSubgroupModule :
    Module CurveExponent SelectedG1AcceptedSubgroup :=
  AddCommGroup.zmodModule (n := blsScalarOrder) (by
    intro point
    apply Subtype.ext
    change blsScalarOrder •
      (point : selectedG1Weierstrass.toProjective.Point) = 0
    exact g1_relation_implies_prime_order_of_small_fixed_kernel
      selectedG1SigmaHom point selected_g1_sigma_hom_cube
      selected_g1_sigma_hom_squared_fixed_three_torsion
      ((selected_g1_mem_accepted_iff_relation point).mp point.property))

noncomputable instance selectedG2AcceptedSubgroupModule :
    Module CurveExponent SelectedG2AcceptedSubgroup :=
  AddCommGroup.zmodModule (n := blsScalarOrder) (by
    intro point
    apply Subtype.ext
    change blsScalarOrder •
      (point : selectedG2Weierstrass.toProjective.Point) = 0
    exact g2_relation_implies_prime_order_of_small_fixed_kernel
      selectedG2PsiHom point selected_g2_psi_hom_six
      selected_g2_psi_hom_squared_neg_fixed_three_torsion
      ((selected_g2_mem_accepted_iff_relation point).mp point.property))

theorem selected_g1_accepted_unsigned_basis
    (point : SelectedG1AcceptedSubgroup) :
    -selectedG1SigmaHom (selectedG1SigmaHom
        (point : selectedG1Weierstrass.toProjective.Point)) =
      (curveRadix ^ 2) •
        (point : selectedG1Weierstrass.toProjective.Point) := by
  have relation :=
    (selected_g1_mem_accepted_iff_relation point).mp point.property
  rw [relation]
  exact neg_neg _

theorem selected_g2_accepted_unsigned_basis
    (point : SelectedG2AcceptedSubgroup) :
    -selectedG2PsiHom
        (point : selectedG2Weierstrass.toProjective.Point) =
      curveRadix •
        (point : selectedG2Weierstrass.toProjective.Point) := by
  have relation :=
    (selected_g2_mem_accepted_iff_relation point).mp point.property
  rw [relation]
  exact neg_neg _

theorem selected_g1_accepted_mult_correct
    (point : SelectedG1AcceptedSubgroup) (scalar nbits : Nat) :
    selectedG1Mult point scalar nbits =
      (curveMaskedScalar scalar nbits : CurveExponent) • point :=
  selected_g1_mult_correct point scalar nbits

theorem selected_g1_accepted_mult_vartime_correct
    (point : SelectedG1AcceptedSubgroup) (scalar nbits : Nat) :
    selectedG1MultVartime point scalar nbits =
      (curveMaskedScalar scalar nbits : CurveExponent) • point :=
  selected_g1_mult_vartime_correct point scalar nbits

theorem selected_g2_accepted_mult_correct
    (point : SelectedG2AcceptedSubgroup) (scalar nbits : Nat) :
    selectedG2Mult point scalar nbits =
      (curveMaskedScalar scalar nbits : CurveExponent) • point :=
  selected_g2_mult_correct point scalar nbits

theorem selected_g2_accepted_mult_vartime_correct
    (point : SelectedG2AcceptedSubgroup) (scalar nbits : Nat) :
    selectedG2MultVartime point scalar nbits =
      (curveMaskedScalar scalar nbits : CurveExponent) • point :=
  selected_g2_mult_vartime_correct point scalar nbits

end CatBLST.AVX.MatrixNoK
