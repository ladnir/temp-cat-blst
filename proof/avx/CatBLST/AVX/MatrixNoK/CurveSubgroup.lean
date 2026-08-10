/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.CurveEndomorphism
import CatBLST.AVX.MatrixNoK.CurveScalarSplit
import Mathlib.GroupTheory.OrderOfElement
import Mathlib.Tactic

set_option maxHeartbeats 800000

/-!
# Reduction of the fast subgroup relation to explicit group-order facts

The source predicates use endomorphism relations instead of multiplying every
candidate by the 255-bit scalar order.  This file separates the inexpensive
algebraic reduction from the genuinely large BLS12-381 point-count theorem.

For G1, the relation and `sigma^3 = 1` imply that `(z^6 + 1)P = 0`.
Because `z^6 + 1 = r(z^2 + 1)` and the G1 cofactor is coprime to
`z^2 + 1`, the known full curve order `h1*r` then implies `[r]P = 0`.
The analogous G2 proof uses `psi^6 = -1`, the same annihilator, and the G2
cofactor.  The theorems below make the full-order and prime-subgroup spanning
facts explicit premises; neither is silently inferred from coordinate formulas.

The concrete G2 cofactor and characteristic-polynomial checks also instantiate
Proposition 5 of El Housni, Guillevic, and Piellard, ePrint 2022/352, which
repairs the circular argument in the original Scott G2 proof.  Their published
BLS12 theorem supplies the external curve-family fact; Lean checks the exact
BLS12-381 seed, cofactor polynomial, and coprimality hypothesis used here.
-/

namespace CatBLST.AVX.MatrixNoK

/-- BLS12-381 G1 cofactor. -/
def g1Cofactor : Nat :=
  0x396c8c005555e1568c00aaab0000aaab

/-- BLS12-381 G2 cofactor. -/
def g2Cofactor : Nat :=
  0x5d543a95414e7f1091d50792876a202cd91de4547085abaa68a205b2e5a7ddfa628f1cb4d9e82ef21537e293a6691ae1616ec6e786f0c70cf1c38e31c7238e5

/-- `chi(u)` for the G2 Frobenius characteristic polynomial
`chi(X) = X^2 - (u + 1)X + p` at the negative BLS12-381 seed `u`. -/
def g2PsiCharacteristicAtSeed : Nat :=
  Generated.target_modulus + curveRadix

theorem target_modulus_as_bls12_polynomial :
    3 * (Generated.target_modulus + curveRadix) =
      (curveRadix + 1) ^ 2 * blsScalarOrder := by
  native_decide

theorem g2_psi_characteristic_at_seed :
    blsZ ^ 2 - (blsZ + 1) * blsZ + (Generated.target_modulus : Int) =
      (g2PsiCharacteristicAtSeed : Int) := by
  native_decide

theorem g2_cofactor_as_bls12_polynomial :
    (9 : Int) * g2Cofactor =
      blsZ ^ 8 - 4 * blsZ ^ 7 + 5 * blsZ ^ 6 - 4 * blsZ ^ 4 +
        6 * blsZ ^ 3 - 4 * blsZ ^ 2 - 4 * blsZ + 13 := by
  native_decide

/-- The concrete hypothesis that excludes every cofactor eigenspace from the
accepted `psi(P) = [u]P` relation in ePrint 2022/352, Proposition 5. -/
theorem g2_scott_characteristic_coprime_cofactor :
    Nat.Coprime g2PsiCharacteristicAtSeed g2Cofactor := by
  native_decide

theorem scalar_order_as_bls_polynomial :
    blsScalarOrder = curveRadix ^ 4 - curveRadix ^ 2 + 1 := by
  native_decide

theorem curve_radix_six_factorization :
    curveRadix ^ 6 + 1 =
      blsScalarOrder * (curveRadix ^ 2 + 1) := by
  native_decide

theorem g1_cofactor_coprime_radix_squared_plus_one :
    Nat.Coprime g1Cofactor (curveRadix ^ 2 + 1) := by
  native_decide

theorem g1_relation_full_order_gcd :
    Nat.gcd (curveRadix ^ 6 + 1) (g1Cofactor * blsScalarOrder) =
      blsScalarOrder := by
  native_decide

theorem g2_cofactor_coprime_radix_squared_plus_one :
    Nat.Coprime g2Cofactor (curveRadix ^ 2 + 1) := by
  native_decide

theorem g2_relation_full_order_gcd :
    Nat.gcd (curveRadix ^ 6 + 1) (g2Cofactor * blsScalarOrder) =
      blsScalarOrder := by
  native_decide

theorem g1_relation_annihilates_radix_six_plus_one
    {G : Type*} [AddCommGroup G]
    (sigma : G →+ G) (point : G)
    (sigmaCube : ∀ value, sigma (sigma (sigma value)) = value)
    (relation : sigma (sigma point) =
      -((curveRadix ^ 2) • point)) :
    (curveRadix ^ 6 + 1) • point = 0 := by
  let a := curveRadix ^ 2
  have sigmaFourth :
      sigma (sigma (sigma (sigma point))) = (a * a) • point := by
    calc
      sigma (sigma (sigma (sigma point))) =
          sigma (sigma (-((a) • point))) := by
        rw [show a = curveRadix ^ 2 by rfl, relation]
      _ = -(a • sigma (sigma point)) := by simp
      _ = -(a • -((a) • point)) := by
        rw [show a = curveRadix ^ 2 by rfl, relation]
      _ = (a * a) • point := by simp [mul_nsmul]
  have sigmaSixth :
      sigma (sigma (sigma (sigma (sigma (sigma point))))) =
        -((a ^ 3) • point) := by
    calc
      sigma (sigma (sigma (sigma (sigma (sigma point))))) =
          sigma (sigma ((a * a) • point)) := by rw [sigmaFourth]
      _ = (a * a) • sigma (sigma point) := by simp
      _ = (a * a) • (-((a) • point)) := by
        rw [show a = curveRadix ^ 2 by rfl, relation]
      _ = -((a ^ 3) • point) := by
        simp [mul_nsmul, pow_succ, mul_assoc]
  have sigmaSixthIdentity :
      sigma (sigma (sigma (sigma (sigma (sigma point))))) = point := by
    calc
      sigma (sigma (sigma (sigma (sigma (sigma point))))) =
          sigma (sigma (sigma point)) := by
        rw [sigmaCube]
      _ = point := sigmaCube point
  have pointRelation : point = -((a ^ 3) • point) := by
    exact sigmaSixthIdentity.symm.trans sigmaSixth
  rw [show curveRadix ^ 6 = a ^ 3 by simp [a]; ring]
  calc
    (a ^ 3 + 1) • point = a ^ 3 • point + point := by
      rw [add_nsmul, one_nsmul]
    _ = a ^ 3 • point + -(a ^ 3 • point) :=
      congrArg (fun value => a ^ 3 • point + value) pointRelation
    _ = 0 := add_neg_cancel _

theorem g1_relation_implies_prime_order
    {G : Type*} [AddCommGroup G] [Finite G]
    (sigma : G →+ G) (point : G)
    (fullOrder : Nat.card G = g1Cofactor * blsScalarOrder)
    (sigmaCube : ∀ value, sigma (sigma (sigma value)) = value)
    (relation : sigma (sigma point) =
      -((curveRadix ^ 2) • point)) :
    blsScalarOrder • point = 0 := by
  have relationAnnihilator :=
    g1_relation_annihilates_radix_six_plus_one sigma point sigmaCube relation
  have orderDividesRelation :
      addOrderOf point ∣ curveRadix ^ 6 + 1 :=
    addOrderOf_dvd_iff_nsmul_eq_zero.mpr relationAnnihilator
  have orderDividesFull :
      addOrderOf point ∣ g1Cofactor * blsScalarOrder := by
    rw [← fullOrder]
    exact addOrderOf_dvd_natCard point
  have orderDividesGcd :
      addOrderOf point ∣
        Nat.gcd (curveRadix ^ 6 + 1)
          (g1Cofactor * blsScalarOrder) :=
    Nat.dvd_gcd orderDividesRelation orderDividesFull
  rw [g1_relation_full_order_gcd] at orderDividesGcd
  exact addOrderOf_dvd_iff_nsmul_eq_zero.mp orderDividesGcd

theorem g2_relation_annihilates_radix_six_plus_one
    {G : Type*} [AddCommGroup G]
    (psi : G →+ G) (point : G)
    (psiSix : ∀ value,
      psi (psi (psi (psi (psi (psi value))))) = -value)
    (relation : psi point = -(curveRadix • point)) :
    (curveRadix ^ 6 + 1) • point = 0 := by
  let a := curveRadix ^ 2
  have psiSquared : psi (psi point) = a • point := by
    calc
      psi (psi point) = psi (-(curveRadix • point)) := by rw [relation]
      _ = -(curveRadix • psi point) := by simp
      _ = -(curveRadix • -(curveRadix • point)) := by rw [relation]
      _ = a • point := by simp [a, mul_nsmul, pow_two]
  have psiFourth :
      psi (psi (psi (psi point))) = (a * a) • point := by
    calc
      psi (psi (psi (psi point))) = psi (psi (a • point)) := by
        rw [psiSquared]
      _ = a • psi (psi point) := by simp
      _ = a • (a • point) := by rw [psiSquared]
      _ = (a * a) • point := by simp [mul_nsmul]
  have psiSixth :
      psi (psi (psi (psi (psi (psi point))))) = (a ^ 3) • point := by
    calc
      psi (psi (psi (psi (psi (psi point))))) =
          psi (psi ((a * a) • point)) := by rw [psiFourth]
      _ = (a * a) • psi (psi point) := by simp
      _ = (a * a) • (a • point) := by rw [psiSquared]
      _ = (a ^ 3) • point := by simp [mul_nsmul, pow_succ, mul_assoc]
  have pointRelation : (a ^ 3) • point = -point := by
    exact psiSixth.symm.trans (psiSix point)
  rw [show curveRadix ^ 6 = a ^ 3 by simp [a]; ring]
  calc
    (a ^ 3 + 1) • point = a ^ 3 • point + point := by
      rw [add_nsmul, one_nsmul]
    _ = -point + point :=
      congrArg (fun value => value + point) pointRelation
    _ = 0 := neg_add_cancel point

theorem g2_relation_implies_prime_order
    {G : Type*} [AddCommGroup G] [Finite G]
    (psi : G →+ G) (point : G)
    (fullOrder : Nat.card G = g2Cofactor * blsScalarOrder)
    (psiSix : ∀ value,
      psi (psi (psi (psi (psi (psi value))))) = -value)
    (relation : psi point = -(curveRadix • point)) :
    blsScalarOrder • point = 0 := by
  have relationAnnihilator :=
    g2_relation_annihilates_radix_six_plus_one psi point psiSix relation
  have orderDividesRelation :
      addOrderOf point ∣ curveRadix ^ 6 + 1 :=
    addOrderOf_dvd_iff_nsmul_eq_zero.mpr relationAnnihilator
  have orderDividesFull :
      addOrderOf point ∣ g2Cofactor * blsScalarOrder := by
    rw [← fullOrder]
    exact addOrderOf_dvd_natCard point
  have orderDividesGcd :
      addOrderOf point ∣
        Nat.gcd (curveRadix ^ 6 + 1)
          (g2Cofactor * blsScalarOrder) :=
    Nat.dvd_gcd orderDividesRelation orderDividesFull
  rw [g2_relation_full_order_gcd] at orderDividesGcd
  exact addOrderOf_dvd_iff_nsmul_eq_zero.mp orderDividesGcd

theorem g1_generator_span_implies_relation
    {G : Type*} [AddCommGroup G]
    (sigma : G →+ G) (generator point : G)
    (generatorRelation : sigma (sigma generator) =
      -((curveRadix ^ 2) • generator))
    (span : ∃ scalar : Nat, point = scalar • generator) :
    sigma (sigma point) = -((curveRadix ^ 2) • point) := by
  rcases span with ⟨scalar, rfl⟩
  simp only [map_nsmul]
  rw [generatorRelation]
  simp only [smul_neg, ← mul_nsmul]
  rw [Nat.mul_comm scalar (curveRadix ^ 2)]

theorem g2_generator_span_implies_relation
    {G : Type*} [AddCommGroup G]
    (psi : G →+ G) (generator point : G)
    (generatorRelation : psi generator = -(curveRadix • generator))
    (span : ∃ scalar : Nat, point = scalar • generator) :
    psi point = -(curveRadix • point) := by
  rcases span with ⟨scalar, rfl⟩
  simp only [map_nsmul]
  rw [generatorRelation]
  simp only [smul_neg, ← mul_nsmul]
  rw [Nat.mul_comm scalar curveRadix]

/-!
These equivalences isolate the two external facts still required for a concrete
BLS12-381 instantiation: the full curve cardinality and that every point killed
by `r` is generated by the selected prime-order generator.  The forward
directions need only the cardinality; the reverse directions need only the
generator relation and spanning theorem.
-/
theorem g1_relation_iff_prime_subgroup
    {G : Type*} [AddCommGroup G] [Finite G]
    (sigma : G →+ G) (generator point : G)
    (fullOrder : Nat.card G = g1Cofactor * blsScalarOrder)
    (sigmaCube : ∀ value, sigma (sigma (sigma value)) = value)
    (generatorRelation : sigma (sigma generator) =
      -((curveRadix ^ 2) • generator))
    (primeSubgroupSpan : ∀ value, blsScalarOrder • value = 0 →
      ∃ scalar : Nat, value = scalar • generator) :
    sigma (sigma point) = -((curveRadix ^ 2) • point) ↔
      blsScalarOrder • point = 0 := by
  constructor
  · exact g1_relation_implies_prime_order sigma point fullOrder sigmaCube
  · intro primeOrder
    exact g1_generator_span_implies_relation sigma generator point
      generatorRelation (primeSubgroupSpan point primeOrder)

theorem g2_relation_iff_prime_subgroup
    {G : Type*} [AddCommGroup G] [Finite G]
    (psi : G →+ G) (generator point : G)
    (fullOrder : Nat.card G = g2Cofactor * blsScalarOrder)
    (psiSix : ∀ value,
      psi (psi (psi (psi (psi (psi value))))) = -value)
    (generatorRelation : psi generator = -(curveRadix • generator))
    (primeSubgroupSpan : ∀ value, blsScalarOrder • value = 0 →
      ∃ scalar : Nat, value = scalar • generator) :
    psi point = -(curveRadix • point) ↔ blsScalarOrder • point = 0 := by
  constructor
  · exact g2_relation_implies_prime_order psi point fullOrder psiSix
  · intro primeOrder
    exact g2_generator_span_implies_relation psi generator point
      generatorRelation (primeSubgroupSpan point primeOrder)

end CatBLST.AVX.MatrixNoK
