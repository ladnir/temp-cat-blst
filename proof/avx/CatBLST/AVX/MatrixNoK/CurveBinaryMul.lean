/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.CurveGroupAdd
import Mathlib.Data.Nat.BinaryRec

/-!
# Proof-only binary scalar multiplication

This is a logarithmic evaluator used only to certify concrete BLS12-381
generator relations.  It is deliberately built from `cppCompleteAdd`, so its
coordinate computation can be reduced by the kernel while its semantic proof
lands in Mathlib's projective group law.  It is not production code and does
not define or recommend a new runtime scalar-multiplication path.
-/

namespace CatBLST.AVX.MatrixNoK

def cppBinaryMul [CommRing K] (b : K) (point : HomogeneousPoint K) :
    Nat → HomogeneousPoint K :=
  Nat.binaryRec' cppProjectiveZero fun bit _ _ accumulated =>
    let doubled := cppCompleteAdd b accumulated accumulated
    if bit then cppCompleteAdd b doubled point else doubled

@[simp] theorem cpp_binary_mul_zero [CommRing K] (b : K)
    (point : HomogeneousPoint K) :
    cppBinaryMul b point 0 = cppProjectiveZero := by
  simp [cppBinaryMul]

theorem cpp_binary_mul_bit [CommRing K] (b : K)
    (point : HomogeneousPoint K) (bit : Bool) (n : Nat)
    (bitValid : n = 0 → bit = true) :
    cppBinaryMul b point (Nat.bit bit n) =
      let doubled := cppCompleteAdd b (cppBinaryMul b point n)
        (cppBinaryMul b point n)
      if bit then cppCompleteAdd b doubled point else doubled := by
  simp only [cppBinaryMul, Nat.binaryRec'_eq bit n bitValid]

theorem cpp_binary_mul_projectively_valid [Field K] [DecidableEq K]
    (b : K) (noTwoTorsion : NoAffineTwoTorsion b) (twoNonzero : (2 : K) ≠ 0)
    (point : HomogeneousPoint K) (pointValid : point.projectivelyValid b) :
    ∀ n, (cppBinaryMul b point n).projectivelyValid b := by
  intro n
  induction n using Nat.binaryRec' with
  | zero => exact projective_zero_valid b
  | bit digit n digitValid ih =>
      rw [cpp_binary_mul_bit b point digit n digitValid]
      have doubledValid := cpp_complete_add_projectively_valid b
        noTwoTorsion twoNonzero _ _ ih ih
      cases digit
      · exact doubledValid
      · exact cpp_complete_add_projectively_valid b
          noTwoTorsion twoNonzero _ _ doubledValid pointValid

theorem to_projective_point_cpp_binary_mul [Field K] [DecidableEq K]
    (b : K) [elliptic : (shortWeierstrass b).IsElliptic]
    (noTwoTorsion : NoAffineTwoTorsion b) (twoNonzero : (2 : K) ≠ 0)
    (point : HomogeneousPoint K) (pointValid : point.projectivelyValid b) :
    ∀ n,
      toProjectivePoint b (cppBinaryMul b point n)
          (cpp_binary_mul_projectively_valid b noTwoTorsion twoNonzero
            point pointValid n) =
        n • toProjectivePoint b point pointValid := by
  intro n
  induction n using Nat.binaryRec' with
  | zero => exact to_projective_point_zero b
  | bit digit n digitValid ih =>
      have accumulatedValid := cpp_binary_mul_projectively_valid b
        noTwoTorsion twoNonzero point pointValid n
      have doubledValid := cpp_complete_add_projectively_valid b
        noTwoTorsion twoNonzero _ _ accumulatedValid accumulatedValid
      have doubledSemantics := projective_complete_add_eq_add b
        noTwoTorsion twoNonzero _ _ accumulatedValid accumulatedValid
      cases digit
      · simp only [cpp_binary_mul_bit b point false n digitValid,
          Bool.false_eq_true, if_false]
        calc
          toProjectivePoint b
              (cppCompleteAdd b (cppBinaryMul b point n)
                (cppBinaryMul b point n)) doubledValid =
              toProjectivePoint b (cppBinaryMul b point n) accumulatedValid +
                toProjectivePoint b (cppBinaryMul b point n) accumulatedValid :=
            doubledSemantics
          _ = n • toProjectivePoint b point pointValid +
                n • toProjectivePoint b point pointValid := by rw [ih]
          _ = (2 * n) • toProjectivePoint b point pointValid := by
            rw [mul_nsmul, two_nsmul, nsmul_add]
          _ = Nat.bit false n • toProjectivePoint b point pointValid := by rfl
      · simp only [cpp_binary_mul_bit b point true n digitValid, if_true]
        have sumValid := cpp_complete_add_projectively_valid b
          noTwoTorsion twoNonzero _ _ doubledValid pointValid
        calc
          toProjectivePoint b
              (cppCompleteAdd b
                (cppCompleteAdd b (cppBinaryMul b point n)
                  (cppBinaryMul b point n)) point) sumValid =
              toProjectivePoint b
                  (cppCompleteAdd b (cppBinaryMul b point n)
                    (cppBinaryMul b point n)) doubledValid +
                toProjectivePoint b point pointValid :=
            projective_complete_add_eq_add b noTwoTorsion twoNonzero
              _ _ doubledValid pointValid
          _ = (toProjectivePoint b (cppBinaryMul b point n) accumulatedValid +
                toProjectivePoint b (cppBinaryMul b point n) accumulatedValid) +
                toProjectivePoint b point pointValid := by rw [doubledSemantics]
          _ = (n • toProjectivePoint b point pointValid +
                n • toProjectivePoint b point pointValid) +
                toProjectivePoint b point pointValid := by rw [ih]
          _ = (2 * n + 1) • toProjectivePoint b point pointValid := by
            rw [add_nsmul, mul_nsmul, two_nsmul, nsmul_add, one_nsmul]
          _ = Nat.bit true n • toProjectivePoint b point pointValid := by rfl

theorem to_selected_g1_point_cpp_binary_mul
    (point : HomogeneousPoint SelectedCurveFp)
    (pointValid : point.projectivelyValid 4) :
    ∀ n,
      toSelectedG1Point (cppBinaryMul 4 point n)
          (cpp_binary_mul_projectively_valid 4
            selected_g1_no_affine_two_torsion selected_g1_two_ne_zero
            point pointValid n) =
        n • toSelectedG1Point point pointValid := by
  intro n
  induction n using Nat.binaryRec' with
  | zero => exact to_selected_g1_point_zero
  | bit digit n digitValid ih =>
      have accumulatedValid := cpp_binary_mul_projectively_valid 4
        selected_g1_no_affine_two_torsion selected_g1_two_ne_zero
        point pointValid n
      have doubledValid := cpp_complete_add_projectively_valid 4
        selected_g1_no_affine_two_torsion selected_g1_two_ne_zero
        _ _ accumulatedValid accumulatedValid
      have doubledSemantics := selected_g1_complete_add_eq_add
        _ _ accumulatedValid accumulatedValid
      cases digit
      · simp only [cpp_binary_mul_bit 4 point false n digitValid,
          Bool.false_eq_true, if_false]
        calc
          toSelectedG1Point
              (cppCompleteAdd 4 (cppBinaryMul 4 point n)
                (cppBinaryMul 4 point n)) doubledValid =
              toSelectedG1Point (cppBinaryMul 4 point n) accumulatedValid +
                toSelectedG1Point (cppBinaryMul 4 point n) accumulatedValid :=
            doubledSemantics
          _ = n • toSelectedG1Point point pointValid +
                n • toSelectedG1Point point pointValid := by rw [ih]
          _ = (2 * n) • toSelectedG1Point point pointValid := by
            rw [mul_nsmul, two_nsmul, nsmul_add]
          _ = Nat.bit false n • toSelectedG1Point point pointValid := by rfl
      · simp only [cpp_binary_mul_bit 4 point true n digitValid, if_true]
        have sumValid := cpp_complete_add_projectively_valid 4
          selected_g1_no_affine_two_torsion selected_g1_two_ne_zero
          _ _ doubledValid pointValid
        calc
          toSelectedG1Point
              (cppCompleteAdd 4
                (cppCompleteAdd 4 (cppBinaryMul 4 point n)
                  (cppBinaryMul 4 point n)) point) sumValid =
              toSelectedG1Point
                  (cppCompleteAdd 4 (cppBinaryMul 4 point n)
                    (cppBinaryMul 4 point n)) doubledValid +
                toSelectedG1Point point pointValid :=
            selected_g1_complete_add_eq_add _ _ doubledValid pointValid
          _ = (toSelectedG1Point (cppBinaryMul 4 point n) accumulatedValid +
                toSelectedG1Point (cppBinaryMul 4 point n) accumulatedValid) +
                toSelectedG1Point point pointValid := by rw [doubledSemantics]
          _ = (n • toSelectedG1Point point pointValid +
                n • toSelectedG1Point point pointValid) +
                toSelectedG1Point point pointValid := by rw [ih]
          _ = (2 * n + 1) • toSelectedG1Point point pointValid := by
            rw [add_nsmul, mul_nsmul, two_nsmul, nsmul_add, one_nsmul]
          _ = Nat.bit true n • toSelectedG1Point point pointValid := by rfl

theorem to_selected_g2_point_cpp_binary_mul
    (point : HomogeneousPoint SelectedMillerFp2)
    (pointValid : point.projectivelyValid selectedG2CurveB) :
    ∀ n,
      toSelectedG2Point (cppBinaryMul selectedG2CurveB point n)
          (cpp_binary_mul_projectively_valid selectedG2CurveB
            selected_g2_no_affine_two_torsion selected_g2_two_ne_zero
            point pointValid n) =
        n • toSelectedG2Point point pointValid := by
  intro n
  induction n using Nat.binaryRec' with
  | zero => exact to_selected_g2_point_zero
  | bit digit n digitValid ih =>
      have accumulatedValid := cpp_binary_mul_projectively_valid selectedG2CurveB
        selected_g2_no_affine_two_torsion selected_g2_two_ne_zero
        point pointValid n
      have doubledValid := cpp_complete_add_projectively_valid selectedG2CurveB
        selected_g2_no_affine_two_torsion selected_g2_two_ne_zero
        _ _ accumulatedValid accumulatedValid
      have doubledSemantics := selected_g2_complete_add_eq_add
        _ _ accumulatedValid accumulatedValid
      cases digit
      · simp only [cpp_binary_mul_bit selectedG2CurveB point false n digitValid,
          Bool.false_eq_true, if_false]
        calc
          toSelectedG2Point
              (cppCompleteAdd selectedG2CurveB
                (cppBinaryMul selectedG2CurveB point n)
                (cppBinaryMul selectedG2CurveB point n)) doubledValid =
              toSelectedG2Point
                  (cppBinaryMul selectedG2CurveB point n) accumulatedValid +
                toSelectedG2Point
                  (cppBinaryMul selectedG2CurveB point n) accumulatedValid :=
            doubledSemantics
          _ = n • toSelectedG2Point point pointValid +
                n • toSelectedG2Point point pointValid := by rw [ih]
          _ = (2 * n) • toSelectedG2Point point pointValid := by
            rw [mul_nsmul, two_nsmul, nsmul_add]
          _ = Nat.bit false n • toSelectedG2Point point pointValid := by rfl
      · simp only [cpp_binary_mul_bit selectedG2CurveB point true n digitValid,
          if_true]
        have sumValid := cpp_complete_add_projectively_valid selectedG2CurveB
          selected_g2_no_affine_two_torsion selected_g2_two_ne_zero
          _ _ doubledValid pointValid
        calc
          toSelectedG2Point
              (cppCompleteAdd selectedG2CurveB
                (cppCompleteAdd selectedG2CurveB
                  (cppBinaryMul selectedG2CurveB point n)
                  (cppBinaryMul selectedG2CurveB point n)) point) sumValid =
              toSelectedG2Point
                  (cppCompleteAdd selectedG2CurveB
                    (cppBinaryMul selectedG2CurveB point n)
                    (cppBinaryMul selectedG2CurveB point n)) doubledValid +
                toSelectedG2Point point pointValid :=
            selected_g2_complete_add_eq_add _ _ doubledValid pointValid
          _ = (toSelectedG2Point
                  (cppBinaryMul selectedG2CurveB point n) accumulatedValid +
                toSelectedG2Point
                  (cppBinaryMul selectedG2CurveB point n) accumulatedValid) +
                toSelectedG2Point point pointValid := by rw [doubledSemantics]
          _ = (n • toSelectedG2Point point pointValid +
                n • toSelectedG2Point point pointValid) +
                toSelectedG2Point point pointValid := by rw [ih]
          _ = (2 * n + 1) • toSelectedG2Point point pointValid := by
            rw [add_nsmul, mul_nsmul, two_nsmul, nsmul_add, one_nsmul]
          _ = Nat.bit true n • toSelectedG2Point point pointValid := by rfl

end CatBLST.AVX.MatrixNoK
