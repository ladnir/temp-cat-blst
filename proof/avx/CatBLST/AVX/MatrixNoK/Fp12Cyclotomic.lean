/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Fp12CyclotomicPredicate

set_option maxHeartbeats 800000
set_option maxRecDepth 10000

/-!
Correctness of the policy-selected MatrixNoK cyclotomic square.

The core algebraic statement is an unconditional residual identity between
three literal flat schedules.  The optimized square differs from the ordinary
square by a fixed diagonal scaling of the cyclotomic-predicate residual.
-/

namespace CatBLST.AVX.MatrixNoK

def cyclotomicResidualCorrection
    (residual : FlatFp12 SelectedBaseField) :
    FlatFp12 SelectedBaseField :=
  ![-2 * residual 0, -2 * residual 1,
    (-2 * frob1_2) * residual 2, (-2 * frob1_2) * residual 3,
    (2 * frob1_3) * residual 4, (2 * frob1_3) * residual 5,
    (2 * frob1_3) * residual 6, (2 * frob1_3) * residual 7,
    -2 * residual 8, -2 * residual 9,
    (-2 * frob1_2) * residual 10,
    (-2 * frob1_2) * residual 11]

def cyclotomicPredicateResidual
    (x : FlatFp12 SelectedBaseField) :
    FlatFp12 SelectedBaseField :=
  cppFp12Frobenius2Flat x -
    cppFp12LiteralMul 0 (cppFp12Frobenius4Flat x) x

theorem cyclotomic_residual_coordinate0
    (x : FlatFp12 SelectedBaseField) :
    cppFp12CyclotomicSquareFlat x 0 =
      cppFp12LiteralSquare 0 x 0 +
        cyclotomicResidualCorrection (cyclotomicPredicateResidual x) 0 := by
  simp [cppFp12CyclotomicSquareFlat, cppFp12CyclotomicFirstSplit,
    cppFp12CyclotomicFirstPolynomials, cppFp12LiteralSquare,
    cyclotomicResidualCorrection, cyclotomicPredicateResidual,
    cppFp12Frobenius2Flat, cppFp12Frobenius4Flat,
    cppFp12LiteralMul, cppFp12LiteralMulJoin,
    cppFp12LiteralMulFirstSplit, cppFp12LiteralMulSecondSplit]
  rw [cyclotomic_fraction_facts.2,
    frobenius2_coefficient_linear_facts.1]
  ring_nf
  linear_combination (-2 * x 0) * cyclotomic_fraction_facts.1

theorem cyclotomic_residual_coordinate1
    (x : FlatFp12 SelectedBaseField) :
    cppFp12CyclotomicSquareFlat x 1 =
      cppFp12LiteralSquare 0 x 1 +
        cyclotomicResidualCorrection (cyclotomicPredicateResidual x) 1 := by
  simp [cppFp12CyclotomicSquareFlat, cppFp12CyclotomicFirstSplit,
    cppFp12CyclotomicFirstPolynomials, cppFp12LiteralSquare,
    cyclotomicResidualCorrection, cyclotomicPredicateResidual,
    cppFp12Frobenius2Flat, cppFp12Frobenius4Flat,
    cppFp12LiteralMul, cppFp12LiteralMulJoin,
    cppFp12LiteralMulFirstSplit, cppFp12LiteralMulSecondSplit]
  rw [cyclotomic_fraction_facts.2,
    frobenius2_coefficient_linear_facts.1]
  ring_nf
  linear_combination (-2 * x 1) * cyclotomic_fraction_facts.1

theorem cyclotomic_residual_coordinate2
    (x : FlatFp12 SelectedBaseField) :
    cppFp12CyclotomicSquareFlat x 2 =
      cppFp12LiteralSquare 0 x 2 +
        cyclotomicResidualCorrection (cyclotomicPredicateResidual x) 2 := by
  simp [cppFp12CyclotomicSquareFlat, cppFp12CyclotomicFirstSplit,
    cppFp12CyclotomicFirstPolynomials, cppFp12LiteralSquare,
    cyclotomicResidualCorrection, cyclotomicPredicateResidual,
    cppFp12Frobenius2Flat, cppFp12Frobenius4Flat,
    cppFp12LiteralMul, cppFp12LiteralMulJoin,
    cppFp12LiteralMulFirstSplit, cppFp12LiteralMulSecondSplit]
  rw [cyclotomic_fraction_facts.2,
    frobenius2_coefficient_linear_facts.1]
  ring_nf
  rw [frobenius2_coefficient_square_facts.2.1,
    frobenius2_coefficient_linear_facts.1]
  ring_nf
  linear_combination (-2 * x 2) * cyclotomic_fraction_facts.1

theorem cyclotomic_residual_coordinate3
    (x : FlatFp12 SelectedBaseField) :
    cppFp12CyclotomicSquareFlat x 3 =
      cppFp12LiteralSquare 0 x 3 +
        cyclotomicResidualCorrection (cyclotomicPredicateResidual x) 3 := by
  simp [cppFp12CyclotomicSquareFlat, cppFp12CyclotomicFirstSplit,
    cppFp12CyclotomicFirstPolynomials, cppFp12LiteralSquare,
    cyclotomicResidualCorrection, cyclotomicPredicateResidual,
    cppFp12Frobenius2Flat, cppFp12Frobenius4Flat,
    cppFp12LiteralMul, cppFp12LiteralMulJoin,
    cppFp12LiteralMulFirstSplit, cppFp12LiteralMulSecondSplit]
  rw [cyclotomic_fraction_facts.2,
    frobenius2_coefficient_linear_facts.1]
  ring_nf
  rw [frobenius2_coefficient_square_facts.2.1,
    frobenius2_coefficient_linear_facts.1]
  ring_nf
  linear_combination (-2 * x 3) * cyclotomic_fraction_facts.1

theorem cyclotomic_residual_coordinate4
    (x : FlatFp12 SelectedBaseField) :
    cppFp12CyclotomicSquareFlat x 4 =
      cppFp12LiteralSquare 0 x 4 +
        cyclotomicResidualCorrection (cyclotomicPredicateResidual x) 4 := by
  simp [cppFp12CyclotomicSquareFlat, cppFp12CyclotomicFirstSplit,
    cppFp12CyclotomicFirstPolynomials, cppFp12LiteralSquare,
    cyclotomicResidualCorrection, cyclotomicPredicateResidual,
    cppFp12Frobenius2Flat, cppFp12Frobenius4Flat,
    cppFp12LiteralMul, cppFp12LiteralMulJoin,
    cppFp12LiteralMulFirstSplit, cppFp12LiteralMulSecondSplit]
  simp only [cyclotomic_fraction_facts.2,
    frobenius2_coefficient_linear_facts.1,
    frobenius2_coefficient_linear_facts.2.2]
  ring_nf
  rw [frobenius2_coefficient_square_facts.2.1,
    frobenius2_coefficient_linear_facts.1]
  ring_nf
  linear_combination (-2 * x 4) * cyclotomic_fraction_facts.1

theorem cyclotomic_residual_coordinate5
    (x : FlatFp12 SelectedBaseField) :
    cppFp12CyclotomicSquareFlat x 5 =
      cppFp12LiteralSquare 0 x 5 +
        cyclotomicResidualCorrection (cyclotomicPredicateResidual x) 5 := by
  simp [cppFp12CyclotomicSquareFlat, cppFp12CyclotomicFirstSplit,
    cppFp12CyclotomicFirstPolynomials, cppFp12LiteralSquare,
    cyclotomicResidualCorrection, cyclotomicPredicateResidual,
    cppFp12Frobenius2Flat, cppFp12Frobenius4Flat,
    cppFp12LiteralMul, cppFp12LiteralMulJoin,
    cppFp12LiteralMulFirstSplit, cppFp12LiteralMulSecondSplit]
  simp only [cyclotomic_fraction_facts.2,
    frobenius2_coefficient_linear_facts.1,
    frobenius2_coefficient_linear_facts.2.2]
  ring_nf
  rw [frobenius2_coefficient_square_facts.2.1,
    frobenius2_coefficient_linear_facts.1]
  ring_nf
  linear_combination (-2 * x 5) * cyclotomic_fraction_facts.1

theorem cyclotomic_residual_coordinate6
    (x : FlatFp12 SelectedBaseField) :
    cppFp12CyclotomicSquareFlat x 6 =
      cppFp12LiteralSquare 0 x 6 +
        cyclotomicResidualCorrection (cyclotomicPredicateResidual x) 6 := by
  simp [cppFp12CyclotomicSquareFlat, cppFp12CyclotomicSecondSplit,
    cppFp12CyclotomicSecondPolynomials, cppFp12LiteralSquare,
    cyclotomicResidualCorrection, cyclotomicPredicateResidual,
    cppFp12Frobenius2Flat, cppFp12Frobenius4Flat,
    cppFp12LiteralMul, cppFp12LiteralMulJoin,
    cppFp12LiteralMulFirstSplit, cppFp12LiteralMulSecondSplit]
  simp only [frobenius2_coefficient_linear_facts.1,
    frobenius2_coefficient_linear_facts.2.1,
    frobenius2_coefficient_linear_facts.2.2]
  ring_nf
  rw [frobenius2_coefficient_square_facts.2.1,
    frobenius2_coefficient_linear_facts.1]
  ring_nf
  linear_combination (2 * x 6) * cyclotomic_fraction_facts.1

theorem cyclotomic_residual_coordinate7
    (x : FlatFp12 SelectedBaseField) :
    cppFp12CyclotomicSquareFlat x 7 =
      cppFp12LiteralSquare 0 x 7 +
        cyclotomicResidualCorrection (cyclotomicPredicateResidual x) 7 := by
  simp [cppFp12CyclotomicSquareFlat, cppFp12CyclotomicSecondSplit,
    cppFp12CyclotomicSecondPolynomials, cppFp12LiteralSquare,
    cyclotomicResidualCorrection, cyclotomicPredicateResidual,
    cppFp12Frobenius2Flat, cppFp12Frobenius4Flat,
    cppFp12LiteralMul, cppFp12LiteralMulJoin,
    cppFp12LiteralMulFirstSplit, cppFp12LiteralMulSecondSplit]
  simp only [frobenius2_coefficient_linear_facts.1,
    frobenius2_coefficient_linear_facts.2.1,
    frobenius2_coefficient_linear_facts.2.2]
  ring_nf
  rw [frobenius2_coefficient_square_facts.2.1,
    frobenius2_coefficient_linear_facts.1]
  ring_nf
  linear_combination (2 * x 7) * cyclotomic_fraction_facts.1

theorem cyclotomic_residual_coordinate8
    (x : FlatFp12 SelectedBaseField) :
    cppFp12CyclotomicSquareFlat x 8 =
      cppFp12LiteralSquare 0 x 8 +
        cyclotomicResidualCorrection (cyclotomicPredicateResidual x) 8 := by
  simp [cppFp12CyclotomicSquareFlat, cppFp12CyclotomicSecondSplit,
    cppFp12CyclotomicSecondPolynomials, cppFp12LiteralSquare,
    cyclotomicResidualCorrection, cyclotomicPredicateResidual,
    cppFp12Frobenius2Flat, cppFp12Frobenius4Flat,
    cppFp12LiteralMul, cppFp12LiteralMulJoin,
    cppFp12LiteralMulFirstSplit, cppFp12LiteralMulSecondSplit]
  simp only [frobenius2_coefficient_linear_facts.1]
  ring_nf
  linear_combination (2 * x 8) * cyclotomic_fraction_facts.1

theorem cyclotomic_residual_coordinate9
    (x : FlatFp12 SelectedBaseField) :
    cppFp12CyclotomicSquareFlat x 9 =
      cppFp12LiteralSquare 0 x 9 +
        cyclotomicResidualCorrection (cyclotomicPredicateResidual x) 9 := by
  simp [cppFp12CyclotomicSquareFlat, cppFp12CyclotomicSecondSplit,
    cppFp12CyclotomicSecondPolynomials, cppFp12LiteralSquare,
    cyclotomicResidualCorrection, cyclotomicPredicateResidual,
    cppFp12Frobenius2Flat, cppFp12Frobenius4Flat,
    cppFp12LiteralMul, cppFp12LiteralMulJoin,
    cppFp12LiteralMulFirstSplit, cppFp12LiteralMulSecondSplit]
  simp only [frobenius2_coefficient_linear_facts.1]
  ring_nf
  linear_combination (2 * x 9) * cyclotomic_fraction_facts.1

theorem cyclotomic_residual_coordinate10
    (x : FlatFp12 SelectedBaseField) :
    cppFp12CyclotomicSquareFlat x 10 =
      cppFp12LiteralSquare 0 x 10 +
        cyclotomicResidualCorrection (cyclotomicPredicateResidual x) 10 := by
  simp [cppFp12CyclotomicSquareFlat, cppFp12CyclotomicSecondSplit,
    cppFp12CyclotomicSecondPolynomials, cppFp12LiteralSquare,
    cyclotomicResidualCorrection, cyclotomicPredicateResidual,
    cppFp12Frobenius2Flat, cppFp12Frobenius4Flat,
    cppFp12LiteralMul, cppFp12LiteralMulJoin,
    cppFp12LiteralMulFirstSplit, cppFp12LiteralMulSecondSplit]
  simp only [frobenius2_coefficient_linear_facts.1,
    frobenius2_coefficient_linear_facts.2.2]
  ring_nf
  rw [frobenius2_coefficient_square_facts.2.1,
    frobenius2_coefficient_linear_facts.1]
  ring_nf
  linear_combination (2 * x 10) * cyclotomic_fraction_facts.1

theorem cyclotomic_residual_coordinate11
    (x : FlatFp12 SelectedBaseField) :
    cppFp12CyclotomicSquareFlat x 11 =
      cppFp12LiteralSquare 0 x 11 +
        cyclotomicResidualCorrection (cyclotomicPredicateResidual x) 11 := by
  simp [cppFp12CyclotomicSquareFlat, cppFp12CyclotomicSecondSplit,
    cppFp12CyclotomicSecondPolynomials, cppFp12LiteralSquare,
    cyclotomicResidualCorrection, cyclotomicPredicateResidual,
    cppFp12Frobenius2Flat, cppFp12Frobenius4Flat,
    cppFp12LiteralMul, cppFp12LiteralMulJoin,
    cppFp12LiteralMulFirstSplit, cppFp12LiteralMulSecondSplit]
  simp only [frobenius2_coefficient_linear_facts.1,
    frobenius2_coefficient_linear_facts.2.2]
  ring_nf
  rw [frobenius2_coefficient_square_facts.2.1,
    frobenius2_coefficient_linear_facts.1]
  ring_nf
  linear_combination (2 * x 11) * cyclotomic_fraction_facts.1

theorem cpp_cyclotomic_square_residual_identity
    (x : FlatFp12 SelectedBaseField) :
    cppFp12CyclotomicSquareFlat x =
      cppFp12LiteralSquare 0 x +
        cyclotomicResidualCorrection (cyclotomicPredicateResidual x) := by
  funext index
  fin_cases index
  · exact cyclotomic_residual_coordinate0 x
  · exact cyclotomic_residual_coordinate1 x
  · exact cyclotomic_residual_coordinate2 x
  · exact cyclotomic_residual_coordinate3 x
  · exact cyclotomic_residual_coordinate4 x
  · exact cyclotomic_residual_coordinate5 x
  · exact cyclotomic_residual_coordinate6 x
  · exact cyclotomic_residual_coordinate7 x
  · exact cyclotomic_residual_coordinate8 x
  · exact cyclotomic_residual_coordinate9 x
  · exact cyclotomic_residual_coordinate10 x
  · exact cyclotomic_residual_coordinate11 x

theorem cpp_cyclotomic_square_flat_correct
    (x : FlatFp12 SelectedBaseField)
    (h : cppFp12IsCyclotomicFlat x) :
    cppFp12CyclotomicSquareFlat x = cppFp12LiteralSquare 0 x := by
  rw [cpp_cyclotomic_square_residual_identity]
  have residualZero : cyclotomicPredicateResidual x = 0 := by
    unfold cyclotomicPredicateResidual cppFp12IsCyclotomicFlat at *
    rw [h, sub_self]
  rw [residualZero]
  funext index
  fin_cases index <;> simp [cyclotomicResidualCorrection]

theorem cpp_cyclotomic_square_correct_of_direct_predicate
    (x : SelectedFp12) (h : cppFp12IsCyclotomic x) :
    cppFp12CyclotomicSquare x = x * x := by
  have flatPredicate : cppFp12IsCyclotomicFlat (unpackFlatFp12 x) :=
    (cpp_fp12_is_cyclotomic_iff_flat x).mp h
  have flatCorrect := cpp_cyclotomic_square_flat_correct
    (unpackFlatFp12 x) flatPredicate
  calc
    cppFp12CyclotomicSquare x =
        packFlatFp12 (cppFp12LiteralSquare 0 (unpackFlatFp12 x)) :=
      congrArg packFlatFp12 flatCorrect
    _ = cppFp12FlatSquare (packFlatFp12 (unpackFlatFp12 x)) :=
      cpp_fp12_literal_square_packs _
    _ = cppFp12FlatSquare x := by rw [pack_unpack_fp12]
    _ = Fp12Value.square x := cpp_fp12_flat_square_correct x
    _ = x * x := rfl

theorem cpp_cyclotomic_square_correct
    (x : SelectedFp12) (h : cppFp12IsCyclotomicSource x) :
    cppFp12CyclotomicSquare x = x * x :=
  cpp_cyclotomic_square_correct_of_direct_predicate x
    ((cpp_fp12_is_cyclotomic_source_iff_direct x).mp h)

end CatBLST.AVX.MatrixNoK
