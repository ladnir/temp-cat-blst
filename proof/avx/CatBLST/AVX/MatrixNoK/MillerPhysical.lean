/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.MillerSchedule
import Mathlib.Tactic

/-!
# MatrixNoK Miller physical bounds and call ledger

Miller `mul_accumulate` coordinates contain two Fp2 products, hence four base
products rather than the two-product envelope used by the original P2 packet.
This file establishes the larger semantic and physical type bounds and records
the exact direct and prepared source call counts.  Reducer composition against
the generalized envelope is the remaining P4 physical obligation.
-/

namespace CatBLST.AVX.MatrixNoK

/-- Four base products of operands bounded by `2m` fit below `16m^2`. -/
def millerFp2AccumulatorBound : Nat :=
  4 * fp2PreparedBound * fp2PreparedBound

theorem selected_miller_fp2_accumulator_bound :
    millerFp2AccumulatorBound = 16 := by
  norm_num [millerFp2AccumulatorBound, fp2PreparedBound, standardBound]

theorem miller_fp2_four_product_accumulator_bound
    (a0 a1 a2 a3 b0 b1 b2 b3 modulus : Nat)
    (a0Fit : a0 ≤ fp2PreparedBound * modulus)
    (a1Fit : a1 ≤ fp2PreparedBound * modulus)
    (a2Fit : a2 ≤ fp2PreparedBound * modulus)
    (a3Fit : a3 ≤ fp2PreparedBound * modulus)
    (b0Fit : b0 ≤ fp2PreparedBound * modulus)
    (b1Fit : b1 ≤ fp2PreparedBound * modulus)
    (b2Fit : b2 ≤ fp2PreparedBound * modulus)
    (b3Fit : b3 ≤ fp2PreparedBound * modulus) :
    a0 * b0 + a1 * b1 + a2 * b2 + a3 * b3 ≤
      millerFp2AccumulatorBound * modulus * modulus := by
  calc
    a0 * b0 + a1 * b1 + a2 * b2 + a3 * b3 ≤
        (fp2PreparedBound * modulus) * (fp2PreparedBound * modulus) +
        (fp2PreparedBound * modulus) * (fp2PreparedBound * modulus) +
        (fp2PreparedBound * modulus) * (fp2PreparedBound * modulus) +
        (fp2PreparedBound * modulus) * (fp2PreparedBound * modulus) := by
      gcongr
    _ = millerFp2AccumulatorBound * modulus * modulus := by
      simp [millerFp2AccumulatorBound]
      ring

theorem miller_fp2_accumulator_fits_delayed_input :
    millerFp2AccumulatorBound ≤ delayedInputBound := by
  norm_num [millerFp2AccumulatorBound, fp2PreparedBound,
    delayedInputBound, standardBound]

theorem miller_fp2_completed_m1_fits_selected_type :
    millerFp2AccumulatorBound + 1 ≤ delayedCombinedTypeBound := by
  norm_num [millerFp2AccumulatorBound, fp2PreparedBound,
    delayedCombinedTypeBound, delayedInputBound, delayedR1TypeBound,
    r1ExactAccumulationBound, multiplicationProductBound, standardBound]

theorem miller_fp2_completed_m2_fits_selected_type :
    millerFp2AccumulatorBound + 1 + 16 ≤ delayedCombinedTypeBound := by
  norm_num [millerFp2AccumulatorBound, fp2PreparedBound,
    delayedCombinedTypeBound, delayedInputBound, delayedR1TypeBound,
    r1ExactAccumulationBound, multiplicationProductBound, standardBound]

theorem miller_fp2_completion_offset_dominates_signed_expression :
    (4 * (fp2BaseRedundantMultiple * Generated.target_modulus) ^ 2 : Int) <
      fp2CompletionOffset := by
  norm_num [fp2CompletionOffset, fp2BaseRedundantMultiple,
    fp2MaxNegativeRns]
  have targetPositive : (0 : Int) < Generated.target_modulus := by
    norm_num [Generated.target_modulus]
  nlinarith [sq_pos_of_pos targetPositive]

theorem miller_fp2_completed_expression_under_five_target_products :
    (fp2CompletionOffset : Int) +
        4 * (fp2BaseRedundantMultiple * Generated.target_modulus) ^ 2 <
      5 * Generated.target_modulus * Generated.M1 := by
  have coefficientBound :
      (3097600 : Int) * Generated.target_modulus < 5 * Generated.M1 := by
    norm_num [Generated.target_modulus, Generated.M1]
  have targetPositive : (0 : Int) < Generated.target_modulus := by
    norm_num [Generated.target_modulus]
  have scaled := mul_lt_mul_of_pos_right coefficientBound targetPositive
  norm_num [fp2CompletionOffset, fp2BaseRedundantMultiple,
    fp2MaxNegativeRns]
  nlinarith

theorem miller_fp2_completed_r1_representative_range
    (value : Int) (residues : Fin 8 → Nat)
    (valueBound : |value| <
      4 * (fp2BaseRedundantMultiple * Generated.target_modulus) ^ 2)
    (residuesFit : ∀ index,
      residues index < 2 * listAt8 Generated.moduli1 index.val)
    (residuesRepresent : ∀ index,
      Nat.ModEq (modulus1At index) (residues index)
        (signedR1ReducedValue (cppFp2CompletedExpression value))) :
    0 ≤ r1SignedValueFromResidues
        (cppFp2CompletedExpression value) residues ∧
      r1SignedValueFromResidues
          (cppFp2CompletedExpression value) residues <
        22 * Generated.target_modulus := by
  have valueBounds := abs_lt.mp valueBound
  have completedNonnegative :
      (0 : Int) ≤ cppFp2CompletedExpression value := by
    unfold cppFp2CompletedExpression
    have dominates := miller_fp2_completion_offset_dominates_signed_expression
    nlinarith
  have completedUpper :
      cppFp2CompletedExpression value <
        5 * Generated.target_modulus * Generated.M1 := by
    unfold cppFp2CompletedExpression
    have upper := miller_fp2_completed_expression_under_five_target_products
    nlinarith
  have liftBoundNat := r1_standard_raw_lift_lt_17_product residues
    (by simpa [modulus1At] using residuesFit)
  have liftBound :
      (r1RawLiftFromResidues residues : Int) < 17 * Generated.M1 := by
    exact_mod_cast liftBoundNat
  have liftNonnegative :
      (0 : Int) ≤ r1RawLiftFromResidues residues := by
    exact_mod_cast Nat.zero_le (r1RawLiftFromResidues residues)
  have reconstruction := r1_signed_value_reconstruction
    (cppFp2CompletedExpression value) residues residuesRepresent
  have modulusPositive : (0 : Int) < Generated.M1 := by
    norm_num [Generated.M1]
  unfold r1SignedValueNumerator at reconstruction
  constructor <;> nlinarith

/-!
Each entry below is the number of base-product terms in one component of an
Fp2 value passed to `batch_reduce_expand`.  Squares and ordinary Fp2 products
have two terms per component; the accumulated point outputs have four; line
scales and normalization-by-one have one.  Integer scalar multiplies and
add/sub preparation introduce no additional base-field products.
-/

def directMillerDoubleRawTerms : List Nat :=
  [2, 2, 2, 2, 2, 1, 1, 2, 4, 2]

def directMillerAddRawTerms : List Nat :=
  [2, 2, 2, 2, 2, 2, 4, 4, 4, 1, 1]

def preparedMillerDoubleRawTerms : List Nat :=
  [2, 2, 2, 2, 2, 4, 4, 4, 1, 1, 1]

def preparedMillerAddRawTerms : List Nat :=
  [2, 2, 2, 2, 2, 2, 4, 4, 4, 1, 1, 1]

theorem every_miller_fp2_reduction_fits_four_product_envelope :
    ∀ terms ∈ directMillerDoubleRawTerms ++ directMillerAddRawTerms ++
        preparedMillerDoubleRawTerms ++ preparedMillerAddRawTerms,
      terms ≤ 4 := by
  native_decide

theorem miller_fp2_batch_permutation_reconstructs
    {R : Type*} (values : List (Fp2Value R))
    (index : Fin values.length) :
    ⟨cppFp2StridedXAt values index, cppFp2StridedYAt values index⟩ =
      values.get index := by
  exact cpp_fp2_strided_reconstruct values index

structure MillerStepPhysicalLedger where
  curveBaseProducts : Nat
  sparseBaseProducts : Nat
  reducedFp2Coordinates : Nat
  normalizedLineCoordinates : Nat
  deriving DecidableEq

/-- Direct `Miller::double_step`, including application of its line. -/
def directMillerDoubleLedger : MillerStepPhysicalLedger :=
  { curveBaseProducts := 37
    sparseBaseProducts := 72
    reducedFp2Coordinates := 10
    normalizedLineCoordinates := 0 }

/-- Direct `Miller::add_step`, including application of its line. -/
def directMillerAddLedger : MillerStepPhysicalLedger :=
  { curveBaseProducts := 52
    sparseBaseProducts := 72
    reducedFp2Coordinates := 11
    normalizedLineCoordinates := 0 }

/-- Prepared doubling plus its later G1 scaling and line application. -/
def preparedMillerDoubleLedger : MillerStepPhysicalLedger :=
  { curveBaseProducts := 43
    sparseBaseProducts := 72
    reducedFp2Coordinates := 11
    normalizedLineCoordinates := 3 }

/-- Prepared addition plus its later G1 scaling and line application. -/
def preparedMillerAddLedger : MillerStepPhysicalLedger :=
  { curveBaseProducts := 58
    sparseBaseProducts := 72
    reducedFp2Coordinates := 12
    normalizedLineCoordinates := 3 }

structure MillerLoopPhysicalLedger where
  doubles : Nat
  additions : Nat
  fp12Squares : Nat
  sparseLineProducts : Nat
  curveBaseProducts : Nat
  sparseBaseProducts : Nat
  reducedFp2Coordinates : Nat
  finalConjugations : Nat
  deriving DecidableEq

def directMillerLoopLedger : MillerLoopPhysicalLedger :=
  { doubles := 63
    additions := 5
    fp12Squares := 62
    sparseLineProducts := 68
    curveBaseProducts := 63 * directMillerDoubleLedger.curveBaseProducts +
      5 * directMillerAddLedger.curveBaseProducts
    sparseBaseProducts := 68 * cppFlatMulXY00Z0Ledger.baseProducts
    reducedFp2Coordinates :=
      63 * directMillerDoubleLedger.reducedFp2Coordinates +
        5 * directMillerAddLedger.reducedFp2Coordinates
    finalConjugations := 1 }

def preparedMillerLoopLedger : MillerLoopPhysicalLedger :=
  { doubles := 63
    additions := 5
    fp12Squares := 62
    sparseLineProducts := 68
    curveBaseProducts := 63 * preparedMillerDoubleLedger.curveBaseProducts +
      5 * preparedMillerAddLedger.curveBaseProducts
    sparseBaseProducts := 68 * cppFlatMulXY00Z0Ledger.baseProducts
    reducedFp2Coordinates :=
      63 * preparedMillerDoubleLedger.reducedFp2Coordinates +
        5 * preparedMillerAddLedger.reducedFp2Coordinates
    finalConjugations := 1 }

theorem direct_miller_loop_ledger_exact :
    directMillerLoopLedger =
      { doubles := 63, additions := 5, fp12Squares := 62,
        sparseLineProducts := 68, curveBaseProducts := 2591,
        sparseBaseProducts := 4896, reducedFp2Coordinates := 685,
        finalConjugations := 1 } := by
  native_decide

theorem prepared_miller_loop_ledger_exact :
    preparedMillerLoopLedger =
      { doubles := 63, additions := 5, fp12Squares := 62,
        sparseLineProducts := 68, curveBaseProducts := 2999,
        sparseBaseProducts := 4896, reducedFp2Coordinates := 753,
        finalConjugations := 1 } := by
  native_decide

end CatBLST.AVX.MatrixNoK
