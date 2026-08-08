/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.SignedPipeline
import Mathlib.Tactic

/-!
# Selected MatrixNoK Fp2 construction

This file models the selected `FP2Ring<FpRingCandidate>` implementation in
`blst/avx/fields/fp2.h`.  It intentionally does not model the comparison-only
`matrixnok::Fp2RingCandidate` microkernel.

The C++ representation is `x + u*y` with `u^2 = -1`.  `prep_left` retains
`x`, `y`, and a prepared `-y`; multiplication then forms `ac - bd` and
`ad + bc` without adding a branch or a runtime abstraction.  Component
reduction is semantically the identity in the target field, even though the
selected implementation uses the delayed `800 + 32 = 832` RNS path.
-/

namespace CatBLST.AVX.MatrixNoK

/-- A semantic Fp2 value, matching the C++ `(x, y)` storage order. -/
structure Fp2Value (R : Type*) where
  x : R
  y : R
deriving DecidableEq

@[ext]
theorem Fp2Value.ext {R : Type*} {a b : Fp2Value R}
    (xEqual : a.x = b.x) (yEqual : a.y = b.y) : a = b := by
  cases a
  cases b
  simp_all

namespace Fp2Value

variable {R : Type*}

def map (f : R → R) (a : Fp2Value R) : Fp2Value R :=
  ⟨f a.x, f a.y⟩

def add [Add R] (a b : Fp2Value R) : Fp2Value R :=
  ⟨a.x + b.x, a.y + b.y⟩

def sub [Sub R] (a b : Fp2Value R) : Fp2Value R :=
  ⟨a.x - b.x, a.y - b.y⟩

def neg [Neg R] (a : Fp2Value R) : Fp2Value R :=
  ⟨-a.x, -a.y⟩

/-- Multiplication in `R[u]/(u^2 + 1)`. -/
def mul [Ring R] (a b : Fp2Value R) : Fp2Value R :=
  ⟨a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x⟩

def square [Ring R] (a : Fp2Value R) : Fp2Value R :=
  ⟨a.x * a.x - a.y * a.y, 2 * (a.x * a.y)⟩

/-- Multiplication by the BLS12-381 tower nonresidue `1 + u`. -/
def mulByNonresidue [Ring R] (a : Fp2Value R) : Fp2Value R :=
  ⟨a.x - a.y, a.x + a.y⟩

/-- Multiplication by `3b = 12(1 + u)` for the G2 curve coefficient. -/
def mul3b [Ring R] (a : Fp2Value R) : Fp2Value R :=
  map (fun value => 12 * value) (mulByNonresidue a)

end Fp2Value

/-- The three values materialized by C++ `FP2::prep_left`. -/
structure PreparedFp2 (R : Type*) where
  x : R
  y : R
  minusY : R

/-- Source-shaped model of `FP2::prep_left`. -/
def cppFp2PrepLeft {R : Type*}
    (prep negate : R → R) (a : Fp2Value R) : PreparedFp2 R :=
  ⟨prep a.x, prep a.y, prep (negate (prep a.y))⟩

/-- Source-shaped model of `FP2N::operator*`. -/
def cppFp2DelayedMul {R : Type*} [Ring R]
    (a : PreparedFp2 R) (b : Fp2Value R) : Fp2Value R :=
  ⟨a.x * b.x + a.minusY * b.y,
    a.x * b.y + a.y * b.x⟩

/-- Source-shaped model of `FP2N::square`. -/
def cppFp2DelayedSquare {R : Type*} [Ring R]
    (a : PreparedFp2 R) : Fp2Value R :=
  ⟨a.x * a.x + a.minusY * a.y, (a.x * a.y) * 2⟩

def PreparedFp2.toValue {R : Type*} (a : PreparedFp2 R) : Fp2Value R :=
  ⟨a.x, a.y⟩

/-- Source-shaped model of `FP2::mul_accumulate`. -/
def cppFp2MulAccumulate {R : Type*} [Ring R]
    (accumulator : Fp2Value R) (a : PreparedFp2 R)
    (b : Fp2Value R) : Fp2Value R :=
  ⟨accumulator.x + a.x * b.x + a.minusY * b.y,
    accumulator.y + a.x * b.y + a.y * b.x⟩

/-- Source-shaped model of `FP2::mul_sub_accumulate`. -/
def cppFp2MulSubAccumulate {R : Type*} [Ring R]
    (negate : R → R) (accumulator : Fp2Value R) (a : PreparedFp2 R)
    (b : Fp2Value R) : Fp2Value R :=
  ⟨accumulator.x + negate a.x * b.x + a.y * b.y,
    accumulator.y + negate a.x * b.y + a.minusY * b.x⟩

def cppFp2MulScalar {R : Type*} [Mul R]
    (a : Fp2Value R) (scalar : R) : Fp2Value R :=
  ⟨a.x * scalar, a.y * scalar⟩

/-- The variadic C++ ready array layout: all `x` values, then all `y` values. -/
def cppFp2StridedLayout {R : Type*} (values : List (Fp2Value R)) : List R :=
  values.map Fp2Value.x ++ values.map Fp2Value.y

def cppFp2StridedXAt {R : Type*} (values : List (Fp2Value R))
    (index : Fin values.length) : R :=
  (cppFp2StridedLayout values).get ⟨index.val, by
    simp [cppFp2StridedLayout]
    omega⟩

def cppFp2StridedYAt {R : Type*} (values : List (Fp2Value R))
    (index : Fin values.length) : R :=
  (cppFp2StridedLayout values).get ⟨index.val + values.length, by
    simp [cppFp2StridedLayout]⟩

theorem cpp_fp2_strided_x {R : Type*} (values : List (Fp2Value R))
    (index : Fin values.length) :
    cppFp2StridedXAt values index = (values.get index).x := by
  simp [cppFp2StridedXAt, cppFp2StridedLayout]

theorem cpp_fp2_strided_y {R : Type*} (values : List (Fp2Value R))
    (index : Fin values.length) :
    cppFp2StridedYAt values index = (values.get index).y := by
  simp [cppFp2StridedYAt, cppFp2StridedLayout]

theorem cpp_fp2_strided_reconstruct {R : Type*}
    (values : List (Fp2Value R)) (index : Fin values.length) :
    ⟨cppFp2StridedXAt values index, cppFp2StridedYAt values index⟩ =
      values.get index := by
  rw [cpp_fp2_strided_x values index, cpp_fp2_strided_y values index]

/--
The selected multiplication construction, parameterized by the base-field
preparation, negation, and delayed-finish operations used by C++.
-/
def cppSelectedFp2Mul {R : Type*} [Ring R]
    (prep negate finish : R → R) (a b : Fp2Value R) : Fp2Value R :=
  Fp2Value.map finish
    (cppFp2DelayedMul (cppFp2PrepLeft prep negate a)
      (Fp2Value.map prep b))

def cppSelectedFp2Square {R : Type*} [Ring R]
    (prep negate finish : R → R) (a : Fp2Value R) : Fp2Value R :=
  Fp2Value.map finish
    (cppFp2DelayedSquare (cppFp2PrepLeft prep negate a))

theorem cpp_selected_fp2_mul_correct {R : Type*} [CommRing R]
    (prep negate finish : R → R)
    (prepCorrect : ∀ value, prep value = value)
    (negateCorrect : ∀ value, negate value = -value)
    (finishCorrect : ∀ value, finish value = value)
    (a b : Fp2Value R) :
    cppSelectedFp2Mul prep negate finish a b = Fp2Value.mul a b := by
  ext
  · simp [cppSelectedFp2Mul, cppFp2DelayedMul, cppFp2PrepLeft,
      Fp2Value.map, Fp2Value.mul, prepCorrect, negateCorrect, finishCorrect]
    ring
  ·
    simp [cppSelectedFp2Mul, cppFp2DelayedMul, cppFp2PrepLeft,
      Fp2Value.map, Fp2Value.mul, prepCorrect, negateCorrect, finishCorrect]

theorem cpp_selected_fp2_square_correct {R : Type*} [CommRing R]
    (prep negate finish : R → R)
    (prepCorrect : ∀ value, prep value = value)
    (negateCorrect : ∀ value, negate value = -value)
    (finishCorrect : ∀ value, finish value = value)
    (a : Fp2Value R) :
    cppSelectedFp2Square prep negate finish a = Fp2Value.square a := by
  ext
  · simp [cppSelectedFp2Square, cppFp2DelayedSquare, cppFp2PrepLeft,
      Fp2Value.map, Fp2Value.square, prepCorrect, negateCorrect,
      finishCorrect]
    ring
  · simp [cppSelectedFp2Square, cppFp2DelayedSquare, cppFp2PrepLeft,
      Fp2Value.map, Fp2Value.square, prepCorrect, negateCorrect,
      finishCorrect]
    ring

theorem cpp_fp2_mul_accumulate_correct {R : Type*} [CommRing R]
    (accumulator b : Fp2Value R) (a : PreparedFp2 R)
    (minusYCorrect : a.minusY = -a.y) :
    cppFp2MulAccumulate accumulator a b =
      Fp2Value.add accumulator (Fp2Value.mul a.toValue b) := by
  ext <;>
    simp [cppFp2MulAccumulate, PreparedFp2.toValue, Fp2Value.add,
      Fp2Value.mul, minusYCorrect] <;>
    ring

theorem cpp_fp2_mul_sub_accumulate_correct {R : Type*} [CommRing R]
    (negate : R → R) (accumulator b : Fp2Value R) (a : PreparedFp2 R)
    (negateCorrect : ∀ value, negate value = -value)
    (minusYCorrect : a.minusY = -a.y) :
    cppFp2MulSubAccumulate negate accumulator a b =
      Fp2Value.sub accumulator (Fp2Value.mul a.toValue b) := by
  ext <;>
    simp [cppFp2MulSubAccumulate, PreparedFp2.toValue, Fp2Value.sub,
      Fp2Value.mul, negateCorrect, minusYCorrect] <;>
    ring

theorem cpp_fp2_mul_scalar_formula {R : Type*} [Mul R]
    (a : Fp2Value R) (scalar : R) :
    cppFp2MulScalar a scalar = Fp2Value.map (fun value => value * scalar) a :=
  rfl

theorem cpp_fp2_mul_by_nonresidue_formula {R : Type*} [CommRing R]
    (a : Fp2Value R) :
    ⟨a.x + -a.y, a.x + a.y⟩ = Fp2Value.mulByNonresidue a := by
  ext <;> simp [Fp2Value.mulByNonresidue, sub_eq_add_neg]

theorem cpp_fp2_mul_3b_formula {R : Type*} [CommRing R]
    (a : Fp2Value R) :
    Fp2Value.map (fun value => value * 12)
        ⟨a.x + -a.y, a.x + a.y⟩ = Fp2Value.mul3b a := by
  ext <;> simp [Fp2Value.mul3b, Fp2Value.mulByNonresidue,
    Fp2Value.map, sub_eq_add_neg] <;> ring

/--
Standard Fp2 inputs remain below `2m`: `prep<4>` returns its existing
`Bounds<0,2>` input unchanged, and `elementwise_moduli<2>() - y` has the same
bound before the second `prep`.
-/
def fp2PreparedBound : Nat := standardBound

/-- Physical lane value created by `elementwise_moduli<2>() - value`. -/
def cppFp2PreparedNegate (value modulus : Nat) : Nat :=
  fp2PreparedBound * modulus - value

theorem cpp_fp2_prepared_negate_bound
    (value modulus : Nat) :
    cppFp2PreparedNegate value modulus ≤ fp2PreparedBound * modulus := by
  exact Nat.sub_le _ _

theorem cpp_fp2_prepared_negate_semantics
    (value modulus : Nat) (_valueFit : value ≤ fp2PreparedBound * modulus) :
    (cppFp2PreparedNegate value modulus : ZMod modulus) = -value := by
  unfold cppFp2PreparedNegate
  rw [Nat.cast_sub _valueFit]
  simp [fp2PreparedBound, standardBound]

/-- Each two-product Fp2 accumulator fits below `8m^2`. -/
def fp2AccumulatorBound : Nat :=
  2 * fp2PreparedBound * fp2PreparedBound

/-- The tight M2 accumulator plus the exact eight-term R1 contribution. -/
def fp2ExactCombinedBound : Nat :=
  fp2AccumulatorBound + r1ExactAccumulationBound

theorem selected_fp2_bound_values :
    fp2PreparedBound = 2 ∧ fp2AccumulatorBound = 8 ∧
      fp2ExactCombinedBound = 24 := by
  norm_num [fp2PreparedBound, fp2AccumulatorBound, fp2ExactCombinedBound,
    r1ExactAccumulationBound, standardBound]

theorem fp2_two_product_accumulator_bound
    (left0 left1 right0 right1 modulus : Nat)
    (left0Fit : left0 ≤ fp2PreparedBound * modulus)
    (left1Fit : left1 ≤ fp2PreparedBound * modulus)
    (right0Fit : right0 ≤ fp2PreparedBound * modulus)
    (right1Fit : right1 ≤ fp2PreparedBound * modulus) :
    left0 * right0 + left1 * right1 ≤
      fp2AccumulatorBound * modulus * modulus := by
  calc
    left0 * right0 + left1 * right1 ≤
        (fp2PreparedBound * modulus) * (fp2PreparedBound * modulus) +
          (fp2PreparedBound * modulus) * (fp2PreparedBound * modulus) := by
      exact Nat.add_le_add
        (Nat.mul_le_mul left0Fit right0Fit)
        (Nat.mul_le_mul left1Fit right1Fit)
    _ = fp2AccumulatorBound * modulus * modulus := by
      simp [fp2AccumulatorBound]
      ring

theorem fp2_accumulator_fits_delayed_input
    (value modulus : Nat)
    (valueFit : value ≤ fp2AccumulatorBound * modulus * modulus) :
    value ≤ delayedInputBound * modulus * modulus := by
  have boundLe : fp2AccumulatorBound ≤ delayedInputBound := by
    norm_num [fp2AccumulatorBound, fp2PreparedBound,
      multiplicationProductBound, standardBound, delayedInputBound]
  have scaled := Nat.mul_le_mul_right (modulus * modulus) boundLe
  have scaled' : fp2AccumulatorBound * modulus * modulus ≤
      delayedInputBound * modulus * modulus := by
    simpa [Nat.mul_assoc] using scaled
  exact valueFit.trans scaled'

theorem fp2_delayed_combined_bound
    (accumulator r1Contribution modulus : Nat)
    (accumulatorFit : accumulator ≤
      fp2AccumulatorBound * modulus * modulus)
    (r1Fit : r1Contribution ≤ delayedR1TypeBound * modulus * modulus) :
    accumulator + r1Contribution ≤
      delayedCombinedTypeBound * modulus * modulus := by
  have accumulatorWide :=
    fp2_accumulator_fits_delayed_input accumulator modulus accumulatorFit
  unfold delayedCombinedTypeBound
  nlinarith

theorem fp2_exact_combined_bound
    (accumulator r1Contribution modulus : Nat)
    (accumulatorFit : accumulator ≤
      fp2AccumulatorBound * modulus * modulus)
    (r1Fit : r1Contribution ≤
      r1ExactAccumulationBound * modulus * modulus) :
    accumulator + r1Contribution ≤
      fp2ExactCombinedBound * modulus * modulus := by
  unfold fp2ExactCombinedBound
  nlinarith

theorem fp2_exact_combined_fits_selected_type
    (value modulus : Nat)
    (valueFit : value ≤ fp2ExactCombinedBound * modulus * modulus) :
    value ≤ delayedCombinedTypeBound * modulus * modulus := by
  have boundLe : fp2ExactCombinedBound ≤ delayedCombinedTypeBound := by
    norm_num [fp2ExactCombinedBound, fp2AccumulatorBound, fp2PreparedBound,
      r1ExactAccumulationBound, delayedCombinedTypeBound, delayedInputBound,
      delayedR1TypeBound, standardBound]
  have scaled := Nat.mul_le_mul_right (modulus * modulus) boundLe
  have scaled' : fp2ExactCombinedBound * modulus * modulus ≤
      delayedCombinedTypeBound * modulus * modulus := by
    simpa [Nat.mul_assoc] using scaled
  exact valueFit.trans scaled'

/--
Selected delayed reduction is semantically transparent after cancelling the
one Montgomery radix factor.  This is the base-field contract composed by the
Fp2 formula theorems above.
-/
theorem selected_delayed_finish_semantics
    (high low value modulus factor : Nat)
    (modulusMember : modulus ∈ allModuli)
    (inverseRelation : modulus * factor % ifmaRadix = 1)
    (wideFits : high * ifmaRadix + low ≤
      delayedCombinedTypeBound * modulus * modulus)
    (wideRepresents :
      (high * ifmaRadix + low : ZMod modulus) = value * ifmaRadix) :
    (selectedReduceMasked high low modulus factor : ZMod modulus) = value := by
  have contract := selected_delayed_reducer_contract high low modulus factor
    modulusMember inverseRelation wideFits
  have equalTimesRadix :
      (selectedReduceMasked high low modulus factor : ZMod modulus) *
          ifmaRadix = (value : ZMod modulus) * ifmaRadix := by
    rw [contract.1, wideRepresents]
  have radixUnit : IsUnit (ifmaRadix : ZMod modulus) := by
    rw [ZMod.isUnit_iff_coprime]
    unfold ifmaRadix
    rw [Nat.coprime_pow_left_iff (by norm_num)]
    exact (Nat.coprime_two_left).2 (all_moduli_odd modulus modulusMember)
  exact radixUnit.mul_right_cancel equalTimesRadix

/-!
## Physical negative-range completion

`FP2::prep_left` gives the real coordinate a negative static RNS lower bound.
Before `batch_reduce` consumes that wide value, `FpRingCandidate::complete`
therefore adds the generated wide `max_neg_offset` in both bases.  The
semantic value below is the exact generator input for that constant; proving
the generated lane encoding and reducer producer correspondence remains the
last source-shaped step.
-/

def fp2BaseRedundantMultiple : Nat := 40

def fp2MaxNegativeRns : Nat := 1932

def fp2CompletionOffset : Nat :=
  fp2MaxNegativeRns * fp2BaseRedundantMultiple *
    fp2BaseRedundantMultiple * Generated.target_modulus *
      Generated.target_modulus

def cppFp2CompletedExpression (value : Int) : Int :=
  value + fp2CompletionOffset

theorem fp2_completion_constant_values :
    fp2BaseRedundantMultiple = 40 ∧ fp2MaxNegativeRns = 1932 := by
  norm_num [fp2BaseRedundantMultiple, fp2MaxNegativeRns]

/--
The generated M1 completion lanes are the exact wide-Montgomery encoding of
the semantic offset. The second radix/premultiplier factor is intentional:
`complete()` adds the constant before the wide M1 reduction removes one
Montgomery radix.
-/
theorem fp2_max_neg_offset_m1_lane_encoding (index : Fin 8) :
    (listAt8 Generated.max_neg_offset_m1 index.val :
        ZMod (modulus1At index)) =
      ((fp2CompletionOffset * ifmaRadix * ifmaRadix *
        Generated.qr_premult * Generated.qr_premult : Nat) :
          ZMod (modulus1At index)) := by
  apply (ZMod.natCast_eq_natCast_iff _ _ _).2
  fin_cases index <;>
    norm_num [listAt8, modulus1At, Generated.max_neg_offset_m1,
      Generated.moduli1, fp2CompletionOffset, fp2BaseRedundantMultiple,
      fp2MaxNegativeRns, Generated.target_modulus, ifmaRadix,
      Generated.qr_premult]

/-- The corresponding generated M2 lanes use the cross-basis wide encoding. -/
theorem fp2_max_neg_offset_m2_lane_encoding (index : Fin 8) :
    (listAt8 Generated.max_neg_offset_m2 index.val :
        ZMod (modulus2At index)) =
      ((fp2CompletionOffset * ifmaRadix * ifmaRadix *
        crossBasisInverse * crossBasisInverse : Nat) :
          ZMod (modulus2At index)) := by
  apply (ZMod.natCast_eq_natCast_iff _ _ _).2
  fin_cases index <;>
    norm_num [listAt8, modulus2At, Generated.max_neg_offset_m2,
      Generated.moduli2, fp2CompletionOffset, fp2BaseRedundantMultiple,
      fp2MaxNegativeRns, Generated.target_modulus, ifmaRadix,
      crossBasisInverse, montgomeryRadix, Generated.convert_to_postmult]

theorem fp2_completed_m1_wide_lane_encoding
    (value : Int) (rawWide : Nat) (index : Fin 8)
    (rawRepresents :
      (rawWide : ZMod (modulus1At index)) =
        value * ifmaRadix * ifmaRadix * Generated.qr_premult *
          Generated.qr_premult) :
    (rawWide + listAt8 Generated.max_neg_offset_m1 index.val :
        ZMod (modulus1At index)) =
      cppFp2CompletedExpression value * ifmaRadix * ifmaRadix *
        Generated.qr_premult * Generated.qr_premult := by
  rw [rawRepresents, fp2_max_neg_offset_m1_lane_encoding index]
  unfold cppFp2CompletedExpression
  push_cast
  ring

theorem fp2_completed_m2_wide_lane_encoding
    (value : Int) (rawWide : Nat) (index : Fin 8)
    (rawRepresents :
      (rawWide : ZMod (modulus2At index)) =
        value * ifmaRadix * ifmaRadix * crossBasisInverse *
          crossBasisInverse) :
    (rawWide + listAt8 Generated.max_neg_offset_m2 index.val :
        ZMod (modulus2At index)) =
      cppFp2CompletedExpression value * ifmaRadix * ifmaRadix *
        crossBasisInverse * crossBasisInverse := by
  rw [rawRepresents, fp2_max_neg_offset_m2_lane_encoding index]
  unfold cppFp2CompletedExpression
  push_cast
  ring

/--
Source-shaped M1 composition for `complete()` followed by the selected wide
reducer. The physical input carries two radix factors; the reducer removes one
and produces exactly the signed-R1 residue expected by the correction-free
reconstruction theorem.
-/
theorem fp2_completed_m1_reducer_represents
    (value : Int) (rawWide high low : Nat) (index : Fin 8)
    (rawRepresents :
      (rawWide : ZMod (modulus1At index)) =
        value * ifmaRadix * ifmaRadix * Generated.qr_premult *
          Generated.qr_premult)
    (wideSplit : high * ifmaRadix + low =
      rawWide + listAt8 Generated.max_neg_offset_m1 index.val)
    (wideFits : high * ifmaRadix + low ≤
      delayedCombinedTypeBound * modulus1At index * modulus1At index) :
    Nat.ModEq (modulus1At index)
      (selectedReduceMasked high low (modulus1At index)
        (reducer1FactorAt index))
      (signedR1ReducedValue (cppFp2CompletedExpression value)) := by
  let modulus := modulus1At index
  have reducer := selected_delayed_reducer_contract high low modulus
    (reducer1FactorAt index) (modulus1_member index)
    (reducer1_factor_relation index) (by simpa [modulus] using wideFits)
  have splitCast := congrArg
    (fun result : Nat => (result : ZMod modulus)) wideSplit
  push_cast at splitCast
  have completedCast := fp2_completed_m1_wide_lane_encoding
    value rawWide index rawRepresents
  have targetCast := signed_r1_reduced_value_lane
    (cppFp2CompletedExpression value) index
  have completedCast' :
      (rawWide + listAt8 Generated.max_neg_offset_m1 index.val :
          ZMod modulus) =
        cppFp2CompletedExpression value * ifmaRadix * ifmaRadix *
          Generated.qr_premult * Generated.qr_premult := by
    simpa [modulus] using completedCast
  have targetCast' :
      (signedR1ReducedValue (cppFp2CompletedExpression value) :
          ZMod modulus) =
        cppFp2CompletedExpression value * ifmaRadix *
          Generated.qr_premult * Generated.qr_premult := by
    simpa [modulus] using targetCast
  have equalTimesRadix :
      (selectedReduceMasked high low modulus (reducer1FactorAt index) :
          ZMod modulus) * ifmaRadix =
        (signedR1ReducedValue (cppFp2CompletedExpression value) :
          ZMod modulus) * ifmaRadix := by
    rw [reducer.1, splitCast, completedCast', targetCast']
    ring
  have radixUnit : IsUnit (ifmaRadix : ZMod modulus) := by
    simpa [modulus] using ifma_radix_is_unit_modulus1 index
  have equalCast := radixUnit.mul_right_cancel equalTimesRadix
  exact (ZMod.natCast_eq_natCast_iff _ _ modulus).1 equalCast

def cppFp2CompletedM1Residue
    (high low : Fin 8 → Nat) (index : Fin 8) : Nat :=
  selectedReduceMasked (high index) (low index) (modulus1At index)
    (reducer1FactorAt index)

theorem fp2_max_neg_offset_m1_lane_bound (index : Fin 8) :
    listAt8 Generated.max_neg_offset_m1 index.val < modulus1At index := by
  fin_cases index <;>
    norm_num [listAt8, modulus1At, Generated.max_neg_offset_m1,
      Generated.moduli1]

theorem fp2_completed_m1_wide_fits
    (rawWide high low : Fin 8 → Nat)
    (rawFit : ∀ index,
      rawWide index ≤ fp2AccumulatorBound * modulus1At index *
        modulus1At index)
    (wideSplit : ∀ index, high index * ifmaRadix + low index =
      rawWide index + listAt8 Generated.max_neg_offset_m1 index.val) :
    ∀ index, high index * ifmaRadix + low index ≤
      delayedCombinedTypeBound * modulus1At index * modulus1At index := by
  intro index
  have offsetBound := fp2_max_neg_offset_m1_lane_bound index
  have modulusPositive : 0 < modulus1At index := by
    fin_cases index <;> norm_num [modulus1At, listAt8, Generated.moduli1]
  have offsetSquare :
      listAt8 Generated.max_neg_offset_m1 index.val ≤
        modulus1At index * modulus1At index := by
    nlinarith
  have rawBound := rawFit index
  have typeBound : fp2AccumulatorBound + 1 ≤ delayedCombinedTypeBound := by
    norm_num [fp2AccumulatorBound, fp2PreparedBound,
      delayedCombinedTypeBound, delayedInputBound, delayedR1TypeBound,
      r1ExactAccumulationBound, multiplicationProductBound, standardBound]
  rw [wideSplit index]
  nlinarith

theorem fp2_completed_m1_residue_represents
    (value : Int) (rawWide high low : Fin 8 → Nat)
    (rawRepresents : ∀ index,
      (rawWide index : ZMod (modulus1At index)) =
        value * ifmaRadix * ifmaRadix * Generated.qr_premult *
          Generated.qr_premult)
    (rawFit : ∀ index,
      rawWide index ≤ fp2AccumulatorBound * modulus1At index *
        modulus1At index)
    (wideSplit : ∀ index, high index * ifmaRadix + low index =
      rawWide index + listAt8 Generated.max_neg_offset_m1 index.val) :
    ∀ index, Nat.ModEq (modulus1At index)
      (cppFp2CompletedM1Residue high low index)
      (signedR1ReducedValue (cppFp2CompletedExpression value)) := by
  intro index
  unfold cppFp2CompletedM1Residue
  exact fp2_completed_m1_reducer_represents value (rawWide index)
    (high index) (low index) index (rawRepresents index) (wideSplit index)
    (fp2_completed_m1_wide_fits rawWide high low rawFit wideSplit index)

theorem fp2_completed_m1_residue_fit
    (rawWide high low : Fin 8 → Nat)
    (rawFit : ∀ index,
      rawWide index ≤ fp2AccumulatorBound * modulus1At index *
        modulus1At index)
    (wideSplit : ∀ index, high index * ifmaRadix + low index =
      rawWide index + listAt8 Generated.max_neg_offset_m1 index.val) :
    ∀ index,
      cppFp2CompletedM1Residue high low index < 2 * modulus1At index := by
  intro index
  have reducer := selected_delayed_reducer_contract
    (high index) (low index) (modulus1At index) (reducer1FactorAt index)
    (modulus1_member index) (reducer1_factor_relation index)
    (fp2_completed_m1_wide_fits rawWide high low rawFit wideSplit index)
  exact reducer.2.2

/-- Physical M2 accumulator after `complete()` and batched R1 change-base. -/
def cppFp2CompletedM2CombinedWideValue
    (rawWide m1Residues : Fin 8 → Nat) (column : Fin 8) : Nat :=
  rawWide column + listAt8 Generated.max_neg_offset_m2 column.val +
    matrixProductValue m1Residues Generated.r1_rns_mat column

def cppFp2CompletedM2Residue
    (rawWide m1Residues : Fin 8 → Nat) (column : Fin 8) : Nat :=
  let wide := cppFp2CompletedM2CombinedWideValue rawWide m1Residues column
  selectedReduceMasked (wide / ifmaRadix) (wide % ifmaRadix)
    (modulus2At column) (reducer2FactorAt column)

theorem fp2_max_neg_offset_m2_lane_bound (index : Fin 8) :
    listAt8 Generated.max_neg_offset_m2 index.val < modulus2At index := by
  fin_cases index <;>
    norm_num [listAt8, modulus2At, Generated.max_neg_offset_m2,
      Generated.moduli2]

theorem fp2_completed_m2_combined_wide_lt_25_square
    (rawWide m1Residues : Fin 8 → Nat)
    (rawFit : ∀ index,
      rawWide index ≤ fp2AccumulatorBound * modulus2At index *
        modulus2At index)
    (m1Fit : ∀ index, m1Residues index < 2 * modulus1At index)
    (column : Fin 8) :
    cppFp2CompletedM2CombinedWideValue rawWide m1Residues column <
      25 * modulus2At column * modulus2At column := by
  have matrixBound := r1_matrix_product_lt_16_square
    m1Residues (by simpa [modulus1At] using m1Fit) column
  have matrixBound' :
      matrixProductValue m1Residues Generated.r1_rns_mat column <
        16 * modulus2At column * modulus2At column := by
    simpa [modulus2At, Nat.pow_two, Nat.mul_assoc] using matrixBound
  have offsetBound := fp2_max_neg_offset_m2_lane_bound column
  have modulusPositive : 0 < modulus2At column := by
    fin_cases column <;> norm_num [modulus2At, listAt8, Generated.moduli2]
  have offsetSquare :
      listAt8 Generated.max_neg_offset_m2 column.val ≤
        modulus2At column * modulus2At column := by
    nlinarith
  have rawBound := rawFit column
  norm_num [fp2AccumulatorBound, fp2PreparedBound, standardBound] at rawBound
  unfold cppFp2CompletedM2CombinedWideValue
  nlinarith

theorem fp2_cross_basis_inverse_lane (index : Fin 8) :
    ((Generated.M1 * crossBasisInverse : Nat) :
        ZMod (modulus2At index)) =
      ((1 : Nat) : ZMod (modulus2At index)) := by
  apply (ZMod.natCast_eq_natCast_iff _ (1 : Nat) _).2
  fin_cases index <;>
    norm_num [modulus2At, listAt8, Generated.moduli2, Generated.M1,
      crossBasisInverse, montgomeryRadix, Generated.convert_to_postmult]

/--
The selected batched M2 path composes the completed physical accumulator with
the R1 change-base contribution and returns the R2 input encoding of the same
correction-free representative.
-/
theorem fp2_completed_m2_residue_represents
    (value : Int) (rawWide m1Residues : Fin 8 → Nat)
    (rawRepresents : ∀ index,
      (rawWide index : ZMod (modulus2At index)) =
        value * ifmaRadix * ifmaRadix * crossBasisInverse *
          crossBasisInverse)
    (rawFit : ∀ index,
      rawWide index ≤ fp2AccumulatorBound * modulus2At index *
        modulus2At index)
    (m1Fit : ∀ index, m1Residues index < 2 * modulus1At index)
    (m1Represent : ∀ index,
      Nat.ModEq (modulus1At index) (m1Residues index)
        (signedR1ReducedValue (cppFp2CompletedExpression value))) :
    ∀ column, Nat.ModEq (modulus2At column)
      (cppFp2CompletedM2Residue rawWide m1Residues column)
      (signedR2InputValue (r1SignedValueFromResidues
        (cppFp2CompletedExpression value) m1Residues)) := by
  intro column
  let modulus := modulus2At column
  let completed := cppFp2CompletedExpression value
  let representative := r1SignedValueFromResidues completed m1Residues
  let wide := cppFp2CompletedM2CombinedWideValue rawWide m1Residues column
  have wideBound : wide ≤ delayedCombinedTypeBound * modulus * modulus := by
    have tight := fp2_completed_m2_combined_wide_lt_25_square
      rawWide m1Residues rawFit m1Fit column
    have typeBound : 25 ≤ delayedCombinedTypeBound := by
      norm_num [delayedCombinedTypeBound, delayedInputBound,
        delayedR1TypeBound, r1ExactAccumulationBound,
        multiplicationProductBound, standardBound]
    dsimp only [wide, modulus]
    nlinarith
  have split := Nat.mod_add_div wide ifmaRadix
  have splitReordered :
      (wide / ifmaRadix) * ifmaRadix + wide % ifmaRadix = wide := by
    simpa [Nat.add_comm, Nat.mul_comm] using split
  have reducer := selected_delayed_reducer_contract
    (wide / ifmaRadix) (wide % ifmaRadix) modulus
    (reducer2FactorAt column) (modulus2_member column)
    (reducer2_factor_relation column) (by
      rw [splitReordered]
      exact wideBound)
  have completedCast := fp2_completed_m2_wide_lane_encoding
    value (rawWide column) column (rawRepresents column)
  have matrixCast :
      (matrixProductValue m1Residues Generated.r1_rns_mat column :
          ZMod modulus) = generatedR1LaneValue m1Residues column := by
    unfold matrixProductValue generatedR1LaneValue modulus modulus2At
    push_cast
    rfl
  have rawMatrix := generated_r1_lane_represents_raw_crt m1Residues column
  have combinedCast :
      (wide : ZMod modulus) =
        (completed + Generated.target_modulus *
          r1RawLiftFromResidues m1Residues) * ifmaRadix * ifmaRadix *
            crossBasisInverse * crossBasisInverse := by
    unfold wide cppFp2CompletedM2CombinedWideValue
    push_cast
    rw [completedCast, matrixCast, rawMatrix]
    unfold completed
    rw [r1_postmult_relation]
    unfold r1RawLiftFromResidues
    push_cast
    simp only [montgomeryRadix, ifmaRadix]
    ring_nf
    ac_rfl
  have reconstruction := r1_signed_value_reconstruction
    completed m1Residues (by simpa [completed] using m1Represent)
  have reconstructionCast :
      (representative : ZMod modulus) * Generated.M1 =
        completed + Generated.target_modulus *
          r1RawLiftFromResidues m1Residues := by
    have castEquality := congrArg
      (fun result : Int => (result : ZMod modulus)) reconstruction
    unfold r1SignedValueNumerator at castEquality
    push_cast at castEquality
    simpa [representative] using castEquality
  have inverseCast := fp2_cross_basis_inverse_lane column
  have targetCast := signed_r2_input_lane representative column
  have splitCast := congrArg
    (fun result : Nat => (result : ZMod modulus)) splitReordered
  push_cast at splitCast
  have outputTimesRadix :
      (cppFp2CompletedM2Residue rawWide m1Residues column : ZMod modulus) *
          ifmaRadix =
        (signedR2InputValue representative : ZMod modulus) * ifmaRadix := by
    unfold cppFp2CompletedM2Residue
    dsimp only
    rw [reducer.1, splitCast, combinedCast, targetCast]
    rw [← reconstructionCast]
    change
      (representative : ZMod modulus) * Generated.M1 * ifmaRadix *
          ifmaRadix * crossBasisInverse * crossBasisInverse =
        representative * ifmaRadix * crossBasisInverse * ifmaRadix
    have inverseCast' :
        (Generated.M1 : ZMod modulus) * crossBasisInverse = 1 := by
      simpa only [Nat.cast_mul, Nat.cast_one] using inverseCast
    rw [show
      (representative : ZMod modulus) * Generated.M1 * ifmaRadix *
          ifmaRadix * crossBasisInverse * crossBasisInverse =
        representative * ifmaRadix * ifmaRadix * crossBasisInverse *
          (Generated.M1 * crossBasisInverse) by ring,
      inverseCast']
    ring
  have radixUnit : IsUnit (ifmaRadix : ZMod modulus) := by
    simpa [modulus] using ifma_radix_is_unit_modulus2 column
  have outputCast := radixUnit.mul_right_cancel outputTimesRadix
  exact (ZMod.natCast_eq_natCast_iff _ _ modulus).1 outputCast

theorem fp2_completed_m2_residue_fit
    (rawWide m1Residues : Fin 8 → Nat)
    (rawFit : ∀ index,
      rawWide index ≤ fp2AccumulatorBound * modulus2At index *
        modulus2At index)
    (m1Fit : ∀ index, m1Residues index < 2 * modulus1At index) :
    ∀ column,
      cppFp2CompletedM2Residue rawWide m1Residues column <
        2 * modulus2At column := by
  intro column
  let modulus := modulus2At column
  let wide := cppFp2CompletedM2CombinedWideValue rawWide m1Residues column
  have tight := fp2_completed_m2_combined_wide_lt_25_square
    rawWide m1Residues rawFit m1Fit column
  have wideBound : wide ≤ delayedCombinedTypeBound * modulus * modulus := by
    have typeBound : 25 ≤ delayedCombinedTypeBound := by
      norm_num [delayedCombinedTypeBound, delayedInputBound,
        delayedR1TypeBound, r1ExactAccumulationBound,
        multiplicationProductBound, standardBound]
    dsimp only [wide, modulus]
    nlinarith
  have split := Nat.mod_add_div wide ifmaRadix
  have splitReordered :
      (wide / ifmaRadix) * ifmaRadix + wide % ifmaRadix = wide := by
    simpa [Nat.add_comm, Nat.mul_comm] using split
  have reducer := selected_delayed_reducer_contract
    (wide / ifmaRadix) (wide % ifmaRadix) modulus
    (reducer2FactorAt column) (modulus2_member column)
    (reducer2_factor_relation column) (by
      rw [splitReordered]
      exact wideBound)
  simpa [cppFp2CompletedM2Residue, wide, modulus] using reducer.2.2

theorem fp2_completion_offset_dominates_signed_expression :
    (2 * (fp2BaseRedundantMultiple * Generated.target_modulus) ^ 2 : Int) <
      fp2CompletionOffset := by
  norm_num [fp2CompletionOffset, fp2BaseRedundantMultiple,
    fp2MaxNegativeRns]
  have targetPositive : (0 : Int) < Generated.target_modulus := by
    norm_num [Generated.target_modulus]
  nlinarith [sq_pos_of_pos targetPositive]

theorem fp2_completed_expression_under_five_target_products :
    (fp2CompletionOffset : Int) +
        2 * (fp2BaseRedundantMultiple * Generated.target_modulus) ^ 2 <
      5 * Generated.target_modulus * Generated.M1 := by
  have coefficientBound :
      (3094400 : Int) * Generated.target_modulus < 5 * Generated.M1 := by
    norm_num [Generated.target_modulus, Generated.M1]
  have targetPositive : (0 : Int) < Generated.target_modulus := by
    norm_num [Generated.target_modulus]
  have scaled := mul_lt_mul_of_pos_right coefficientBound targetPositive
  norm_num [fp2CompletionOffset, fp2BaseRedundantMultiple,
    fp2MaxNegativeRns]
  nlinarith

/--
The physical completion offset converts every selected two-product Fp2
coordinate into a nonnegative correction-free R1 representative below `22p`.
This is deliberately about the completed source value, not the wider abstract
signed R1 domain retained in `SignedPipeline`.
-/
theorem fp2_completed_r1_representative_range
    (value : Int) (residues : Fin 8 → Nat)
    (valueBound : |value| <
      2 * (fp2BaseRedundantMultiple * Generated.target_modulus) ^ 2)
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
    have dominates := fp2_completion_offset_dominates_signed_expression
    nlinarith
  have completedUpper :
      cppFp2CompletedExpression value <
        5 * Generated.target_modulus * Generated.M1 := by
    unfold cppFp2CompletedExpression
    have upper := fp2_completed_expression_under_five_target_products
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

theorem fp2_completed_r2_representative_is_exact
    (value : Int) (m1Residues m2Residues : Fin 8 → Nat)
    (valueBound : |value| <
      2 * (fp2BaseRedundantMultiple * Generated.target_modulus) ^ 2)
    (m1ResiduesFit : ∀ index,
      m1Residues index < 2 * listAt8 Generated.moduli1 index.val)
    (m1ResiduesRepresent : ∀ index,
      Nat.ModEq (modulus1At index) (m1Residues index)
        (signedR1ReducedValue (cppFp2CompletedExpression value)))
    (m2ResiduesFit : ∀ index,
      m2Residues index < 2 * listAt8 Generated.moduli2 index.val)
    (m2ResiduesRepresent : ∀ index,
      Nat.ModEq (modulus2At index) (m2Residues index)
        (signedR2InputValue (r1SignedValueFromResidues
          (cppFp2CompletedExpression value) m1Residues))) :
    correctedCrtRepresentative m2Residues r2ShiftedFactor Generated.M2 =
      r1SignedValueFromResidues
        (cppFp2CompletedExpression value) m1Residues := by
  let representative := r1SignedValueFromResidues
    (cppFp2CompletedExpression value) m1Residues
  have representativeRange := fp2_completed_r1_representative_range
    value m1Residues valueBound m1ResiduesFit m1ResiduesRepresent
  have representativeNonnegative : (0 : Int) ≤ representative := by
    simpa [representative] using representativeRange.1
  have representativeBound :
      representative < 22 * Generated.target_modulus := by
    simpa [representative] using representativeRange.2
  have representativeNatBound :
      representative.toNat < 22 * Generated.target_modulus := by
    rw [← Int.ofNat_lt]
    simpa [Int.toNat_of_nonneg representativeNonnegative] using
      representativeBound
  have representedNat : ∀ index,
      Nat.ModEq (listAt8 Generated.moduli2 index.val)
        (m2Residues index) (r2InputValue representative.toNat) := by
    intro index
    simpa [representative, modulus2At,
      signed_r2_input_of_nonnegative representative
        representativeNonnegative] using m2ResiduesRepresent index
  have exact := r2_standard_corrected_representative_is_exact_22
    representative.toNat m2Residues representativeNatBound
      m2ResiduesFit representedNat
  simpa [representative, Int.toNat_of_nonneg representativeNonnegative] using
    exact

/--
Physical completed M2 reduction followed by the fixed R2 correction returns
the exact selected Fp2 representative. The remaining producer obligations are
the raw two-product lane relations and their already modeled `8m²` bounds.
-/
theorem fp2_completed_physical_r2_representative_is_exact
    (value : Int) (rawWide m1Residues : Fin 8 → Nat)
    (valueBound : |value| <
      2 * (fp2BaseRedundantMultiple * Generated.target_modulus) ^ 2)
    (rawRepresents : ∀ index,
      (rawWide index : ZMod (modulus2At index)) =
        value * ifmaRadix * ifmaRadix * crossBasisInverse *
          crossBasisInverse)
    (rawFit : ∀ index,
      rawWide index ≤ fp2AccumulatorBound * modulus2At index *
        modulus2At index)
    (m1Fit : ∀ index, m1Residues index < 2 * modulus1At index)
    (m1Represent : ∀ index,
      Nat.ModEq (modulus1At index) (m1Residues index)
        (signedR1ReducedValue (cppFp2CompletedExpression value))) :
    correctedCrtRepresentative
        (cppFp2CompletedM2Residue rawWide m1Residues)
        r2ShiftedFactor Generated.M2 =
      r1SignedValueFromResidues
        (cppFp2CompletedExpression value) m1Residues := by
  apply fp2_completed_r2_representative_is_exact value m1Residues
    (cppFp2CompletedM2Residue rawWide m1Residues) valueBound m1Fit
    m1Represent
  · exact fp2_completed_m2_residue_fit rawWide m1Residues rawFit m1Fit
  · exact fp2_completed_m2_residue_represents value rawWide m1Residues
      rawRepresents rawFit m1Fit m1Represent

/--
End-to-end selected completion/reduction packet for one Fp2 coordinate. The
inputs are the raw two-product wide lanes before `complete()`; all generated
offset, M1 reduction, R1 change-base, M2 reduction, and R2 correction steps are
then source-shaped definitions rather than assumed residue relations.
-/
theorem fp2_selected_completed_reduction_is_exact
    (value : Int)
    (rawM1 m1High m1Low rawM2 : Fin 8 → Nat)
    (valueBound : |value| <
      2 * (fp2BaseRedundantMultiple * Generated.target_modulus) ^ 2)
    (rawM1Represents : ∀ index,
      (rawM1 index : ZMod (modulus1At index)) =
        value * ifmaRadix * ifmaRadix * Generated.qr_premult *
          Generated.qr_premult)
    (rawM1Fit : ∀ index,
      rawM1 index ≤ fp2AccumulatorBound * modulus1At index *
        modulus1At index)
    (m1WideSplit : ∀ index,
      m1High index * ifmaRadix + m1Low index =
        rawM1 index + listAt8 Generated.max_neg_offset_m1 index.val)
    (rawM2Represents : ∀ index,
      (rawM2 index : ZMod (modulus2At index)) =
        value * ifmaRadix * ifmaRadix * crossBasisInverse *
          crossBasisInverse)
    (rawM2Fit : ∀ index,
      rawM2 index ≤ fp2AccumulatorBound * modulus2At index *
        modulus2At index) :
    let m1Residues := cppFp2CompletedM1Residue m1High m1Low
    correctedCrtRepresentative
        (cppFp2CompletedM2Residue rawM2 m1Residues)
        r2ShiftedFactor Generated.M2 =
      r1SignedValueFromResidues
        (cppFp2CompletedExpression value) m1Residues := by
  dsimp only
  let m1Residues := cppFp2CompletedM1Residue m1High m1Low
  have m1Fit := fp2_completed_m1_residue_fit
    rawM1 m1High m1Low rawM1Fit m1WideSplit
  have m1Represent := fp2_completed_m1_residue_represents value rawM1
    m1High m1Low rawM1Represents rawM1Fit m1WideSplit
  exact fp2_completed_physical_r2_representative_is_exact value rawM2
    m1Residues valueBound rawM2Represents rawM2Fit
    (by simpa [m1Residues] using m1Fit)
    (by simpa [m1Residues] using m1Represent)

end CatBLST.AVX.MatrixNoK
