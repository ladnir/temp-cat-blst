/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Fp2
import Mathlib.Tactic

/-!
# Selected MatrixNoK Fp12 construction

This file models the general multiplication and squaring schedules in
`blst/avx/fields/fp12.h`. The selected C++ layout consists of six consecutive
Fp2 coefficients, split into two Fp6 halves:

`[a0, a1, a2, b0, b1, b2]`.

The defining tower is `Fp2[v]/(v^3 - (1 + u))`, followed by
`Fp6[w]/(w^2 - v)`. The explicit `cppFp12FlatMul` and
`cppFp12FlatSquare` definitions mirror the flattened C++ coefficient
schedules. Component reduction is composed from the accepted selected-Fp2
packet and is therefore semantically transparent here.

Frobenius maps and cyclotomic-only operations deliberately remain outside
this packet: they depend on generated constants and narrower input domains.
-/

namespace CatBLST.AVX.MatrixNoK

namespace Fp2Value

variable {R : Type*}

def double [Add R] (a : Fp2Value R) : Fp2Value R :=
  add a a

end Fp2Value

/-- Three Fp2 coefficients in increasing powers of `v`. -/
structure Fp6Value (R : Type*) where
  c0 : Fp2Value R
  c1 : Fp2Value R
  c2 : Fp2Value R
deriving DecidableEq

@[ext]
theorem Fp6Value.ext {R : Type*} {a b : Fp6Value R}
    (c0Equal : a.c0 = b.c0) (c1Equal : a.c1 = b.c1)
    (c2Equal : a.c2 = b.c2) : a = b := by
  cases a
  cases b
  simp_all

namespace Fp6Value

variable {R : Type*}

def add [Add R] (a b : Fp6Value R) : Fp6Value R :=
  ⟨Fp2Value.add a.c0 b.c0, Fp2Value.add a.c1 b.c1,
    Fp2Value.add a.c2 b.c2⟩

/-- Multiplication by `v`, using `v^3 = 1 + u`. -/
def mulByV [Ring R] (a : Fp6Value R) : Fp6Value R :=
  ⟨Fp2Value.mulByNonresidue a.c2, a.c0, a.c1⟩

/-- Multiplication in `Fp2[v]/(v^3 - (1 + u))`. -/
def mul [CommRing R] (a b : Fp6Value R) : Fp6Value R :=
  ⟨Fp2Value.add (Fp2Value.mul a.c0 b.c0)
      (Fp2Value.mulByNonresidue
        (Fp2Value.add (Fp2Value.mul a.c1 b.c2)
          (Fp2Value.mul a.c2 b.c1))),
    Fp2Value.add
      (Fp2Value.add (Fp2Value.mul a.c0 b.c1)
        (Fp2Value.mul a.c1 b.c0))
      (Fp2Value.mulByNonresidue (Fp2Value.mul a.c2 b.c2)),
    Fp2Value.add
      (Fp2Value.add (Fp2Value.mul a.c0 b.c2)
        (Fp2Value.mul a.c1 b.c1))
      (Fp2Value.mul a.c2 b.c0)⟩

end Fp6Value

/-- Two Fp6 coefficients in increasing powers of `w`. -/
structure Fp12Value (R : Type*) where
  a : Fp6Value R
  b : Fp6Value R
deriving DecidableEq

@[ext]
theorem Fp12Value.ext {R : Type*} {a b : Fp12Value R}
    (aEqual : a.a = b.a) (bEqual : a.b = b.b) : a = b := by
  cases a
  cases b
  simp_all

namespace Fp12Value

variable {R : Type*}

/-- Multiplication in `Fp6[w]/(w^2 - v)`. -/
def mul [CommRing R] (x y : Fp12Value R) : Fp12Value R :=
  ⟨Fp6Value.add (Fp6Value.mul x.a y.a)
      (Fp6Value.mulByV (Fp6Value.mul x.b y.b)),
    Fp6Value.add (Fp6Value.mul x.a y.b) (Fp6Value.mul x.b y.a)⟩

def square [CommRing R] (x : Fp12Value R) : Fp12Value R :=
  mul x x

/-- Exact C++ coefficient order, before expanding each Fp2 into `(x,y)`. -/
def flatPairs (x : Fp12Value R) : List (Fp2Value R) :=
  [x.a.c0, x.a.c1, x.a.c2, x.b.c0, x.b.c1, x.b.c2]

/-- Exact twelve-scalar C++ order. -/
def flatCoefficients (x : Fp12Value R) : List R :=
  (flatPairs x).flatMap (fun coefficient => [coefficient.x, coefficient.y])

end Fp12Value

def fp2Sum3 {R : Type*} [Add R]
    (a b c : Fp2Value R) : Fp2Value R :=
  Fp2Value.add (Fp2Value.add a b) c

def fp2Sum5 {R : Type*} [Add R]
    (a b c d e : Fp2Value R) : Fp2Value R :=
  Fp2Value.add (Fp2Value.add (Fp2Value.add (Fp2Value.add a b) c) d) e

def fp2Sum6 {R : Type*} [Add R]
    (a b c d e f : Fp2Value R) : Fp2Value R :=
  Fp2Value.add (fp2Sum5 a b c d e) f

/-!
## Source-shaped flattened schedules
-/

/--
The six Fp2 coefficient formulas obtained by pairing the twelve scalar
polynomials in C++ `fp12_flat_mul`.
-/
def cppFp12FlatMul {R : Type*} [CommRing R]
    (x y : Fp12Value R) : Fp12Value R :=
  ⟨⟨
    Fp2Value.add (Fp2Value.mul x.a.c0 y.a.c0)
      (Fp2Value.mulByNonresidue (fp2Sum5
        (Fp2Value.mul x.a.c1 y.a.c2)
        (Fp2Value.mul x.a.c2 y.a.c1)
        (Fp2Value.mul x.b.c0 y.b.c2)
        (Fp2Value.mul x.b.c1 y.b.c1)
        (Fp2Value.mul x.b.c2 y.b.c0))),
    fp2Sum6
      (Fp2Value.mul x.a.c0 y.a.c1)
      (Fp2Value.mul x.a.c1 y.a.c0)
      (Fp2Value.mulByNonresidue (Fp2Value.mul x.a.c2 y.a.c2))
      (Fp2Value.mul x.b.c0 y.b.c0)
      (Fp2Value.mulByNonresidue (Fp2Value.mul x.b.c1 y.b.c2))
      (Fp2Value.mulByNonresidue (Fp2Value.mul x.b.c2 y.b.c1)),
    fp2Sum6
      (Fp2Value.mul x.a.c0 y.a.c2)
      (Fp2Value.mul x.a.c1 y.a.c1)
      (Fp2Value.mul x.a.c2 y.a.c0)
      (Fp2Value.mul x.b.c0 y.b.c1)
      (Fp2Value.mul x.b.c1 y.b.c0)
      (Fp2Value.mulByNonresidue (Fp2Value.mul x.b.c2 y.b.c2))⟩,
   ⟨
    fp2Sum6
      (Fp2Value.mul x.a.c0 y.b.c0)
      (Fp2Value.mulByNonresidue (Fp2Value.mul x.a.c1 y.b.c2))
      (Fp2Value.mulByNonresidue (Fp2Value.mul x.a.c2 y.b.c1))
      (Fp2Value.mul x.b.c0 y.a.c0)
      (Fp2Value.mulByNonresidue (Fp2Value.mul x.b.c1 y.a.c2))
      (Fp2Value.mulByNonresidue (Fp2Value.mul x.b.c2 y.a.c1)),
    fp2Sum6
      (Fp2Value.mul x.a.c0 y.b.c1)
      (Fp2Value.mul x.a.c1 y.b.c0)
      (Fp2Value.mulByNonresidue (Fp2Value.mul x.a.c2 y.b.c2))
      (Fp2Value.mul x.b.c0 y.a.c1)
      (Fp2Value.mul x.b.c1 y.a.c0)
      (Fp2Value.mulByNonresidue (Fp2Value.mul x.b.c2 y.a.c2)),
    fp2Sum6
      (Fp2Value.mul x.a.c0 y.b.c2)
      (Fp2Value.mul x.a.c1 y.b.c1)
      (Fp2Value.mul x.a.c2 y.b.c0)
      (Fp2Value.mul x.b.c0 y.a.c2)
      (Fp2Value.mul x.b.c1 y.a.c1)
      (Fp2Value.mul x.b.c2 y.a.c0)⟩⟩

/-- The optimized C++ square schedule, expressed as its six Fp2 outputs. -/
def cppFp12FlatSquare {R : Type*} [CommRing R]
    (x : Fp12Value R) : Fp12Value R :=
  ⟨⟨
    fp2Sum3
      (Fp2Value.square x.a.c0)
      (Fp2Value.mulByNonresidue
        (Fp2Value.double (Fp2Value.mul x.a.c1 x.a.c2)))
      (Fp2Value.mulByNonresidue
        (Fp2Value.add
          (Fp2Value.double (Fp2Value.mul x.b.c0 x.b.c2))
          (Fp2Value.square x.b.c1))),
    fp2Sum3
      (Fp2Value.double (Fp2Value.mul x.a.c0 x.a.c1))
      (Fp2Value.mulByNonresidue (Fp2Value.square x.a.c2))
      (Fp2Value.add (Fp2Value.square x.b.c0)
        (Fp2Value.mulByNonresidue
          (Fp2Value.double (Fp2Value.mul x.b.c1 x.b.c2)))),
    fp2Sum3
      (Fp2Value.add
        (Fp2Value.double (Fp2Value.mul x.a.c0 x.a.c2))
        (Fp2Value.square x.a.c1))
      (Fp2Value.double (Fp2Value.mul x.b.c0 x.b.c1))
      (Fp2Value.mulByNonresidue (Fp2Value.square x.b.c2))⟩,
   ⟨
    Fp2Value.double (Fp2Value.add (Fp2Value.mul x.a.c0 x.b.c0)
      (Fp2Value.mulByNonresidue (Fp2Value.add
        (Fp2Value.mul x.a.c1 x.b.c2)
        (Fp2Value.mul x.a.c2 x.b.c1)))),
    Fp2Value.double (Fp2Value.add
      (Fp2Value.add (Fp2Value.mul x.a.c0 x.b.c1)
        (Fp2Value.mul x.a.c1 x.b.c0))
      (Fp2Value.mulByNonresidue (Fp2Value.mul x.a.c2 x.b.c2))),
    Fp2Value.double (fp2Sum3
      (Fp2Value.mul x.a.c0 x.b.c2)
      (Fp2Value.mul x.a.c1 x.b.c1)
      (Fp2Value.mul x.a.c2 x.b.c0))⟩⟩

theorem cpp_fp12_flat_mul_correct {R : Type*} [CommRing R]
    (x y : Fp12Value R) :
    cppFp12FlatMul x y = Fp12Value.mul x y := by
  ext <;>
    simp [cppFp12FlatMul, Fp12Value.mul, Fp6Value.add, Fp6Value.mul,
      Fp6Value.mulByV, fp2Sum5, fp2Sum6, Fp2Value.add,
      Fp2Value.mul, Fp2Value.mulByNonresidue] <;>
    ring

theorem cpp_fp12_flat_square_correct {R : Type*} [CommRing R]
    (x : Fp12Value R) :
    cppFp12FlatSquare x = Fp12Value.square x := by
  ext <;>
    simp [cppFp12FlatSquare, Fp12Value.square, Fp12Value.mul,
      Fp6Value.add, Fp6Value.mul, Fp6Value.mulByV, fp2Sum3,
      Fp2Value.double, Fp2Value.add, Fp2Value.mul, Fp2Value.square,
      Fp2Value.mulByNonresidue] <;>
    ring

theorem cpp_fp12_flat_pair_count {R : Type*} (x : Fp12Value R) :
    (Fp12Value.flatPairs x).length = 6 := by
  rfl

theorem cpp_fp12_flat_coefficient_count {R : Type*} (x : Fp12Value R) :
    (Fp12Value.flatCoefficients x).length = 12 := by
  simp [Fp12Value.flatCoefficients, Fp12Value.flatPairs]

/-!
## Static bound ledger

The largest general multiplication coordinate contains 22 signed scalar
products. A coordinate performs at most ten direct unsigned subtraction units
and at most 22 unsigned addition units. Standard inputs are below `2m`, so the
wide offset occupies less than 81 units and additions contribute at most 88.
The square schedule doubles an offset-bearing expression in four coordinates,
giving the overall physical bound 242. Both paths fit the actual general-field
`MAX_ADD=800` reduction boundary. The offset's lower 80 units dominate all 40
subtraction units. The version-3 offset bounds are the
exact template values selected by `Fp12ConstantsCandidate`.
-/

def fp12MaximumSignedProducts : Nat := 22
def fp12MaximumNegativeProducts : Nat := 10
def fp12MaximumPositiveProducts : Nat := 22
def fp12ProductAccumulatorBound : Nat :=
  fp12MaximumSignedProducts * standardBound * standardBound
def fp12ReadyBound : Nat := 800
def fp12Version3OffsetElementBound : Nat :=
  10 * standardBound * multiplicationProductBound
def fp12Version3OffsetRnsBound : Nat :=
  10 * standardBound * standardBound * 6 * 6
def fp12PhysicalNegativeBound : Nat :=
  fp12MaximumNegativeProducts * standardBound * standardBound
def fp12PhysicalPositiveBound : Nat :=
  fp12MaximumPositiveProducts * standardBound * standardBound
def fp12PhysicalCompletedBound : Nat :=
  242

theorem selected_fp12_bound_values :
    fp12ProductAccumulatorBound = 88 ∧
      fp12Version3OffsetElementBound = 80 ∧
      fp12Version3OffsetRnsBound = 1440 ∧
      fp12PhysicalNegativeBound = 40 ∧
      fp12PhysicalPositiveBound = 88 ∧
      fp12PhysicalCompletedBound = 242 := by
  norm_num [fp12ProductAccumulatorBound, fp12MaximumSignedProducts,
    fp12Version3OffsetElementBound, fp12Version3OffsetRnsBound,
    fp12PhysicalNegativeBound, fp12PhysicalPositiveBound,
    fp12PhysicalCompletedBound, fp12MaximumNegativeProducts,
    fp12MaximumPositiveProducts, standardBound, multiplicationProductBound]

theorem fp12_product_accumulator_fits_ready :
    fp12PhysicalNegativeBound ≤ fp12Version3OffsetElementBound ∧
      fp12PhysicalCompletedBound ≤ fp12ReadyBound := by
  norm_num [fp12PhysicalNegativeBound, fp12Version3OffsetElementBound,
    fp12PhysicalCompletedBound, fp12PhysicalPositiveBound,
    fp12MaximumNegativeProducts, fp12MaximumPositiveProducts,
    fp12ReadyBound, standardBound, multiplicationProductBound]

/-!
## Generated wide offset

The source stores each offset as a high/low radix-`2^52` pair. Extra multiples
of the auxiliary modulus make unsigned subtraction safe; modulo that modulus,
each pair is exactly the wide selected encoding of
`9 * 1440 * p^2`, which is zero in the target field.
-/

def fp12SemanticOffset : Nat :=
  9 * fp12Version3OffsetRnsBound * Generated.target_modulus *
    Generated.target_modulus

def cppFp12OffsetM1Wide (index : Fin 8) : Nat :=
  listAt8 Generated.fp12_offset_m1_high index.val * ifmaRadix +
    listAt8 Generated.fp12_offset_m1_low index.val

def cppFp12OffsetM2Wide (index : Fin 8) : Nat :=
  listAt8 Generated.fp12_offset_m2_high index.val * ifmaRadix +
    listAt8 Generated.fp12_offset_m2_low index.val

theorem fp12_offset_m1_lane_encoding (index : Fin 8) :
    (cppFp12OffsetM1Wide index : ZMod (modulus1At index)) =
      ((fp12SemanticOffset * ifmaRadix * ifmaRadix *
        Generated.qr_premult * Generated.qr_premult : Nat) :
          ZMod (modulus1At index)) := by
  apply (ZMod.natCast_eq_natCast_iff _ _ _).2
  fin_cases index <;>
    norm_num [cppFp12OffsetM1Wide, listAt8, modulus1At,
      Generated.fp12_offset_m1_high, Generated.fp12_offset_m1_low,
      Generated.moduli1, fp12SemanticOffset,
      fp12Version3OffsetRnsBound, standardBound, Generated.target_modulus, ifmaRadix,
      Generated.qr_premult]

theorem fp12_offset_m2_lane_encoding (index : Fin 8) :
    (cppFp12OffsetM2Wide index : ZMod (modulus2At index)) =
      ((fp12SemanticOffset * ifmaRadix * ifmaRadix *
        crossBasisInverse * crossBasisInverse : Nat) :
          ZMod (modulus2At index)) := by
  apply (ZMod.natCast_eq_natCast_iff _ _ _).2
  fin_cases index <;>
    norm_num [cppFp12OffsetM2Wide, listAt8, modulus2At,
      Generated.fp12_offset_m2_high, Generated.fp12_offset_m2_low,
      Generated.moduli2, fp12SemanticOffset,
      fp12Version3OffsetRnsBound, standardBound, Generated.target_modulus, ifmaRadix,
      crossBasisInverse, montgomeryRadix, Generated.convert_to_postmult]

theorem fp12_offset_m1_lane_bounds (index : Fin 8) :
    fp12Version3OffsetElementBound * modulus1At index * modulus1At index ≤
        cppFp12OffsetM1Wide index ∧
      cppFp12OffsetM1Wide index <
        (fp12Version3OffsetElementBound + 1) *
          modulus1At index * modulus1At index := by
  fin_cases index <;>
    norm_num [cppFp12OffsetM1Wide, listAt8, modulus1At,
      Generated.fp12_offset_m1_high, Generated.fp12_offset_m1_low,
      Generated.moduli1, fp12Version3OffsetElementBound, standardBound,
      multiplicationProductBound, ifmaRadix]

theorem fp12_offset_m2_lane_bounds (index : Fin 8) :
    fp12Version3OffsetElementBound * modulus2At index * modulus2At index ≤
        cppFp12OffsetM2Wide index ∧
      cppFp12OffsetM2Wide index <
        (fp12Version3OffsetElementBound + 1) *
          modulus2At index * modulus2At index := by
  fin_cases index <;>
    norm_num [cppFp12OffsetM2Wide, listAt8, modulus2At,
      Generated.fp12_offset_m2_high, Generated.fp12_offset_m2_low,
      Generated.moduli2, fp12Version3OffsetElementBound, standardBound,
      multiplicationProductBound, ifmaRadix]

/-!
## Physical unsigned schedule and M1 reduction

The flattened source evaluates `offset + positive - negative` in unsigned
wide lanes. The generated lower bound proves that subtraction cannot
underflow. These definitions aggregate a single scalar coordinate; the formula
theorems above establish the exact positive and negative terms for all twelve
coordinates.
-/

def cppFp12PhysicalWide
    (offset positive negative : Nat) : Nat :=
  offset + positive - negative

def cppFp12CompletedExpression (value : Int) : Int :=
  value + fp12SemanticOffset

theorem fp12_physical_wide_bounds
    (offset positive negative modulus : Nat)
    (offsetLower : fp12Version3OffsetElementBound * modulus * modulus ≤
      offset)
    (offsetUpper : offset <
      (fp12Version3OffsetElementBound + 1) * modulus * modulus)
    (positiveFit : positive ≤
      fp12PhysicalPositiveBound * modulus * modulus)
    (negativeFit : negative ≤
      fp12PhysicalNegativeBound * modulus * modulus) :
    negative ≤ offset + positive ∧
      cppFp12PhysicalWide offset positive negative <
        fp12PhysicalCompletedBound * modulus * modulus := by
  have negativeDominated : negative ≤ offset := by
    have ledger : fp12PhysicalNegativeBound ≤
        fp12Version3OffsetElementBound := by
      exact fp12_product_accumulator_fits_ready.1
    have scaled := Nat.mul_le_mul_right (modulus * modulus) ledger
    have scaled' : fp12PhysicalNegativeBound * modulus * modulus ≤
        fp12Version3OffsetElementBound * modulus * modulus := by
      simpa [Nat.mul_assoc] using scaled
    exact negativeFit.trans (scaled'.trans offsetLower)
  constructor
  · exact negativeDominated.trans (Nat.le_add_right offset positive)
  · unfold cppFp12PhysicalWide
    have subtractionBound : offset + positive - negative ≤
        offset + positive := Nat.sub_le _ _
    have ledger : fp12Version3OffsetElementBound + 1 +
        fp12PhysicalPositiveBound = 169 := by
      norm_num [fp12Version3OffsetElementBound, fp12PhysicalPositiveBound,
        fp12MaximumPositiveProducts, standardBound,
        multiplicationProductBound]
    rw [show fp12PhysicalCompletedBound = 242 by rfl]
    nlinarith

theorem fp12_physical_m1_lane_encoding
    (value : Int) (positive negative : Nat) (index : Fin 8)
    (positiveValue negativeValue : Int)
    (valueFormula : value = positiveValue - negativeValue)
    (positiveRepresents :
      (positive : ZMod (modulus1At index)) =
        positiveValue * ifmaRadix * ifmaRadix *
          Generated.qr_premult * Generated.qr_premult)
    (negativeRepresents :
      (negative : ZMod (modulus1At index)) =
        negativeValue * ifmaRadix * ifmaRadix *
          Generated.qr_premult * Generated.qr_premult)
    (noUnderflow : negative ≤ cppFp12OffsetM1Wide index + positive) :
    (cppFp12PhysicalWide (cppFp12OffsetM1Wide index) positive negative :
        ZMod (modulus1At index)) =
      cppFp12CompletedExpression value * ifmaRadix * ifmaRadix *
        Generated.qr_premult * Generated.qr_premult := by
  unfold cppFp12PhysicalWide
  rw [Nat.cast_sub noUnderflow, Nat.cast_add, positiveRepresents,
    negativeRepresents, fp12_offset_m1_lane_encoding index]
  unfold cppFp12CompletedExpression
  rw [valueFormula]
  push_cast
  ring

theorem fp12_physical_m2_lane_encoding
    (value : Int) (positive negative : Nat) (index : Fin 8)
    (positiveValue negativeValue : Int)
    (valueFormula : value = positiveValue - negativeValue)
    (positiveRepresents :
      (positive : ZMod (modulus2At index)) =
        positiveValue * ifmaRadix * ifmaRadix *
          crossBasisInverse * crossBasisInverse)
    (negativeRepresents :
      (negative : ZMod (modulus2At index)) =
        negativeValue * ifmaRadix * ifmaRadix *
          crossBasisInverse * crossBasisInverse)
    (noUnderflow : negative ≤ cppFp12OffsetM2Wide index + positive) :
    (cppFp12PhysicalWide (cppFp12OffsetM2Wide index) positive negative :
        ZMod (modulus2At index)) =
      cppFp12CompletedExpression value * ifmaRadix * ifmaRadix *
        crossBasisInverse * crossBasisInverse := by
  unfold cppFp12PhysicalWide
  rw [Nat.cast_sub noUnderflow, Nat.cast_add, positiveRepresents,
    negativeRepresents, fp12_offset_m2_lane_encoding index]
  unfold cppFp12CompletedExpression
  rw [valueFormula]
  push_cast
  ring

def cppFp12M1Residue
    (high low : Fin 8 → Nat) (index : Fin 8) : Nat :=
  selectedReduceMasked (high index) (low index) (modulus1At index)
    (reducer1FactorAt index)

theorem fp12_m1_reducer_represents
    (value : Int) (rawWide high low : Nat) (index : Fin 8)
    (rawRepresents :
      (rawWide : ZMod (modulus1At index)) =
        cppFp12CompletedExpression value * ifmaRadix * ifmaRadix *
          Generated.qr_premult * Generated.qr_premult)
    (wideSplit : high * ifmaRadix + low = rawWide)
    (wideFits : high * ifmaRadix + low ≤
      fp12ReadyBound * modulus1At index * modulus1At index) :
    Nat.ModEq (modulus1At index)
      (selectedReduceMasked high low (modulus1At index)
        (reducer1FactorAt index))
      (signedR1ReducedValue (cppFp12CompletedExpression value)) := by
  let modulus := modulus1At index
  have selectedFits : high * ifmaRadix + low ≤
      delayedCombinedTypeBound * modulus * modulus := by
    have ledger : fp12ReadyBound ≤ delayedCombinedTypeBound := by
      norm_num [fp12ReadyBound, delayedCombinedTypeBound,
        delayedInputBound, delayedR1TypeBound, r1ExactAccumulationBound,
        multiplicationProductBound, standardBound]
    have scaled := Nat.mul_le_mul_right (modulus * modulus) ledger
    have scaled' : fp12ReadyBound * modulus * modulus ≤
        delayedCombinedTypeBound * modulus * modulus := by
      simpa [Nat.mul_assoc] using scaled
    have wideFits' : high * ifmaRadix + low ≤
        fp12ReadyBound * modulus * modulus := by
      simpa [modulus] using wideFits
    exact wideFits'.trans scaled'
  have reducer := selected_delayed_reducer_contract high low modulus
    (reducer1FactorAt index) (modulus1_member index)
    (reducer1_factor_relation index) selectedFits
  have splitCast := congrArg
    (fun result : Nat => (result : ZMod modulus)) wideSplit
  push_cast at splitCast
  have targetCast := signed_r1_reduced_value_lane
    (cppFp12CompletedExpression value) index
  have rawCast : (rawWide : ZMod modulus) =
      cppFp12CompletedExpression value * ifmaRadix * ifmaRadix *
        Generated.qr_premult * Generated.qr_premult := by
    simpa [modulus] using rawRepresents
  have targetCast' :
      (signedR1ReducedValue (cppFp12CompletedExpression value) :
          ZMod modulus) =
        cppFp12CompletedExpression value * ifmaRadix *
          Generated.qr_premult * Generated.qr_premult := by
    simpa [modulus] using targetCast
  have equalTimesRadix :
      (selectedReduceMasked high low modulus (reducer1FactorAt index) :
          ZMod modulus) * ifmaRadix =
        (signedR1ReducedValue (cppFp12CompletedExpression value) :
          ZMod modulus) * ifmaRadix := by
    rw [reducer.1, splitCast, rawCast, targetCast']
    ring
  have radixUnit : IsUnit (ifmaRadix : ZMod modulus) := by
    simpa [modulus] using ifma_radix_is_unit_modulus1 index
  have equalCast := radixUnit.mul_right_cancel equalTimesRadix
  exact (ZMod.natCast_eq_natCast_iff _ _ modulus).1 equalCast

theorem fp12_m1_residues_represent
    (value : Int) (rawWide high low : Fin 8 → Nat)
    (rawRepresents : ∀ index,
      (rawWide index : ZMod (modulus1At index)) =
        cppFp12CompletedExpression value * ifmaRadix * ifmaRadix *
          Generated.qr_premult * Generated.qr_premult)
    (wideSplit : ∀ index,
      high index * ifmaRadix + low index = rawWide index)
    (wideFits : ∀ index, rawWide index ≤
      fp12ReadyBound * modulus1At index * modulus1At index) :
    ∀ index, Nat.ModEq (modulus1At index)
      (cppFp12M1Residue high low index)
      (signedR1ReducedValue (cppFp12CompletedExpression value)) := by
  intro index
  unfold cppFp12M1Residue
  exact fp12_m1_reducer_represents value (rawWide index) (high index)
    (low index) index (rawRepresents index) (wideSplit index)
    (by rw [wideSplit index]; exact wideFits index)

theorem fp12_m1_residues_fit
    (rawWide high low : Fin 8 → Nat)
    (wideSplit : ∀ index,
      high index * ifmaRadix + low index = rawWide index)
    (wideFits : ∀ index, rawWide index ≤
      fp12ReadyBound * modulus1At index * modulus1At index) :
    ∀ index, cppFp12M1Residue high low index <
      2 * modulus1At index := by
  intro index
  have ledger : fp12ReadyBound ≤ delayedCombinedTypeBound := by
    norm_num [fp12ReadyBound, delayedCombinedTypeBound,
      delayedInputBound, delayedR1TypeBound, r1ExactAccumulationBound,
      multiplicationProductBound, standardBound]
  have scaled := Nat.mul_le_mul_right
    (modulus1At index * modulus1At index) ledger
  have selectedFits : high index * ifmaRadix + low index ≤
      delayedCombinedTypeBound * modulus1At index * modulus1At index := by
    rw [wideSplit index]
    exact (wideFits index).trans (by
      simpa [Nat.mul_assoc] using scaled)
  have reducer := selected_delayed_reducer_contract
    (high index) (low index) (modulus1At index) (reducer1FactorAt index)
    (modulus1_member index) (reducer1_factor_relation index) selectedFits
  simpa [cppFp12M1Residue] using reducer.2.2

theorem fp12_m1_residues_represent_completed
    (completed : Int) (rawWide high low : Fin 8 → Nat)
    (rawRepresents : ∀ index,
      (rawWide index : ZMod (modulus1At index)) =
        completed * ifmaRadix * ifmaRadix *
          Generated.qr_premult * Generated.qr_premult)
    (wideSplit : ∀ index,
      high index * ifmaRadix + low index = rawWide index)
    (wideFits : ∀ index, rawWide index ≤
      fp12ReadyBound * modulus1At index * modulus1At index) :
    ∀ index, Nat.ModEq (modulus1At index)
      (cppFp12M1Residue high low index)
      (signedR1ReducedValue completed) := by
  let value := completed - fp12SemanticOffset
  have represented := fp12_m1_residues_represent value rawWide high low
    (by
      intro index
      simpa [value, cppFp12CompletedExpression] using rawRepresents index)
    wideSplit wideFits
  simpa [value, cppFp12CompletedExpression] using represented

def cppFp12M2CombinedWide
    (rawWide m1Residues : Fin 8 → Nat) (column : Fin 8) : Nat :=
  rawWide column +
    matrixProductValue m1Residues Generated.r1_rns_mat column

def cppFp12M2Residue
    (rawWide m1Residues : Fin 8 → Nat) (column : Fin 8) : Nat :=
  let wide := cppFp12M2CombinedWide rawWide m1Residues column
  selectedReduceMasked (wide / ifmaRadix) (wide % ifmaRadix)
    (modulus2At column) (reducer2FactorAt column)

theorem fp12_m2_combined_wide_lt_258_square
    (rawWide m1Residues : Fin 8 → Nat)
    (rawFit : ∀ index, rawWide index <
      fp12PhysicalCompletedBound * modulus2At index * modulus2At index)
    (m1Fit : ∀ index, m1Residues index < 2 * modulus1At index)
    (column : Fin 8) :
    cppFp12M2CombinedWide rawWide m1Residues column <
      258 * modulus2At column * modulus2At column := by
  have matrixBound := r1_matrix_product_lt_16_square
    m1Residues (by simpa [modulus1At] using m1Fit) column
  have matrixBound' :
      matrixProductValue m1Residues Generated.r1_rns_mat column <
        16 * modulus2At column * modulus2At column := by
    simpa [modulus2At, Nat.pow_two, Nat.mul_assoc] using matrixBound
  have rawBound := rawFit column
  have physicalLedger : fp12PhysicalCompletedBound = 242 := by
    norm_num [fp12PhysicalCompletedBound,
      fp12Version3OffsetElementBound, fp12PhysicalPositiveBound,
      fp12MaximumPositiveProducts, standardBound,
      multiplicationProductBound]
  rw [physicalLedger] at rawBound
  unfold cppFp12M2CombinedWide
  nlinarith

theorem fp12_m2_residue_represents
    (value : Int) (rawWide m1Residues : Fin 8 → Nat)
    (rawRepresents : ∀ index,
      (rawWide index : ZMod (modulus2At index)) =
        cppFp12CompletedExpression value * ifmaRadix * ifmaRadix *
          crossBasisInverse * crossBasisInverse)
    (rawFit : ∀ index, rawWide index <
      fp12PhysicalCompletedBound * modulus2At index * modulus2At index)
    (m1Fit : ∀ index, m1Residues index < 2 * modulus1At index)
    (m1Represent : ∀ index,
      Nat.ModEq (modulus1At index) (m1Residues index)
        (signedR1ReducedValue (cppFp12CompletedExpression value))) :
    ∀ column, Nat.ModEq (modulus2At column)
      (cppFp12M2Residue rawWide m1Residues column)
      (signedR2InputValue (r1SignedValueFromResidues
        (cppFp12CompletedExpression value) m1Residues)) := by
  intro column
  let modulus := modulus2At column
  let completed := cppFp12CompletedExpression value
  let representative := r1SignedValueFromResidues completed m1Residues
  let wide := cppFp12M2CombinedWide rawWide m1Residues column
  have wideBound : wide ≤ delayedCombinedTypeBound * modulus * modulus := by
    have tight := fp12_m2_combined_wide_lt_258_square
      rawWide m1Residues rawFit m1Fit column
    have typeBound : 258 ≤ delayedCombinedTypeBound := by
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
    unfold wide cppFp12M2CombinedWide
    push_cast
    rw [rawRepresents column, matrixCast, rawMatrix]
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
      (cppFp12M2Residue rawWide m1Residues column : ZMod modulus) *
          ifmaRadix =
        (signedR2InputValue representative : ZMod modulus) *
          ifmaRadix := by
    unfold cppFp12M2Residue
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

theorem fp12_m2_residue_fit
    (rawWide m1Residues : Fin 8 → Nat)
    (rawFit : ∀ index, rawWide index <
      fp12PhysicalCompletedBound * modulus2At index * modulus2At index)
    (m1Fit : ∀ index, m1Residues index < 2 * modulus1At index) :
    ∀ column, cppFp12M2Residue rawWide m1Residues column <
      2 * modulus2At column := by
  intro column
  let modulus := modulus2At column
  let wide := cppFp12M2CombinedWide rawWide m1Residues column
  have tight := fp12_m2_combined_wide_lt_258_square
    rawWide m1Residues rawFit m1Fit column
  have wideBound : wide ≤ delayedCombinedTypeBound * modulus * modulus := by
    have typeBound : 258 ≤ delayedCombinedTypeBound := by
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
  simpa [cppFp12M2Residue, wide, modulus] using reducer.2.2

theorem fp12_m2_residue_represents_completed
    (completed : Int) (rawWide m1Residues : Fin 8 → Nat)
    (rawRepresents : ∀ index,
      (rawWide index : ZMod (modulus2At index)) =
        completed * ifmaRadix * ifmaRadix *
          crossBasisInverse * crossBasisInverse)
    (rawFit : ∀ index, rawWide index <
      fp12PhysicalCompletedBound * modulus2At index * modulus2At index)
    (m1Fit : ∀ index, m1Residues index < 2 * modulus1At index)
    (m1Represent : ∀ index,
      Nat.ModEq (modulus1At index) (m1Residues index)
        (signedR1ReducedValue completed)) :
    ∀ column, Nat.ModEq (modulus2At column)
      (cppFp12M2Residue rawWide m1Residues column)
      (signedR2InputValue
        (r1SignedValueFromResidues completed m1Residues)) := by
  let value := completed - fp12SemanticOffset
  have represented := fp12_m2_residue_represents value rawWide m1Residues
    (by
      intro index
      simpa [value, cppFp12CompletedExpression] using rawRepresents index)
    rawFit m1Fit (by
      intro index
      simpa [value, cppFp12CompletedExpression] using m1Represent index)
  simpa [value, cppFp12CompletedExpression] using represented

/-!
## Completed representative range and R2 correction
-/

def fp12NegativeRnsExpressionBound : Nat :=
  fp12MaximumNegativeProducts * standardBound * standardBound * 6 * 6

def fp12PositiveRnsExpressionBound : Nat :=
  fp12MaximumPositiveProducts * standardBound * standardBound * 6 * 6

theorem selected_fp12_rns_expression_bounds :
    fp12NegativeRnsExpressionBound = 1440 ∧
      fp12PositiveRnsExpressionBound = 3168 := by
  norm_num [fp12NegativeRnsExpressionBound,
    fp12PositiveRnsExpressionBound, fp12MaximumNegativeProducts,
    fp12MaximumPositiveProducts, standardBound]

theorem fp12_completed_expression_range
    (value : Int)
    (valueLower :
      -(fp12NegativeRnsExpressionBound *
        Generated.target_modulus * Generated.target_modulus : Nat) ≤ value)
    (valueUpper : value <
      fp12PositiveRnsExpressionBound *
        Generated.target_modulus * Generated.target_modulus) :
    0 ≤ cppFp12CompletedExpression value ∧
      cppFp12CompletedExpression value <
        Generated.target_modulus * Generated.M1 := by
  have negativeLedger : fp12NegativeRnsExpressionBound = 1440 := by
    exact selected_fp12_rns_expression_bounds.1
  have positiveLedger : fp12PositiveRnsExpressionBound = 3168 := by
    exact selected_fp12_rns_expression_bounds.2
  have targetPositive : (0 : Int) < Generated.target_modulus := by
    norm_num [Generated.target_modulus]
  have offsetDominates :
      (fp12NegativeRnsExpressionBound * Generated.target_modulus *
          Generated.target_modulus : Nat) ≤ fp12SemanticOffset := by
    rw [negativeLedger]
    norm_num [fp12SemanticOffset, fp12Version3OffsetRnsBound,
      standardBound]
    nlinarith [Nat.zero_le
      (Generated.target_modulus * Generated.target_modulus)]
  have completedUpper :
      (fp12SemanticOffset : Int) +
          fp12PositiveRnsExpressionBound * Generated.target_modulus *
            Generated.target_modulus <
        Generated.target_modulus * Generated.M1 := by
    rw [positiveLedger]
    norm_num [fp12SemanticOffset, fp12Version3OffsetRnsBound,
      standardBound, Generated.target_modulus, Generated.M1]
  unfold cppFp12CompletedExpression
  constructor
  · have offsetDominatesInt :
        (fp12NegativeRnsExpressionBound * Generated.target_modulus *
          Generated.target_modulus : Nat) ≤ (fp12SemanticOffset : Int) := by
      exact_mod_cast offsetDominates
    nlinarith
  · nlinarith

theorem fp12_completed_r1_representative_range
    (value : Int) (residues : Fin 8 → Nat)
    (valueLower :
      -(fp12NegativeRnsExpressionBound *
        Generated.target_modulus * Generated.target_modulus : Nat) ≤ value)
    (valueUpper : value <
      fp12PositiveRnsExpressionBound *
        Generated.target_modulus * Generated.target_modulus)
    (residuesFit : ∀ index,
      residues index < 2 * listAt8 Generated.moduli1 index.val)
    (residuesRepresent : ∀ index,
      Nat.ModEq (modulus1At index) (residues index)
        (signedR1ReducedValue (cppFp12CompletedExpression value))) :
    0 ≤ r1SignedValueFromResidues
        (cppFp12CompletedExpression value) residues ∧
      r1SignedValueFromResidues
          (cppFp12CompletedExpression value) residues <
        18 * Generated.target_modulus := by
  have completedRange := fp12_completed_expression_range
    value valueLower valueUpper
  have liftBoundNat := r1_standard_raw_lift_lt_17_product residues
    residuesFit
  have liftBound :
      (r1RawLiftFromResidues residues : Int) < 17 * Generated.M1 := by
    exact_mod_cast liftBoundNat
  have liftNonnegative :
      (0 : Int) ≤ r1RawLiftFromResidues residues := by
    exact_mod_cast Nat.zero_le (r1RawLiftFromResidues residues)
  have reconstruction := r1_signed_value_reconstruction
    (cppFp12CompletedExpression value) residues residuesRepresent
  have modulusPositive : (0 : Int) < Generated.M1 := by
    norm_num [Generated.M1]
  unfold r1SignedValueNumerator at reconstruction
  constructor <;> nlinarith

theorem fp12_completed_r2_representative_is_exact
    (value : Int) (m1Residues m2Residues : Fin 8 → Nat)
    (valueLower :
      -(fp12NegativeRnsExpressionBound *
        Generated.target_modulus * Generated.target_modulus : Nat) ≤ value)
    (valueUpper : value <
      fp12PositiveRnsExpressionBound *
        Generated.target_modulus * Generated.target_modulus)
    (m1ResiduesFit : ∀ index,
      m1Residues index < 2 * listAt8 Generated.moduli1 index.val)
    (m1ResiduesRepresent : ∀ index,
      Nat.ModEq (modulus1At index) (m1Residues index)
        (signedR1ReducedValue (cppFp12CompletedExpression value)))
    (m2ResiduesFit : ∀ index,
      m2Residues index < 2 * listAt8 Generated.moduli2 index.val)
    (m2ResiduesRepresent : ∀ index,
      Nat.ModEq (modulus2At index) (m2Residues index)
        (signedR2InputValue (r1SignedValueFromResidues
          (cppFp12CompletedExpression value) m1Residues))) :
    correctedCrtRepresentative m2Residues r2ShiftedFactor Generated.M2 =
      r1SignedValueFromResidues
        (cppFp12CompletedExpression value) m1Residues := by
  let representative := r1SignedValueFromResidues
    (cppFp12CompletedExpression value) m1Residues
  have representativeRange := fp12_completed_r1_representative_range
    value m1Residues valueLower valueUpper m1ResiduesFit
      m1ResiduesRepresent
  have representativeNonnegative : (0 : Int) ≤ representative := by
    simpa [representative] using representativeRange.1
  have representativeBound :
      representative < 18 * Generated.target_modulus := by
    simpa [representative] using representativeRange.2
  have representativeNatBound :
      representative.toNat < 22 * Generated.target_modulus := by
    rw [← Int.ofNat_lt]
    simp [Int.toNat_of_nonneg representativeNonnegative]
    nlinarith
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

def fp12MaximumCompletedRnsBound : Nat := 5760

theorem fp12_maximum_completed_value_below_r1_dynamic_range :
    fp12MaximumCompletedRnsBound * Generated.target_modulus *
        Generated.target_modulus <
      Generated.target_modulus * Generated.M1 := by
  norm_num [fp12MaximumCompletedRnsBound, Generated.target_modulus,
    Generated.M1]

theorem fp12_physical_r1_representative_range
    (completed : Int) (residues : Fin 8 → Nat)
    (completedNonnegative : 0 ≤ completed)
    (completedBound : completed <
      fp12MaximumCompletedRnsBound * Generated.target_modulus *
        Generated.target_modulus)
    (residuesFit : ∀ index,
      residues index < 2 * listAt8 Generated.moduli1 index.val)
    (residuesRepresent : ∀ index,
      Nat.ModEq (modulus1At index) (residues index)
        (signedR1ReducedValue completed)) :
    0 ≤ r1SignedValueFromResidues completed residues ∧
      r1SignedValueFromResidues completed residues <
        18 * Generated.target_modulus := by
  have dynamicRangeNat := fp12_maximum_completed_value_below_r1_dynamic_range
  have dynamicRange :
      (fp12MaximumCompletedRnsBound * Generated.target_modulus *
          Generated.target_modulus : Nat) <
        (Generated.target_modulus * Generated.M1 : Int) := by
    exact_mod_cast dynamicRangeNat
  have completedUpper : completed <
      Generated.target_modulus * Generated.M1 :=
    completedBound.trans dynamicRange
  have liftBoundNat := r1_standard_raw_lift_lt_17_product residues
    residuesFit
  have liftBound :
      (r1RawLiftFromResidues residues : Int) < 17 * Generated.M1 := by
    exact_mod_cast liftBoundNat
  have liftNonnegative :
      (0 : Int) ≤ r1RawLiftFromResidues residues := by
    exact_mod_cast Nat.zero_le (r1RawLiftFromResidues residues)
  have reconstruction := r1_signed_value_reconstruction
    completed residues residuesRepresent
  have modulusPositive : (0 : Int) < Generated.M1 := by
    norm_num [Generated.M1]
  unfold r1SignedValueNumerator at reconstruction
  constructor <;> nlinarith

theorem fp12_physical_r2_representative_is_exact
    (completed : Int) (m1Residues m2Residues : Fin 8 → Nat)
    (completedNonnegative : 0 ≤ completed)
    (completedBound : completed <
      fp12MaximumCompletedRnsBound * Generated.target_modulus *
        Generated.target_modulus)
    (m1ResiduesFit : ∀ index,
      m1Residues index < 2 * listAt8 Generated.moduli1 index.val)
    (m1ResiduesRepresent : ∀ index,
      Nat.ModEq (modulus1At index) (m1Residues index)
        (signedR1ReducedValue completed))
    (m2ResiduesFit : ∀ index,
      m2Residues index < 2 * listAt8 Generated.moduli2 index.val)
    (m2ResiduesRepresent : ∀ index,
      Nat.ModEq (modulus2At index) (m2Residues index)
        (signedR2InputValue
          (r1SignedValueFromResidues completed m1Residues))) :
    correctedCrtRepresentative m2Residues r2ShiftedFactor Generated.M2 =
      r1SignedValueFromResidues completed m1Residues := by
  let representative := r1SignedValueFromResidues completed m1Residues
  have representativeRange := fp12_physical_r1_representative_range
    completed m1Residues completedNonnegative completedBound
      m1ResiduesFit m1ResiduesRepresent
  have representativeNonnegative : (0 : Int) ≤ representative := by
    simpa [representative] using representativeRange.1
  have representativeBound :
      representative < 18 * Generated.target_modulus := by
    simpa [representative] using representativeRange.2
  have representativeNatBound :
      representative.toNat < 22 * Generated.target_modulus := by
    rw [← Int.ofNat_lt]
    simp [Int.toNat_of_nonneg representativeNonnegative]
    nlinarith
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

/-!
End-to-end reduction for one flattened Fp12 scalar coordinate. The formula
packet supplies `value`; the physical unsigned packet supplies the raw M1/M2
wide lanes and their `242m^2` bounds. This theorem then composes M1 reduction,
R1 change-base, M2 reduction, and fixed R2 correction without assumed
intermediate residues.
-/
theorem fp12_selected_physical_reduction_is_exact
    (value : Int)
    (rawM1 m1High m1Low rawM2 : Fin 8 → Nat)
    (valueLower :
      -(fp12NegativeRnsExpressionBound *
        Generated.target_modulus * Generated.target_modulus : Nat) ≤ value)
    (valueUpper : value <
      fp12PositiveRnsExpressionBound *
        Generated.target_modulus * Generated.target_modulus)
    (rawM1Represents : ∀ index,
      (rawM1 index : ZMod (modulus1At index)) =
        cppFp12CompletedExpression value * ifmaRadix * ifmaRadix *
          Generated.qr_premult * Generated.qr_premult)
    (rawM1Fit : ∀ index, rawM1 index ≤
      fp12ReadyBound * modulus1At index * modulus1At index)
    (m1WideSplit : ∀ index,
      m1High index * ifmaRadix + m1Low index = rawM1 index)
    (rawM2Represents : ∀ index,
      (rawM2 index : ZMod (modulus2At index)) =
        cppFp12CompletedExpression value * ifmaRadix * ifmaRadix *
          crossBasisInverse * crossBasisInverse)
    (rawM2Fit : ∀ index, rawM2 index <
      fp12PhysicalCompletedBound * modulus2At index * modulus2At index) :
    let m1Residues := cppFp12M1Residue m1High m1Low
    correctedCrtRepresentative
        (cppFp12M2Residue rawM2 m1Residues)
        r2ShiftedFactor Generated.M2 =
      r1SignedValueFromResidues
        (cppFp12CompletedExpression value) m1Residues := by
  dsimp only
  let m1Residues := cppFp12M1Residue m1High m1Low
  have m1Fit := fp12_m1_residues_fit rawM1 m1High m1Low
    m1WideSplit rawM1Fit
  have m1Represent := fp12_m1_residues_represent value rawM1
    m1High m1Low rawM1Represents m1WideSplit rawM1Fit
  have m2Fit := fp12_m2_residue_fit rawM2 m1Residues rawM2Fit
    (by simpa [m1Residues, modulus1At] using m1Fit)
  have m2Represent := fp12_m2_residue_represents value rawM2
    m1Residues rawM2Represents rawM2Fit
    (by simpa [m1Residues, modulus1At] using m1Fit)
    (by simpa [m1Residues] using m1Represent)
  exact fp12_completed_r2_representative_is_exact value m1Residues
    (cppFp12M2Residue rawM2 m1Residues) valueLower valueUpper
    (by simpa [m1Residues, modulus1At] using m1Fit)
    (by simpa [m1Residues] using m1Represent)
    (by simpa [modulus2At] using m2Fit)
    (by simpa [m1Residues] using m2Represent)

theorem fp12_selected_redundant_reduction_is_exact
    (completed : Int)
    (rawM1 m1High m1Low rawM2 : Fin 8 → Nat)
    (completedNonnegative : 0 ≤ completed)
    (completedBound : completed <
      fp12MaximumCompletedRnsBound * Generated.target_modulus *
        Generated.target_modulus)
    (rawM1Represents : ∀ index,
      (rawM1 index : ZMod (modulus1At index)) =
        completed * ifmaRadix * ifmaRadix *
          Generated.qr_premult * Generated.qr_premult)
    (rawM1Fit : ∀ index, rawM1 index ≤
      fp12ReadyBound * modulus1At index * modulus1At index)
    (m1WideSplit : ∀ index,
      m1High index * ifmaRadix + m1Low index = rawM1 index)
    (rawM2Represents : ∀ index,
      (rawM2 index : ZMod (modulus2At index)) =
        completed * ifmaRadix * ifmaRadix *
          crossBasisInverse * crossBasisInverse)
    (rawM2Fit : ∀ index, rawM2 index <
      fp12PhysicalCompletedBound * modulus2At index * modulus2At index) :
    let m1Residues := cppFp12M1Residue m1High m1Low
    correctedCrtRepresentative
        (cppFp12M2Residue rawM2 m1Residues)
        r2ShiftedFactor Generated.M2 =
      r1SignedValueFromResidues completed m1Residues := by
  dsimp only
  let m1Residues := cppFp12M1Residue m1High m1Low
  have m1Fit := fp12_m1_residues_fit rawM1 m1High m1Low
    m1WideSplit rawM1Fit
  have m1Represent := fp12_m1_residues_represent_completed completed
    rawM1 m1High m1Low rawM1Represents m1WideSplit rawM1Fit
  have m2Fit := fp12_m2_residue_fit rawM2 m1Residues rawM2Fit
    (by simpa [m1Residues, modulus1At] using m1Fit)
  have m2Represent := fp12_m2_residue_represents_completed completed
    rawM2 m1Residues rawM2Represents rawM2Fit
    (by simpa [m1Residues, modulus1At] using m1Fit)
    (by simpa [m1Residues] using m1Represent)
  exact fp12_physical_r2_representative_is_exact completed m1Residues
    (cppFp12M2Residue rawM2 m1Residues) completedNonnegative completedBound
    (by simpa [m1Residues, modulus1At] using m1Fit)
    (by simpa [m1Residues] using m1Represent)
    (by simpa [modulus2At] using m2Fit)
    (by simpa [m1Residues] using m2Represent)

end CatBLST.AVX.MatrixNoK
