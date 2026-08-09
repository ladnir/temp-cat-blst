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

private def fp2Sum3 {R : Type*} [Add R]
    (a b c : Fp2Value R) : Fp2Value R :=
  Fp2Value.add (Fp2Value.add a b) c

private def fp2Sum5 {R : Type*} [Add R]
    (a b c d e : Fp2Value R) : Fp2Value R :=
  Fp2Value.add (Fp2Value.add (Fp2Value.add (Fp2Value.add a b) c) d) e

private def fp2Sum6 {R : Type*} [Add R]
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
products. Standard inputs have component bound two, so its absolute product
sum is at most 88 modulus-squared units, below the C++ `ready<200>` boundary.
The version-3 offset bounds are the exact template values selected by
`Fp12ConstantsCandidate`.
-/

def fp12MaximumSignedProducts : Nat := 22
def fp12ProductAccumulatorBound : Nat :=
  fp12MaximumSignedProducts * standardBound * standardBound
def fp12ReadyBound : Nat := 200
def fp12Version3OffsetElementBound : Nat :=
  10 * standardBound * multiplicationProductBound
def fp12Version3OffsetRnsBound : Nat :=
  10 * standardBound * standardBound * 6 * 6

theorem selected_fp12_bound_values :
    fp12ProductAccumulatorBound = 88 ∧
      fp12Version3OffsetElementBound = 80 ∧
      fp12Version3OffsetRnsBound = 1440 := by
  norm_num [fp12ProductAccumulatorBound, fp12MaximumSignedProducts,
    fp12Version3OffsetElementBound, fp12Version3OffsetRnsBound,
    standardBound, multiplicationProductBound]

theorem fp12_product_accumulator_fits_ready :
    fp12ProductAccumulatorBound ≤ fp12ReadyBound := by
  norm_num [fp12ProductAccumulatorBound, fp12MaximumSignedProducts,
    fp12ReadyBound, standardBound]

end CatBLST.AVX.MatrixNoK
