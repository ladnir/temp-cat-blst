/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Fp12Frobenius3
import CatBLST.AVX.MatrixNoK.Fp12Coordinates
import Mathlib.Tactic

set_option maxHeartbeats 800000
set_option maxRecDepth 10000

/-!
# Selected MatrixNoK cyclotomic square

This file transcribes the two no-inline halves selected by
the selected split schedule. The semantic model retains the
pre-reduction polynomials, their respective factors of three and six, and the
source `concat` order.

The operation is deliberately not a general Fp12 square.  Its exact domain is
the same Frobenius relation tested by scalar `Fp12::is_cyclotomic`:

`Frob^2(x) = Frob^4(x) * x`.

Zero satisfies this relation as it does in the production predicate.  GT is a
strictly narrower, nonzero prime-order subgroup and is not assumed here.
-/

namespace CatBLST.AVX.MatrixNoK

open Fp12Concrete

/-- Canonical field value loaded by `Fp12Constants::one_third()`. -/
def cyclotomicOneThird : SelectedBaseField :=
  0x11560bf17baa99bc32126fced787c88f984f87adf7ae0c7f9a208c6b4f20a4181472aaa9cb8d555526a9ffffffffc71d

/-- Canonical field value loaded by `Fp12Constants::neg_two_thirds()`. -/
def cyclotomicNegTwoThird : SelectedBaseField :=
  0x11560bf17baa99bc32126fced787c88f984f87adf7ae0c7f9a208c6b4f20a4181472aaa9cb8d555526a9ffffffffc71c

theorem cyclotomic_fraction_facts :
    (3 : SelectedBaseField) * cyclotomicOneThird = 1 ∧
      cyclotomicNegTwoThird = -(2 : SelectedBaseField) * cyclotomicOneThird := by
  native_decide

/-- The six unscaled polynomials in the first selected no-inline half. -/
def cppFp12CyclotomicFirstPolynomials
    (x : FlatFp12 SelectedBaseField) : Fin 6 → SelectedBaseField :=
  let minusX1 := -x 1
  let minusX3 := -x 3
  let minusX5 := -x 5
  let minusX7 := -x 7
  let minusX9 := -x 9
  let minusX11 := -x 11
  let shared0 := x 8 * x 8 + x 9 * minusX9
  let shared2 := x 4 * x 4 + x 5 * minusX5
  let shared4 := x 10 * x 10 + x 11 * minusX11
  ![
    shared0 + (x 8 * minusX9) * 2 + x 0 * x 0 +
      x 0 * cyclotomicNegTwoThird + x 1 * minusX1,
    shared0 + (x 0 * x 1 + x 8 * x 9) * 2 +
      x 1 * cyclotomicNegTwoThird,
    shared2 + (x 4 * minusX5) * 2 +
      x 2 * cyclotomicNegTwoThird + x 6 * x 6 + x 7 * minusX7,
    shared2 + (x 6 * x 7 + x 4 * x 5) * 2 +
      x 3 * cyclotomicNegTwoThird,
    shared4 + (x 10 * minusX11) * 2 + x 2 * x 2 +
      x 3 * minusX3 + x 4 * cyclotomicNegTwoThird,
    shared4 + (x 2 * x 3 + x 10 * x 11) * 2 +
      x 5 * cyclotomicNegTwoThird]

/-- The six unscaled polynomials in the second selected no-inline half. -/
def cppFp12CyclotomicSecondPolynomials
    (x : FlatFp12 SelectedBaseField) : Fin 6 → SelectedBaseField :=
  let minusX7 := -x 7
  let minusX9 := -x 9
  let minusX11 := -x 11
  let sharedPositive := x 2 * x 10 + x 3 * minusX11
  let sharedNegative := x 2 * x 11 + x 3 * x 10
  ![
    sharedPositive - sharedNegative + x 6 * cyclotomicOneThird,
    sharedPositive + sharedNegative + x 7 * cyclotomicOneThird,
    x 0 * x 8 + x 1 * minusX9 + x 8 * cyclotomicOneThird,
    x 0 * x 9 + x 1 * x 8 + x 9 * cyclotomicOneThird,
    x 4 * x 6 + x 5 * minusX7 + x 10 * cyclotomicOneThird,
    x 4 * x 7 + x 5 * x 6 + x 11 * cyclotomicOneThird]

/-- Field semantics of `fp12_flat_cyclotomic_sqr_first_split`. -/
def cppFp12CyclotomicFirstSplit
    (x : FlatFp12 SelectedBaseField) : Fin 6 → SelectedBaseField :=
  fun i => 3 * cppFp12CyclotomicFirstPolynomials x i

/-- Field semantics of `fp12_flat_cyclotomic_sqr_second_split`. -/
def cppFp12CyclotomicSecondSplit
    (x : FlatFp12 SelectedBaseField) : Fin 6 → SelectedBaseField :=
  fun i => 6 * cppFp12CyclotomicSecondPolynomials x i

/-- Exact coordinate order produced by the two split kernels in their final owner. -/
def cppFp12CyclotomicSquareFlat
    (x : FlatFp12 SelectedBaseField) : FlatFp12 SelectedBaseField :=
  let first := cppFp12CyclotomicFirstSplit x
  let second := cppFp12CyclotomicSecondSplit x
  ![first 0, first 1, first 2, first 3, first 4, first 5,
    second 0, second 1, second 2, second 3, second 4, second 5]

def cppFp12CyclotomicSquare (x : SelectedFp12) : SelectedFp12 :=
  packFlatFp12 (cppFp12CyclotomicSquareFlat (unpackFlatFp12 x))

/-!
The direct structure views below are definitionally equal to the flattened
source transcription.  They keep later coordinate proofs small without
replacing or abstracting the source-shaped definitions above.
-/

def cppCyclotomicP0 (x : SelectedFp12) : SelectedBaseField :=
  let minusX1 := -x.a.c0.y
  let minusX9 := -x.b.c1.y
  let shared := x.b.c1.x * x.b.c1.x + x.b.c1.y * minusX9
  shared + (x.b.c1.x * minusX9) * 2 + x.a.c0.x * x.a.c0.x +
    x.a.c0.x * cyclotomicNegTwoThird + x.a.c0.y * minusX1

def cppCyclotomicP1 (x : SelectedFp12) : SelectedBaseField :=
  let minusX9 := -x.b.c1.y
  let shared := x.b.c1.x * x.b.c1.x + x.b.c1.y * minusX9
  shared + (x.a.c0.x * x.a.c0.y + x.b.c1.x * x.b.c1.y) * 2 +
    x.a.c0.y * cyclotomicNegTwoThird

def cppCyclotomicP2 (x : SelectedFp12) : SelectedBaseField :=
  let minusX5 := -x.a.c2.y
  let minusX7 := -x.b.c0.y
  let shared := x.a.c2.x * x.a.c2.x + x.a.c2.y * minusX5
  shared + (x.a.c2.x * minusX5) * 2 +
    x.a.c1.x * cyclotomicNegTwoThird +
    x.b.c0.x * x.b.c0.x + x.b.c0.y * minusX7

def cppCyclotomicP3 (x : SelectedFp12) : SelectedBaseField :=
  let minusX5 := -x.a.c2.y
  let shared := x.a.c2.x * x.a.c2.x + x.a.c2.y * minusX5
  shared + (x.b.c0.x * x.b.c0.y + x.a.c2.x * x.a.c2.y) * 2 +
    x.a.c1.y * cyclotomicNegTwoThird

def cppCyclotomicP4 (x : SelectedFp12) : SelectedBaseField :=
  let minusX3 := -x.a.c1.y
  let minusX11 := -x.b.c2.y
  let shared := x.b.c2.x * x.b.c2.x + x.b.c2.y * minusX11
  shared + (x.b.c2.x * minusX11) * 2 +
    x.a.c1.x * x.a.c1.x + x.a.c1.y * minusX3 +
    x.a.c2.x * cyclotomicNegTwoThird

def cppCyclotomicP5 (x : SelectedFp12) : SelectedBaseField :=
  let minusX11 := -x.b.c2.y
  let shared := x.b.c2.x * x.b.c2.x + x.b.c2.y * minusX11
  shared + (x.a.c1.x * x.a.c1.y + x.b.c2.x * x.b.c2.y) * 2 +
    x.a.c2.y * cyclotomicNegTwoThird

def cppCyclotomicP6 (x : SelectedFp12) : SelectedBaseField :=
  let minusX11 := -x.b.c2.y
  let positive := x.a.c1.x * x.b.c2.x + x.a.c1.y * minusX11
  let negative := x.a.c1.x * x.b.c2.y + x.a.c1.y * x.b.c2.x
  positive - negative + x.b.c0.x * cyclotomicOneThird

def cppCyclotomicP7 (x : SelectedFp12) : SelectedBaseField :=
  let minusX11 := -x.b.c2.y
  let positive := x.a.c1.x * x.b.c2.x + x.a.c1.y * minusX11
  let negative := x.a.c1.x * x.b.c2.y + x.a.c1.y * x.b.c2.x
  positive + negative + x.b.c0.y * cyclotomicOneThird

def cppCyclotomicP8 (x : SelectedFp12) : SelectedBaseField :=
  let minusX9 := -x.b.c1.y
  x.a.c0.x * x.b.c1.x + x.a.c0.y * minusX9 +
    x.b.c1.x * cyclotomicOneThird

def cppCyclotomicP9 (x : SelectedFp12) : SelectedBaseField :=
  x.a.c0.x * x.b.c1.y + x.a.c0.y * x.b.c1.x +
    x.b.c1.y * cyclotomicOneThird

def cppCyclotomicP10 (x : SelectedFp12) : SelectedBaseField :=
  let minusX7 := -x.b.c0.y
  x.a.c2.x * x.b.c0.x + x.a.c2.y * minusX7 +
    x.b.c2.x * cyclotomicOneThird

def cppCyclotomicP11 (x : SelectedFp12) : SelectedBaseField :=
  x.a.c2.x * x.b.c0.y + x.a.c2.y * x.b.c0.x +
    x.b.c2.y * cyclotomicOneThird

def cppFp12CyclotomicSquareDirect (x : SelectedFp12) : SelectedFp12 :=
  ⟨⟨⟨3 * cppCyclotomicP0 x, 3 * cppCyclotomicP1 x⟩,
      ⟨3 * cppCyclotomicP2 x, 3 * cppCyclotomicP3 x⟩,
      ⟨3 * cppCyclotomicP4 x, 3 * cppCyclotomicP5 x⟩⟩,
    ⟨⟨6 * cppCyclotomicP6 x, 6 * cppCyclotomicP7 x⟩,
      ⟨6 * cppCyclotomicP8 x, 6 * cppCyclotomicP9 x⟩,
      ⟨6 * cppCyclotomicP10 x, 6 * cppCyclotomicP11 x⟩⟩⟩

theorem cpp_cyclotomic_square_eq_direct (x : SelectedFp12) :
    cppFp12CyclotomicSquare x = cppFp12CyclotomicSquareDirect x := by
  ext <;> rfl

def cppFp12Frobenius2Direct (x : SelectedFp12) : SelectedFp12 :=
  ⟨⟨⟨x.a.c0.x, x.a.c0.y⟩,
      ⟨x.a.c1.x * frob1_4, x.a.c1.y * frob1_4⟩,
      ⟨x.a.c2.x * frob1_2, x.a.c2.y * frob1_2⟩⟩,
    ⟨⟨x.b.c0.x * frob2_3, x.b.c0.y * frob2_3⟩,
      ⟨-x.b.c1.x, -x.b.c1.y⟩,
      ⟨x.b.c2.x * frob1_3, x.b.c2.y * frob1_3⟩⟩⟩

theorem cpp_frobenius2_eq_direct (x : SelectedFp12) :
    cppFp12Frobenius2 x = cppFp12Frobenius2Direct x := by
  ext <;> rfl

@[simp] theorem unpack_frobenius2_direct
    (x : SelectedFp12) (index : Fin 12) :
    unpackFlatFp12 (cppFp12Frobenius2Direct x) index =
      cppFp12Frobenius2Flat (unpackFlatFp12 x) index := by
  fin_cases index <;> rfl

theorem frobenius2_coefficient_square_facts :
    frob1_4 ^ 2 = frob1_2 ∧
      frob1_2 ^ 2 = frob1_4 ∧
      frob2_3 ^ 2 = frob1_4 ∧
      frob1_3 ^ 2 = frob1_2 := by
  native_decide

theorem frobenius2_coefficient_linear_facts :
    frob1_4 = -1 - frob1_2 ∧
      frob2_3 = -frob1_2 ∧
      frob1_3 = 1 + frob1_2 := by
  native_decide

theorem mul_frob1_4_twice (value : SelectedBaseField) :
    (value * frob1_4) * frob1_4 = value * frob1_2 := by
  rw [mul_assoc, ← pow_two, frobenius2_coefficient_square_facts.1]

theorem mul_frob1_2_twice (value : SelectedBaseField) :
    (value * frob1_2) * frob1_2 = value * frob1_4 := by
  rw [mul_assoc, ← pow_two, frobenius2_coefficient_square_facts.2.1]

theorem mul_frob2_3_twice (value : SelectedBaseField) :
    (value * frob2_3) * frob2_3 = value * frob1_4 := by
  rw [mul_assoc, ← pow_two,
    frobenius2_coefficient_square_facts.2.2.1]

theorem mul_frob1_3_twice (value : SelectedBaseField) :
    (value * frob1_3) * frob1_3 = value * frob1_2 := by
  rw [mul_assoc, ← pow_two,
    frobenius2_coefficient_square_facts.2.2.2]

def cppFp12Frobenius4Direct (x : SelectedFp12) : SelectedFp12 :=
  ⟨⟨⟨x.a.c0.x, x.a.c0.y⟩,
      ⟨x.a.c1.x * frob1_2, x.a.c1.y * frob1_2⟩,
      ⟨x.a.c2.x * frob1_4, x.a.c2.y * frob1_4⟩⟩,
    ⟨⟨x.b.c0.x * frob1_4, x.b.c0.y * frob1_4⟩,
      ⟨x.b.c1.x, x.b.c1.y⟩,
      ⟨x.b.c2.x * frob1_2, x.b.c2.y * frob1_2⟩⟩⟩

def cppFp12Frobenius4Flat (x : FlatFp12 SelectedBaseField) :
    FlatFp12 SelectedBaseField :=
  ![x 0, x 1,
    x 2 * frob1_2, x 3 * frob1_2,
    x 4 * frob1_4, x 5 * frob1_4,
    x 6 * frob1_4, x 7 * frob1_4,
    x 8, x 9,
    x 10 * frob1_2, x 11 * frob1_2]

@[simp] theorem unpack_frobenius4_direct
    (x : SelectedFp12) (index : Fin 12) :
    unpackFlatFp12 (cppFp12Frobenius4Direct x) index =
      cppFp12Frobenius4Flat (unpackFlatFp12 x) index := by
  fin_cases index <;> rfl

@[simp] theorem frobenius2_direct_a_c0 (x : SelectedFp12) :
    (cppFp12Frobenius2Direct x).a.c0 =
      ⟨x.a.c0.x, x.a.c0.y⟩ := rfl
@[simp] theorem frobenius2_direct_a_c1 (x : SelectedFp12) :
    (cppFp12Frobenius2Direct x).a.c1 =
      ⟨x.a.c1.x * frob1_4, x.a.c1.y * frob1_4⟩ := rfl
@[simp] theorem frobenius2_direct_a_c2 (x : SelectedFp12) :
    (cppFp12Frobenius2Direct x).a.c2 =
      ⟨x.a.c2.x * frob1_2, x.a.c2.y * frob1_2⟩ := rfl
@[simp] theorem frobenius2_direct_b_c0 (x : SelectedFp12) :
    (cppFp12Frobenius2Direct x).b.c0 =
      ⟨x.b.c0.x * frob2_3, x.b.c0.y * frob2_3⟩ := rfl
@[simp] theorem frobenius2_direct_b_c1 (x : SelectedFp12) :
    (cppFp12Frobenius2Direct x).b.c1 =
      ⟨-x.b.c1.x, -x.b.c1.y⟩ := rfl
@[simp] theorem frobenius2_direct_b_c2 (x : SelectedFp12) :
    (cppFp12Frobenius2Direct x).b.c2 =
      ⟨x.b.c2.x * frob1_3, x.b.c2.y * frob1_3⟩ := rfl

@[simp] theorem frobenius4_direct_a_c0 (x : SelectedFp12) :
    (cppFp12Frobenius4Direct x).a.c0 =
      ⟨x.a.c0.x, x.a.c0.y⟩ := rfl
@[simp] theorem frobenius4_direct_a_c1 (x : SelectedFp12) :
    (cppFp12Frobenius4Direct x).a.c1 =
      ⟨x.a.c1.x * frob1_2, x.a.c1.y * frob1_2⟩ := rfl
@[simp] theorem frobenius4_direct_a_c2 (x : SelectedFp12) :
    (cppFp12Frobenius4Direct x).a.c2 =
      ⟨x.a.c2.x * frob1_4, x.a.c2.y * frob1_4⟩ := rfl
@[simp] theorem frobenius4_direct_b_c0 (x : SelectedFp12) :
    (cppFp12Frobenius4Direct x).b.c0 =
      ⟨x.b.c0.x * frob1_4, x.b.c0.y * frob1_4⟩ := rfl
@[simp] theorem frobenius4_direct_b_c1 (x : SelectedFp12) :
    (cppFp12Frobenius4Direct x).b.c1 =
      ⟨x.b.c1.x, x.b.c1.y⟩ := rfl
@[simp] theorem frobenius4_direct_b_c2 (x : SelectedFp12) :
    (cppFp12Frobenius4Direct x).b.c2 =
      ⟨x.b.c2.x * frob1_2, x.b.c2.y * frob1_2⟩ := rfl

@[simp] theorem frobenius4_direct_a (x : SelectedFp12) :
    (cppFp12Frobenius4Direct x).a =
      ⟨⟨x.a.c0.x, x.a.c0.y⟩,
        ⟨x.a.c1.x * frob1_2, x.a.c1.y * frob1_2⟩,
        ⟨x.a.c2.x * frob1_4, x.a.c2.y * frob1_4⟩⟩ := rfl

@[simp] theorem frobenius4_direct_b (x : SelectedFp12) :
    (cppFp12Frobenius4Direct x).b =
      ⟨⟨x.b.c0.x * frob1_4, x.b.c0.y * frob1_4⟩,
        ⟨x.b.c1.x, x.b.c1.y⟩,
        ⟨x.b.c2.x * frob1_2, x.b.c2.y * frob1_2⟩⟩ := rfl

end CatBLST.AVX.MatrixNoK
