/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Fp12
import Mathlib.Tactic

/-!
# Literal selected Fp12 scalar schedules

This file transcribes the twelve scalar polynomials in C++ `fp12_flat_mul`
and `fp12_flat_sqr`. It is deliberately separate from the tower-level proof:
the theorems below mechanically establish that pairing adjacent scalar outputs
reconstructs the six Fp2 coefficients already proved in `Fp12.lean`.
-/

namespace CatBLST.AVX.MatrixNoK

abbrev FlatFp12 (R : Type*) := Fin 12 → R

def packFlatFp12 {R : Type*} (x : FlatFp12 R) : Fp12Value R :=
  ⟨⟨⟨x 0, x 1⟩, ⟨x 2, x 3⟩, ⟨x 4, x 5⟩⟩,
   ⟨⟨x 6, x 7⟩, ⟨x 8, x 9⟩, ⟨x 10, x 11⟩⟩⟩

/-- Literal grouping and operation order of C++ `fp12_flat_mul`. -/
def cppFp12LiteralMul {R : Type*} [CommRing R]
    (offset : R) (x y : FlatFp12 R) : FlatFp12 R :=
  let minusY1 := -y 1
  let minusY3 := -y 3
  let minusY5 := -y 5
  let minusY7 := -y 7
  let minusY9 := -y 9
  let minusY11 := -y 11
  let pos0 := offset + x 2 * y 4 + x 3 * minusY5 + x 4 * y 2 +
    x 5 * minusY3 + x 6 * y 10 + x 7 * minusY11 + x 8 * y 8 +
    x 9 * minusY9 + x 10 * y 6 + x 11 * minusY7
  let neg0 := x 2 * y 5 + x 3 * y 4 + x 4 * y 3 + x 5 * y 2 +
    x 6 * y 11 + x 7 * y 10 + x 8 * y 9 + x 9 * y 8 +
    x 10 * y 7 + x 11 * y 6
  let pos2 := offset + x 4 * y 4 + x 5 * minusY5 + x 8 * y 10 +
    x 9 * minusY11 + x 10 * y 8 + x 11 * minusY9
  let neg2 := x 4 * y 5 + x 5 * y 4 + x 8 * y 11 + x 9 * y 10 +
    x 10 * y 9 + x 11 * y 8
  let pos4 := offset + x 10 * y 10 + x 11 * minusY11
  let neg4 := x 10 * y 11 + x 11 * y 10
  let pos6 := offset + x 2 * y 10 + x 3 * minusY11 + x 4 * y 8 +
    x 5 * minusY9 + x 8 * y 4 + x 9 * minusY5 + x 10 * y 2 +
    x 11 * minusY3
  let neg6 := x 2 * y 11 + x 3 * y 10 + x 4 * y 9 + x 5 * y 8 +
    x 8 * y 5 + x 9 * y 4 + x 10 * y 3 + x 11 * y 2
  let pos8 := offset + x 4 * y 10 + x 5 * minusY11 + x 10 * y 4 +
    x 11 * minusY5
  let neg8 := x 4 * y 11 + x 5 * y 10 + x 10 * y 5 + x 11 * y 4
  ![
    pos0 - neg0 + x 0 * y 0 + x 1 * minusY1,
    pos0 + neg0 + x 0 * y 1 + x 1 * y 0,
    pos2 - neg2 + x 0 * y 2 + x 1 * minusY3 + x 2 * y 0 +
      x 3 * minusY1 + x 6 * y 6 + x 7 * minusY7,
    pos2 + neg2 + x 0 * y 3 + x 1 * y 2 + x 2 * y 1 + x 3 * y 0 +
      x 6 * y 7 + x 7 * y 6,
    pos4 - neg4 + x 0 * y 4 + x 1 * minusY5 + x 2 * y 2 +
      x 3 * minusY3 + x 4 * y 0 + x 5 * minusY1 + x 6 * y 8 +
      x 7 * minusY9 + x 8 * y 6 + x 9 * minusY7,
    pos4 + neg4 + x 0 * y 5 + x 1 * y 4 + x 2 * y 3 + x 3 * y 2 +
      x 4 * y 1 + x 5 * y 0 + x 6 * y 9 + x 7 * y 8 + x 8 * y 7 +
      x 9 * y 6,
    pos6 - neg6 + x 0 * y 6 + x 1 * minusY7 + x 6 * y 0 +
      x 7 * minusY1,
    pos6 + neg6 + x 0 * y 7 + x 1 * y 6 + x 6 * y 1 + x 7 * y 0,
    pos8 - neg8 + x 0 * y 8 + x 1 * minusY9 + x 2 * y 6 +
      x 3 * minusY7 + x 6 * y 2 + x 7 * minusY3 + x 8 * y 0 +
      x 9 * minusY1,
    pos8 + neg8 + x 0 * y 9 + x 1 * y 8 + x 2 * y 7 + x 3 * y 6 +
      x 6 * y 3 + x 7 * y 2 + x 8 * y 1 + x 9 * y 0,
    x 0 * y 10 + x 1 * minusY11 + x 2 * y 8 + x 3 * minusY9 +
      x 4 * y 6 + x 5 * minusY7 + x 6 * y 4 + x 7 * minusY5 +
      x 8 * y 2 + x 9 * minusY3 + x 10 * y 0 + x 11 * minusY1,
    x 0 * y 11 + x 1 * y 10 + x 2 * y 9 + x 3 * y 8 + x 4 * y 7 +
      x 5 * y 6 + x 6 * y 5 + x 7 * y 4 + x 8 * y 3 + x 9 * y 2 +
      x 10 * y 1 + x 11 * y 0]

theorem cpp_fp12_literal_mul_packs {R : Type*} [CommRing R]
    (x y : FlatFp12 R) :
    packFlatFp12 (cppFp12LiteralMul 0 x y) =
      cppFp12FlatMul (packFlatFp12 x) (packFlatFp12 y) := by
  ext <;>
    simp [packFlatFp12, cppFp12LiteralMul, cppFp12FlatMul, fp2Sum5,
      fp2Sum6, Fp2Value.add, Fp2Value.mul,
      Fp2Value.mulByNonresidue] <;>
    ring

/-- Literal grouping and operation order of C++ `fp12_flat_sqr`. -/
def cppFp12LiteralSquare {R : Type*} [CommRing R]
    (offset : R) (x : FlatFp12 R) : FlatFp12 R :=
  let minusX1 := -x 1
  let minusX3 := -x 3
  let minusX5 := -x 5
  let minusX7 := -x 7
  let minusX9 := -x 9
  let minusX11 := -x 11
  let pos0 := offset +
    (x 2 * x 4 + x 3 * minusX5 + x 6 * x 10 + x 7 * minusX11) * 2 +
    x 8 * x 8 + x 9 * minusX9
  let neg0 :=
    (x 2 * x 5 + x 3 * x 4 + x 6 * x 11 + x 7 * x 10 + x 8 * x 9) * 2
  let pos2 := offset + (x 8 * x 10 + x 9 * minusX11) * 2 +
    x 4 * x 4 + x 5 * minusX5
  let neg2 := (x 4 * x 5 + x 8 * x 11 + x 9 * x 10) * 2
  let pos4 := x 10 * x 10 + x 11 * minusX11
  let pos6 := offset + x 2 * x 10 + x 3 * minusX11 +
    x 4 * x 8 + x 5 * minusX9
  let neg6 := x 2 * x 11 + x 3 * x 10 + x 4 * x 9 + x 5 * x 8
  let pos8 := offset + x 4 * x 10 + x 5 * minusX11
  let neg8 := x 4 * x 11 + x 5 * x 10
  ![
    pos0 - neg0 + x 0 * x 0 + x 1 * minusX1,
    pos0 + neg0 + x 0 * x 1 * 2,
    pos2 - neg2 + (x 0 * x 2 + x 1 * minusX3) * 2 +
      x 6 * x 6 + x 7 * minusX7,
    pos2 + neg2 + (x 0 * x 3 + x 1 * x 2 + x 6 * x 7) * 2,
    pos4 + (x 0 * x 4 + x 1 * minusX5 + x 6 * x 8 +
      x 7 * minusX9 + x 10 * minusX11) * 2 +
      x 2 * x 2 + x 3 * minusX3,
    pos4 + (x 0 * x 5 + x 1 * x 4 + x 2 * x 3 + x 6 * x 9 +
      x 7 * x 8 + x 10 * x 11) * 2,
    (pos6 - neg6 + x 0 * x 6 + x 1 * minusX7) * 2,
    (pos6 + neg6 + x 0 * x 7 + x 1 * x 6) * 2,
    (pos8 - neg8 + x 0 * x 8 + x 1 * minusX9 + x 2 * x 6 +
      x 3 * minusX7) * 2,
    (pos8 + neg8 + x 0 * x 9 + x 1 * x 8 + x 2 * x 7 + x 3 * x 6) * 2,
    (x 0 * x 10 + x 1 * minusX11 + x 2 * x 8 + x 3 * minusX9 +
      x 4 * x 6 + x 5 * minusX7) * 2,
    (x 0 * x 11 + x 1 * x 10 + x 2 * x 9 + x 3 * x 8 +
      x 4 * x 7 + x 5 * x 6) * 2]

theorem cpp_fp12_literal_square_packs {R : Type*} [CommRing R]
    (x : FlatFp12 R) :
    packFlatFp12 (cppFp12LiteralSquare 0 x) =
      cppFp12FlatSquare (packFlatFp12 x) := by
  ext <;>
    simp [packFlatFp12, cppFp12LiteralSquare, cppFp12FlatSquare,
      fp2Sum3, Fp2Value.double, Fp2Value.add, Fp2Value.mul,
      Fp2Value.square, Fp2Value.mulByNonresidue] <;>
    ring

def fp12MulUnsignedAdditionUnits : FlatFp12 Nat :=
  ![12, 22, 12, 18, 12, 14, 12, 20, 12, 16, 12, 12]

def fp12MulUnsignedSubtractionUnits : FlatFp12 Nat :=
  ![10, 0, 6, 0, 2, 0, 8, 0, 4, 0, 0, 0]

def fp12MulOffsetMultiplicity : FlatFp12 Nat :=
  ![1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0]

def fp12SquareUnsignedAdditionUnits : FlatFp12 Nat :=
  ![12, 22, 12, 18, 14, 14, 12, 20, 12, 16, 12, 12]

def fp12SquareUnsignedSubtractionUnits : FlatFp12 Nat :=
  ![10, 0, 6, 0, 0, 0, 8, 0, 4, 0, 0, 0]

def fp12SquareOffsetMultiplicity : FlatFp12 Nat :=
  ![1, 1, 1, 1, 0, 0, 2, 2, 2, 2, 0, 0]

theorem fp12_literal_mul_unit_bounds (index : Fin 12) :
    fp12MulUnsignedAdditionUnits index ≤ 22 ∧
      fp12MulUnsignedSubtractionUnits index ≤ 10 ∧
      fp12MulOffsetMultiplicity index ≤ 1 := by
  fin_cases index <;>
    norm_num [fp12MulUnsignedAdditionUnits,
      fp12MulUnsignedSubtractionUnits, fp12MulOffsetMultiplicity]

theorem fp12_literal_square_unit_bounds (index : Fin 12) :
    fp12SquareUnsignedAdditionUnits index ≤ 22 ∧
      fp12SquareUnsignedSubtractionUnits index ≤ 10 ∧
      fp12SquareOffsetMultiplicity index ≤ 2 := by
  fin_cases index <;>
    norm_num [fp12SquareUnsignedAdditionUnits,
      fp12SquareUnsignedSubtractionUnits, fp12SquareOffsetMultiplicity]

end CatBLST.AVX.MatrixNoK
