/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Fp12Coordinates
import CatBLST.AVX.MatrixNoK.MillerDouble
import Mathlib.Tactic

set_option maxHeartbeats 800000

/-!
# Literal sparse Fp12 Miller-line multiplication

This file transcribes `FP12Ring::flat_mulxy00z0`.  The source multiplies a
flat Fp12 accumulator by the sparse line element whose nonzero Fp2
coordinates are `c0`, `c1`, and `c4`.  The optimized polynomial schedule is
proved equal to the already accepted dense Fp12 multiplication; it is not
replaced by that dense implementation in C++.
-/

namespace CatBLST.AVX.MatrixNoK

/-- Flat Fp12 layout of the Miller line `(c, x, 0, 0, y, 0)`. -/
def sparseMillerLineFlat {R : Type*} [Zero R]
    (constantCoord xCoord yCoord : Fp2Value R) : FlatFp12 R :=
  ![constantCoord.x, constantCoord.y,
    xCoord.x, xCoord.y,
    0, 0,
    0, 0,
    yCoord.x, yCoord.y,
    0, 0]

/--
Literal grouping and operation order of C++ `flat_mulxy00z0`.  `offset`
models the positive RNS offset inserted before source-level subtractions; its
canonical field value is zero.
-/
def cppFp12LiteralMulXY00Z0 {R : Type*} [CommRing R]
    (offset : R) (x : FlatFp12 R)
    (constantCoord xCoord yCoord : Fp2Value R) : FlatFp12 R :=
  let y0 := constantCoord.x
  let y1 := constantCoord.y
  let y2 := xCoord.x
  let y3 := xCoord.y
  let y8 := yCoord.x
  let y9 := yCoord.y
  let minusY1 := -y1
  let minusY3 := -y3
  let minusY9 := -y9
  let pos0 := offset + x 4 * y2 + x 5 * minusY3 +
    x 8 * y8 + x 9 * minusY9
  let neg0 := x 4 * y3 + x 5 * y2 + x 8 * y9 + x 9 * y8
  let pos2 := offset + x 10 * y8 + x 11 * minusY9
  let neg2 := x 10 * y9 + x 11 * y8
  let pos6 := offset + x 4 * y8 + x 5 * minusY9 +
    x 10 * y2 + x 11 * minusY3
  let neg6 := x 4 * y9 + x 5 * y8 + x 10 * y3 + x 11 * y2
  ![
    pos0 - neg0 + x 0 * y0 + x 1 * minusY1,
    pos0 + neg0 + x 0 * y1 + x 1 * y0,
    pos2 - neg2 + x 0 * y2 + x 1 * minusY3 +
      x 2 * y0 + x 3 * minusY1,
    pos2 + neg2 + x 0 * y3 + x 1 * y2 + x 2 * y1 + x 3 * y0,
    x 2 * y2 + x 3 * minusY3 + x 4 * y0 + x 5 * minusY1 +
      x 6 * y8 + x 7 * minusY9,
    x 2 * y3 + x 3 * y2 + x 4 * y1 + x 5 * y0 +
      x 6 * y9 + x 7 * y8,
    pos6 - neg6 + x 6 * y0 + x 7 * minusY1,
    pos6 + neg6 + x 6 * y1 + x 7 * y0,
    x 0 * y8 + x 1 * minusY9 + x 6 * y2 + x 7 * minusY3 +
      x 8 * y0 + x 9 * minusY1,
    x 0 * y9 + x 1 * y8 + x 6 * y3 + x 7 * y2 +
      x 8 * y1 + x 9 * y0,
    x 2 * y8 + x 3 * minusY9 + x 8 * y2 + x 9 * minusY3 +
      x 10 * y0 + x 11 * minusY1,
    x 2 * y9 + x 3 * y8 + x 8 * y3 + x 9 * y2 +
      x 10 * y1 + x 11 * y0]

theorem cpp_fp12_literal_mulxy00z0_matches_dense
    {R : Type*} [CommRing R] (x : FlatFp12 R)
    (constantCoord xCoord yCoord : Fp2Value R) :
    cppFp12LiteralMulXY00Z0 0 x constantCoord xCoord yCoord =
      cppFp12LiteralMul 0 x
        (sparseMillerLineFlat constantCoord xCoord yCoord) := by
  funext index
  fin_cases index <;>
    simp [cppFp12LiteralMulXY00Z0, cppFp12LiteralMul,
      sparseMillerLineFlat]

theorem cpp_fp12_literal_mulxy00z0_packs
    {R : Type*} [CommRing R] (x : FlatFp12 R)
    (constantCoord xCoord yCoord : Fp2Value R) :
    packFlatFp12
        (cppFp12LiteralMulXY00Z0 0 x constantCoord xCoord yCoord) =
      packFlatFp12 x *
        packFlatFp12 (sparseMillerLineFlat constantCoord xCoord yCoord) := by
  rw [cpp_fp12_literal_mulxy00z0_matches_dense,
    cpp_fp12_literal_mul_packs, cpp_fp12_flat_mul_correct]
  rfl

/-- The accepted Fp12 meaning of a source-shaped Miller line. -/
def MillerLine.toSparseFp12 (line : MillerLine SelectedMillerFp2) :
    SelectedFp12 :=
  packFlatFp12 (sparseMillerLineFlat line.c line.x line.y)

theorem selected_cpp_flat_mulxy00z0_correct
    (accumulator : SelectedFp12) (line : MillerLine SelectedMillerFp2) :
    packFlatFp12
        (cppFp12LiteralMulXY00Z0 0 (unpackFlatFp12 accumulator)
          line.c line.x line.y) =
      accumulator * line.toSparseFp12 := by
  rw [cpp_fp12_literal_mulxy00z0_packs, pack_unpack_fp12]
  rfl

/-!
## Physical sparse-kernel envelope

The selected MatrixNoK standard representation has element bound two and the
Miller/Fp12 kernels allow RNS bounds through six.  Consequently one source
product contributes at most `2 * 2 * 6 * 6 = 144` RNS units.  The three
subtractive coordinate pairs use the same 1440-unit offset accepted in P3.
-/

def sparseMulUnsignedAdditionUnits : FlatFp12 Nat :=
  ![6, 10, 6, 8, 6, 6, 6, 10, 6, 6, 6, 6]

def sparseMulUnsignedSubtractionUnits : FlatFp12 Nat :=
  ![4, 0, 2, 0, 0, 0, 4, 0, 0, 0, 0, 0]

def sparseMulOffsetMultiplicity : FlatFp12 Nat :=
  ![1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0]

def sparseMulRnsLower (index : Fin 12) : Nat :=
  sparseMulOffsetMultiplicity index * fp12Version3OffsetRnsBound -
    sparseMulUnsignedSubtractionUnits index * fp12ProductRnsUnit

def sparseMulRnsUpper (index : Fin 12) : Nat :=
  sparseMulOffsetMultiplicity index * fp12Version3OffsetRnsBound +
    sparseMulUnsignedAdditionUnits index * fp12ProductRnsUnit

theorem sparse_mul_unit_bounds (index : Fin 12) :
    sparseMulUnsignedAdditionUnits index <= 10 /\
      sparseMulUnsignedSubtractionUnits index <= 4 /\
      sparseMulOffsetMultiplicity index <= 1 := by
  fin_cases index <;>
    norm_num [sparseMulUnsignedAdditionUnits,
      sparseMulUnsignedSubtractionUnits, sparseMulOffsetMultiplicity]

theorem sparse_mul_rns_bounds (index : Fin 12) :
    sparseMulUnsignedSubtractionUnits index * fp12ProductRnsUnit <=
        sparseMulOffsetMultiplicity index * fp12Version3OffsetRnsBound /\
      sparseMulRnsUpper index <= 2880 := by
  fin_cases index <;>
    norm_num [sparseMulUnsignedSubtractionUnits,
      sparseMulUnsignedAdditionUnits, sparseMulOffsetMultiplicity,
      sparseMulRnsUpper, fp12ProductRnsUnit, fp12MaximumInputRnsBound,
      fp12Version3OffsetRnsBound, standardBound]

theorem sparse_mul_element_bounds (index : Fin 12) :
    sparseMulUnsignedSubtractionUnits index * 4 <=
        sparseMulOffsetMultiplicity index * 80 /\
      sparseMulOffsetMultiplicity index * 81 +
          sparseMulUnsignedAdditionUnits index * 4 <= 121 := by
  fin_cases index <;>
    norm_num [sparseMulUnsignedSubtractionUnits,
      sparseMulUnsignedAdditionUnits, sparseMulOffsetMultiplicity]

/-- Source operation ledger after accounting for the three shared pairs. -/
structure SparseMulCallLedger where
  baseProducts : Nat
  sourceNegations : Nat
  wideOffsets : Nat
  reducedCoordinates : Nat
  deriving DecidableEq

def cppFlatMulXY00Z0Ledger : SparseMulCallLedger :=
  { baseProducts := 72
    sourceNegations := 3
    wideOffsets := 3
    reducedCoordinates := 12 }

theorem cpp_flat_mulxy00z0_ledger_exact :
    cppFlatMulXY00Z0Ledger =
      { baseProducts := 72
        sourceNegations := 3
        wideOffsets := 3
        reducedCoordinates := 12 } := rfl

end CatBLST.AVX.MatrixNoK
