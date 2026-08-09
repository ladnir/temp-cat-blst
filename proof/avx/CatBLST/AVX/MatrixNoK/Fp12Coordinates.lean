/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Fp12Frobenius
import Mathlib.Tactic

/-! Coordinate lemmas shared by the selected Fp12 source-schedule proofs. -/

namespace CatBLST.AVX.MatrixNoK

open Fp12Concrete

theorem unpack_flat_fp12_injective :
    Function.Injective (unpackFlatFp12 : SelectedFp12 → FlatFp12 SelectedBaseField) := by
  intro x y equal
  rw [← pack_unpack_fp12 x, ← pack_unpack_fp12 y, equal]

@[simp] theorem unpack_flat_fp12_zero (index : Fin 12) :
    unpackFlatFp12 (0 : SelectedFp12) index = 0 := by
  fin_cases index <;> rfl

@[simp] theorem unpack_flat_fp12_add
    (x y : SelectedFp12) (index : Fin 12) :
    unpackFlatFp12 (x + y) index =
      unpackFlatFp12 x index + unpackFlatFp12 y index := by
  fin_cases index <;> rfl

@[simp] theorem unpack_flat_fp12_zero_fun :
    unpackFlatFp12 (0 : SelectedFp12) = 0 := by
  funext index
  exact unpack_flat_fp12_zero index

@[simp] theorem unpack_flat_fp12_add_fun (x y : SelectedFp12) :
    unpackFlatFp12 (x + y) = unpackFlatFp12 x + unpackFlatFp12 y := by
  funext index
  exact unpack_flat_fp12_add x y index

@[simp] theorem unpack_flat_fp12_scalar_mul
    (coefficient : SelectedBaseField) (x : SelectedFp12) (index : Fin 12) :
    unpackFlatFp12 (scalar coefficient * x) index =
      coefficient * unpackFlatFp12 x index := by
  fin_cases index <;>
    simp [unpackFlatFp12, scalar, Fp12Value.mul, Fp6Value.add,
      Fp6Value.mul, Fp6Value.mulByV, Fp2Value.add, Fp2Value.mul,
      Fp2Value.mulByNonresidue]

theorem unpack_flat_fp12_scalar_mul_fun
    (coefficient : SelectedBaseField) (x : SelectedFp12) :
    unpackFlatFp12 (scalar coefficient * x) =
      fun index => coefficient * unpackFlatFp12 x index := by
  funext index
  exact unpack_flat_fp12_scalar_mul coefficient x index

end CatBLST.AVX.MatrixNoK
