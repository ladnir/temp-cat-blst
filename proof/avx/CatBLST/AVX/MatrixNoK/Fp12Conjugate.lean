/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Fp12Frobenius3
import Mathlib.Tactic

set_option maxHeartbeats 800000

/-!
# Selected MatrixNoK Fp12 conjugation

The source keeps the first six flat coefficients and applies
`standard_negate` to the last six.  This is a general Fp12 operation: it is
Frobenius to the sixth power.  It becomes multiplicative inversion only after
the operand is known to be cyclotomic and nonzero.
-/

namespace CatBLST.AVX.MatrixNoK

def cppFp12ConjugateFlat (x : FlatFp12 SelectedBaseField) :
    FlatFp12 SelectedBaseField :=
  ![x 0, x 1, x 2, x 3, x 4, x 5,
    -x 6, -x 7, -x 8, -x 9, -x 10, -x 11]

def cppFp12Conjugate (x : SelectedFp12) : SelectedFp12 :=
  packFlatFp12 (cppFp12ConjugateFlat (unpackFlatFp12 x))

theorem cpp_frobenius3_twice_basis_eq_conjugate :
    forall index : Fin 12,
      cppFp12Frobenius3 (cppFp12Frobenius3 (fp12Basis index)) =
        cppFp12Conjugate (fp12Basis index) := by
  native_decide

def selectedFrobenius3 (x : SelectedFp12) : SelectedFp12 :=
  frobenius SelectedFp12 Generated.target_modulus
    (frobenius SelectedFp12 Generated.target_modulus
      (frobenius SelectedFp12 Generated.target_modulus x))

def selectedFrobenius6 (x : SelectedFp12) : SelectedFp12 :=
  selectedFrobenius3 (selectedFrobenius3 x)

theorem selected_frobenius3_eq_cpp (x : SelectedFp12) :
    selectedFrobenius3 x = cppFp12Frobenius3 x :=
  (cpp_selected_fp12_frobenius3_iterate x).symm

theorem selected_frobenius6_mul (x y : SelectedFp12) :
    selectedFrobenius6 (x * y) =
      selectedFrobenius6 x * selectedFrobenius6 y := by
  simp only [selectedFrobenius6, selectedFrobenius3, map_mul]

theorem selected_frobenius6_sum (f : Fin 12 -> SelectedFp12) :
    selectedFrobenius6 (Finset.univ.sum f) =
      Finset.univ.sum (fun i => selectedFrobenius6 (f i)) := by
  simp only [selectedFrobenius6, selectedFrobenius3, map_sum]

theorem selected_frobenius6_scalar (value : SelectedBaseField) :
    selectedFrobenius6 (Fp12Concrete.scalar value) =
      Fp12Concrete.scalar value := by
  simp [selectedFrobenius6, selectedFrobenius3,
    selected_frobenius1_scalar]

theorem selected_frobenius6_basis (index : Fin 12) :
    selectedFrobenius6 (fp12Basis index) =
      cppFp12Conjugate (fp12Basis index) := by
  rw [selectedFrobenius6, selected_frobenius3_eq_cpp,
    selected_frobenius3_eq_cpp]
  exact cpp_frobenius3_twice_basis_eq_conjugate index

theorem cpp_conjugate_basis_expansion (x : SelectedFp12) :
    cppFp12Conjugate x =
      Finset.univ.sum (fun index : Fin 12 =>
        Fp12Concrete.scalar (unpackFlatFp12 x index) *
          cppFp12Conjugate (fp12Basis index)) := by
  ext <;>
    simp [cppFp12Conjugate, cppFp12ConjugateFlat,
      fp12Basis, unpackFlatFp12, packFlatFp12, Fin.sum_univ_succ,
      Fp12Concrete.scalar, Fp12Concrete.add, Fp12Value.mul,
      Fp6Value.add, Fp6Value.mul, Fp6Value.mulByV,
      Fp2Value.add, Fp2Value.mul, Fp2Value.mulByNonresidue]

theorem cpp_fp12_conjugate_correct (x : SelectedFp12) :
    cppFp12Conjugate x = x ^ (Generated.target_modulus ^ 6) := by
  rw [cpp_conjugate_basis_expansion]
  have mappedExpansion :
      selectedFrobenius6 (fp12BasisExpansion x) =
        Finset.univ.sum (fun index : Fin 12 =>
          Fp12Concrete.scalar (unpackFlatFp12 x index) *
            cppFp12Conjugate (fp12Basis index)) := by
    rw [fp12BasisExpansion, selected_frobenius6_sum]
    apply Finset.sum_congr rfl
    intro index _
    rw [selected_frobenius6_mul, selected_frobenius6_scalar,
      selected_frobenius6_basis]
  rw [fp12_basis_expansion] at mappedExpansion
  rw [mappedExpansion.symm]
  simp [selectedFrobenius6, selectedFrobenius3, frobenius_def,
    pow_mul, pow_succ]

end CatBLST.AVX.MatrixNoK
