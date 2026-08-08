/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Generated
import Mathlib.Tactic

/-!
# MatrixNoK fixed-parameter invariants

The generated module binds the exact C++ literals. This file proves the
mathematical properties consumed by the BLS12-381 MatrixNoK construction.
Primality is deliberately not required: pairwise coprimality and oddness are
the actual RNS preconditions.
-/

namespace CatBLST.AVX.MatrixNoK

def allModuli : List Nat := Generated.moduli1 ++ Generated.moduli2

theorem dimensions :
    Generated.BLS12_381_RNS_DIM = 8 ∧
    Generated.BLS12_381_PADDED_DIM = 8 ∧
    Generated.BLS12_381_CONVERT_TO_DIGITS = 8 ∧
    Generated.BLS12_381_CONVERT_TO_RADIX_BITS = 50 := by
  decide

theorem basis_shapes :
    Generated.moduli1.length = 8 ∧
    Generated.moduli2.length = 8 ∧
    Generated.moduli1_padded = Generated.moduli1 ∧
    Generated.moduli2_padded = Generated.moduli2 := by
  decide

theorem all_moduli_distinct : allModuli.Nodup := by
  decide

theorem all_moduli_odd : ∀ modulus ∈ allModuli, Odd modulus := by
  norm_num [allModuli, Generated.moduli1, Generated.moduli2]

theorem all_moduli_lt_radix :
    ∀ modulus ∈ allModuli,
      modulus < 2 ^ Generated.BLS12_381_CONVERT_TO_RADIX_BITS := by
  norm_num [allModuli, Generated.moduli1, Generated.moduli2,
    Generated.BLS12_381_CONVERT_TO_RADIX_BITS]

theorem all_moduli_pairwise_coprime : allModuli.Pairwise Nat.Coprime := by
  norm_num [allModuli, Generated.moduli1, Generated.moduli2, Nat.Coprime]

theorem M1_eq_product : Generated.M1 = Generated.moduli1.prod := by
  norm_num [Generated.M1, Generated.moduli1]

theorem M2_eq_product : Generated.M2 = Generated.moduli2.prod := by
  norm_num [Generated.M2, Generated.moduli2]

theorem combined_modulus_eq_product :
    Generated.M1M2 = Generated.M1 * Generated.M2 := by
  norm_num [Generated.M1M2, Generated.M1, Generated.M2]

theorem basis_dynamic_range :
    6 * Generated.target_modulus < Generated.M1 ∧
      6 * Generated.target_modulus < Generated.M2 := by
  norm_num [Generated.target_modulus, Generated.M1, Generated.M2]

theorem total_rotation_inverse :
    Generated.total_rotation * Generated.total_inv_rotation % Generated.M1M2 = 1 := by
  norm_num [Generated.total_rotation, Generated.total_inv_rotation, Generated.M1M2]

theorem target_inverse_factor :
    Generated.M1 * Generated.inv_factor_target % Generated.target_modulus = 1 := by
  norm_num [Generated.M1, Generated.inv_factor_target, Generated.target_modulus]

def montgomeryRadix : Nat := 2 ^ 52

/-- The M1 inverse modulo M2 encoded inside the convert-to postmultiplier. -/
def crossBasisInverse : Nat :=
  Generated.convert_to_postmult / (montgomeryRadix * montgomeryRadix)

theorem cross_basis_inverse_relation :
    Generated.M1 * crossBasisInverse % Generated.M2 = 1 := by
  norm_num [Generated.M1, Generated.M2, Generated.convert_to_postmult,
    crossBasisInverse, montgomeryRadix]

theorem r1_premult_relation :
    (Generated.qr_premult * Generated.qr_premult * Generated.r1_premult *
      (Generated.M1 - Generated.target_modulus) * montgomeryRadix) %
        Generated.M1 = 1 := by
  norm_num [Generated.qr_premult, Generated.r1_premult, Generated.M1,
    Generated.target_modulus, montgomeryRadix]

theorem r1_postmult_relation :
    Generated.r1_postmult =
      Generated.target_modulus * crossBasisInverse * crossBasisInverse *
        montgomeryRadix * montgomeryRadix := by
  norm_num [Generated.r1_postmult, Generated.target_modulus,
    Generated.convert_to_postmult, crossBasisInverse, montgomeryRadix]

theorem r2_premult_relation :
    Generated.r2_premult * montgomeryRadix % Generated.M2 =
      Generated.M1 % Generated.M2 := by
  norm_num [Generated.r2_premult, Generated.M1, Generated.M2,
    montgomeryRadix]

theorem r2_postmult_relation :
    Generated.r2_postmult =
      montgomeryRadix * montgomeryRadix * Generated.qr_premult := by
  norm_num [Generated.r2_postmult, Generated.qr_premult, montgomeryRadix]

theorem conversion_multiplier_relations :
    Generated.convert_to_premult = Generated.M1 ∧
    Generated.convert_to_postmult =
      crossBasisInverse * montgomeryRadix * montgomeryRadix ∧
    Generated.convert_from_premult = Generated.r2_premult ∧
    Generated.convert_from_postmult = Generated.inv_factor_target := by
  norm_num [Generated.convert_to_premult, Generated.M1,
    Generated.convert_to_postmult, crossBasisInverse, montgomeryRadix,
    Generated.convert_from_premult, Generated.r2_premult,
    Generated.convert_from_postmult, Generated.inv_factor_target]

theorem reducer1_montgomery_factors :
    List.Forall₂ (fun modulus factor =>
      modulus * factor % montgomeryRadix = 1)
      Generated.moduli1 Generated.reducer1_mont_factor := by
  norm_num [Generated.moduli1, Generated.reducer1_mont_factor, montgomeryRadix]

theorem reducer2_montgomery_factors :
    List.Forall₂ (fun modulus factor =>
      modulus * factor % montgomeryRadix = 1)
      Generated.moduli2 Generated.reducer2_mont_factor := by
  norm_num [Generated.moduli2, Generated.reducer2_mont_factor, montgomeryRadix]

def isMatrix8 (matrix : List (List Nat)) : Prop :=
  matrix.length = 8 ∧ ∀ row ∈ matrix, row.length = 8

def entriesBelow (bound : Nat) (matrix : List (List Nat)) : Prop :=
  ∀ row ∈ matrix, ∀ value ∈ row, value < bound

def allFixedMatrices : List (List (List Nat)) := [
  Generated.r1_perm_mat_fixed,
  Generated.r1_rns_mat,
  Generated.r2_rns_mat,
  Generated.convert_to_rns_mat,
  Generated.convert_from_rns_mat
]

def allCorrections : List (List Nat) := [
  Generated.r1_correction,
  Generated.r2_correction,
  Generated.convert_to_correction,
  Generated.convert_from_correction
]

def allQuotientEstimates : List (List Nat) := [
  Generated.r1_shifted_quotient_estimations,
  Generated.r2_shifted_quotient_estimations,
  Generated.convert_to_shifted_quotient_estimations,
  Generated.convert_from_shifted_quotient_estimations
]

theorem fixed_matrix_shapes : ∀ matrix ∈ allFixedMatrices, isMatrix8 matrix := by
  unfold allFixedMatrices isMatrix8
  decide

theorem fixed_matrix_entries_fit_digit :
    ∀ matrix ∈ allFixedMatrices, entriesBelow montgomeryRadix matrix := by
  unfold allFixedMatrices entriesBelow montgomeryRadix
  decide

theorem correction_shapes_and_ranges :
    (∀ correction ∈ allCorrections, correction.length = 8) ∧
    (∀ correction ∈ allCorrections, ∀ value ∈ correction,
      value < montgomeryRadix) := by
  decide

theorem quotient_estimate_shapes_and_ranges :
    (∀ estimates ∈ allQuotientEstimates, estimates.length = 8) ∧
    (∀ estimates ∈ allQuotientEstimates, ∀ value ∈ estimates,
      value < 2 ^ 64) := by
  decide

theorem reducer_factor_shapes_and_ranges :
    Generated.reducer1_mont_factor.length = 8 ∧
    Generated.reducer2_mont_factor.length = 8 ∧
    (∀ value ∈ Generated.reducer1_mont_factor, value < montgomeryRadix) ∧
    (∀ value ∈ Generated.reducer2_mont_factor, value < montgomeryRadix) := by
  decide

def matrixEntry (matrix : List (List Nat)) (row column : Nat) : Nat :=
  (matrix.getD row []).getD column 0

/-- Source lane selected by cyclic diagonal `diagonal` at output `column`. -/
def r1SourceRow (diagonal column : Fin 8) : Fin 8 := column - diagonal

theorem r1SourceRow_bijective (diagonal : Fin 8) :
    Function.Bijective (r1SourceRow diagonal) := by
  constructor
  · intro left right equal
    have shifted := congrArg (fun value => value + diagonal) equal
    simpa [r1SourceRow] using shifted
  · intro value
    refine ⟨value + diagonal, ?_⟩
    simp [r1SourceRow]

theorem r1_permutation_shapes :
    Generated.r1_perm_mat_fixed.length = 8 ∧
    (∀ row ∈ Generated.r1_perm_mat_fixed, row.length = 8) ∧
    Generated.r1_rns_mat.length = 8 ∧
    (∀ row ∈ Generated.r1_rns_mat, row.length = 8) := by
  decide

theorem r1_permutation_matches_base :
    ∀ diagonal column : Fin 8,
      matrixEntry Generated.r1_perm_mat_fixed diagonal.val column.val =
        matrixEntry Generated.r1_rns_mat
          (r1SourceRow diagonal column).val column.val := by
  decide

theorem r1_correction_is_skipped : Generated.r1_skip_correction = true := by
  decide

end CatBLST.AVX.MatrixNoK
