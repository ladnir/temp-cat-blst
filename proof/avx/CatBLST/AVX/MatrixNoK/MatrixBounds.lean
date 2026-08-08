/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.BaseField
import CatBLST.AVX.MatrixNoK.Reducer
import Mathlib.Tactic

/-!
# MatrixNoK fixed-matrix accumulator bounds

These are bounds for the actual checked-in matrices, not an asymptotic CRT
estimate.  They justify the `16` correction-free R1 contribution and the `32`
R2 expansion type used by `FpRingCandidate`.
-/

namespace CatBLST.AVX.MatrixNoK

open scoped BigOperators

def matrixProductValue
    (residues : Fin 8 → Nat) (matrix : List (List Nat))
    (column : Fin 8) : Nat :=
  ∑ row, residues row * matrixEntry matrix row.val column.val

/-- One output lane of the rotate-and-accumulate R1 source schedule. -/
def r1CyclicLaneValue
    (residues : Fin 8 → Nat) (column : Fin 8) : Nat :=
  ∑ diagonal,
    residues (r1SourceRow diagonal column) *
      matrixEntry Generated.r1_perm_mat_fixed diagonal.val column.val

def r1SourceEquiv (column : Fin 8) : Fin 8 ≃ Fin 8 where
  toFun diagonal := r1SourceRow diagonal column
  invFun source := column - source
  left_inv diagonal := by simp [r1SourceRow]
  right_inv source := by simp [r1SourceRow]

/--
The fixed seven-rotation/eight-accumulation source schedule is exactly the
ordinary generated R1 matrix-vector product, for arbitrary lane values.
-/
theorem r1_cyclic_lane_eq_matrix_product
    (residues : Fin 8 → Nat) (column : Fin 8) :
    r1CyclicLaneValue residues column =
      matrixProductValue residues Generated.r1_rns_mat column := by
  unfold r1CyclicLaneValue matrixProductValue
  calc
    _ = ∑ diagonal,
        residues (r1SourceRow diagonal column) *
          matrixEntry Generated.r1_rns_mat
            (r1SourceRow diagonal column).val column.val := by
      apply Finset.sum_congr rfl
      intro diagonal _
      rw [r1_permutation_matches_base]
    _ = _ := (r1SourceEquiv column).sum_comp
      (fun source => residues source *
        matrixEntry Generated.r1_rns_mat source.val column.val)

/--
The batch C++ schedule shares the outer diagonal loop but applies the same
rotate-and-accumulate recurrence independently to every batch item.
-/
def r1BatchCyclicLaneValue {batch : Nat}
    (residues : Fin batch → Fin 8 → Nat)
    (item : Fin batch) (column : Fin 8) : Nat :=
  r1CyclicLaneValue (residues item) column

theorem r1_batch_cyclic_lane_eq_matrix_product {batch : Nat}
    (residues : Fin batch → Fin 8 → Nat)
    (item : Fin batch) (column : Fin 8) :
    r1BatchCyclicLaneValue residues item column =
      matrixProductValue (residues item) Generated.r1_rns_mat column := by
  exact r1_cyclic_lane_eq_matrix_product (residues item) column

def lowMatrixAccumulator
    (residues : Fin 8 → Nat) (matrix : List (List Nat))
    (column : Fin 8) : Nat :=
  ∑ row, lowDigitProduct (residues row)
    (matrixEntry matrix row.val column.val)

def highMatrixAccumulator
    (residues : Fin 8 → Nat) (matrix : List (List Nat))
    (column : Fin 8) : Nat :=
  ∑ row, highDigitProduct (residues row)
    (matrixEntry matrix row.val column.val)

/--
Recombining the low and high IFMA accumulators gives the ordinary matrix
product when both source lanes and matrix digits fit the radix-`2^52` input
domain.  This is the source-level correspondence bridge for the two
`madd52lo`/`madd52hi` streams used by every fixed MatrixNoK matrix.
-/
theorem matrix_accumulators_reconstruct_product
    (residues : Fin 8 → Nat) (matrix : List (List Nat))
    (column : Fin 8)
    (residuesFit : ∀ row : Fin 8, residues row < ifmaRadix)
    (matrixFits : ∀ row : Fin 8,
      matrixEntry matrix row.val column.val < ifmaRadix) :
    lowMatrixAccumulator residues matrix column +
        ifmaRadix * highMatrixAccumulator residues matrix column =
      matrixProductValue residues matrix column := by
  unfold lowMatrixAccumulator highMatrixAccumulator matrixProductValue
  rw [Finset.mul_sum, ← Finset.sum_add_distrib]
  apply Finset.sum_congr rfl
  intro row _
  rw [digit_product_reconstruction]
  rw [Nat.mod_eq_of_lt (residuesFit row)]
  rw [Nat.mod_eq_of_lt (matrixFits row)]

theorem generated_r1_matrix_entry_lt_radix (row column : Fin 8) :
    matrixEntry Generated.r1_rns_mat row.val column.val < ifmaRadix := by
  fin_cases row <;> fin_cases column <;> decide

theorem generated_r2_matrix_entry_lt_radix (row column : Fin 8) :
    matrixEntry Generated.r2_rns_mat row.val column.val < ifmaRadix := by
  fin_cases row <;> fin_cases column <;> decide

/-- The fused M2 product and correction-free R1 streams reconstruct the
ordinary product plus the generated R1 matrix product exactly. -/
theorem r1_fused_accumulators_reconstruct
    (left right residues : Fin 8 → Nat) (column : Fin 8)
    (leftFits : left column < ifmaRadix)
    (rightFits : right column < ifmaRadix)
    (residuesFit : ∀ row, residues row < ifmaRadix) :
    lowDigitProduct (left column) (right column) +
          lowMatrixAccumulator residues Generated.r1_rns_mat column +
        ifmaRadix *
          (highDigitProduct (left column) (right column) +
            highMatrixAccumulator residues Generated.r1_rns_mat column) =
      left column * right column +
        matrixProductValue residues Generated.r1_rns_mat column := by
  have productReconstruction :=
    digit_product_reconstruction (left column) (right column)
  rw [Nat.mod_eq_of_lt leftFits, Nat.mod_eq_of_lt rightFits] at productReconstruction
  have matrixReconstruction := matrix_accumulators_reconstruct_product
    residues Generated.r1_rns_mat column residuesFit
    (fun row => generated_r1_matrix_entry_lt_radix row column)
  nlinarith

/-- The fused M2 product plus all eight R1 IFMA additions cannot wrap. -/
theorem r1_fused_ifma_accumulators_fit_word
    (left right residues : Fin 8 → Nat) (column : Fin 8) :
    lowDigitProduct (left column) (right column) +
          lowMatrixAccumulator residues Generated.r1_rns_mat column <
        uint64Modulus ∧
      highDigitProduct (left column) (right column) +
          highMatrixAccumulator residues Generated.r1_rns_mat column <
        uint64Modulus := by
  constructor
  · apply exact_multiplication_r1_accumulator_fits_word
    · exact (digit_products_fit (left column) (right column)).1
    · intro index
      exact (digit_products_fit (residues index)
        (matrixEntry Generated.r1_rns_mat index.val column.val)).1
  · apply exact_multiplication_r1_accumulator_fits_word
    · exact (digit_products_fit (left column) (right column)).2
    · intro index
      exact (digit_products_fit (residues index)
        (matrixEntry Generated.r1_rns_mat index.val column.val)).2

theorem r1_matrix_product_lt_16_square
    (residues : Fin 8 → Nat)
    (residuesFit : ∀ row,
      residues row < 2 * listAt8 Generated.moduli1 row.val)
    (column : Fin 8) :
    matrixProductValue residues Generated.r1_rns_mat column <
      16 * (listAt8 Generated.moduli2 column.val) ^ 2 := by
  have h0 := residuesFit (0 : Fin 8)
  have h1 := residuesFit (1 : Fin 8)
  have h2 := residuesFit (2 : Fin 8)
  have h3 := residuesFit (3 : Fin 8)
  have h4 := residuesFit (4 : Fin 8)
  have h5 := residuesFit (5 : Fin 8)
  have h6 := residuesFit (6 : Fin 8)
  have h7 := residuesFit (7 : Fin 8)
  fin_cases column <;>
    simp [matrixProductValue, matrixEntry, listAt8,
      Generated.moduli1, Generated.moduli2, Generated.r1_rns_mat,
      Fin.sum_univ_succ] at h0 h1 h2 h3 h4 h5 h6 h7 ⊢ <;>
    omega

theorem r2_matrix_product_lt_16_square
    (residues : Fin 8 → Nat)
    (residuesFit : ∀ row,
      residues row < 2 * listAt8 Generated.moduli2 row.val)
    (column : Fin 8) :
    matrixProductValue residues Generated.r2_rns_mat column <
      16 * (listAt8 Generated.moduli1 column.val) ^ 2 := by
  have h0 := residuesFit (0 : Fin 8)
  have h1 := residuesFit (1 : Fin 8)
  have h2 := residuesFit (2 : Fin 8)
  have h3 := residuesFit (3 : Fin 8)
  have h4 := residuesFit (4 : Fin 8)
  have h5 := residuesFit (5 : Fin 8)
  have h6 := residuesFit (6 : Fin 8)
  have h7 := residuesFit (7 : Fin 8)
  fin_cases column <;>
    simp [matrixProductValue, matrixEntry, listAt8,
      Generated.moduli1, Generated.moduli2, Generated.r2_rns_mat,
      Fin.sum_univ_succ] at h0 h1 h2 h3 h4 h5 h6 h7 ⊢ <;>
    omega

theorem generated_r2_estimate_lt_word
    (row : Fin 8) :
    listAt8 Generated.r2_shifted_quotient_estimations row.val <
      uint64Modulus := by
  fin_cases row <;> decide

theorem generated_r2_quotient_lt_2pow55
    (residues : Fin 8 → Nat)
    (residuesFit : ∀ row,
      residues row < 2 * listAt8 Generated.moduli2 row.val) :
    generatedR2EstimatedQuotient residues < 2 ^ 55 := by
  have scalarsFit : ∀ row, residues row < ifmaRadix := by
    intro row
    have modulusFit :
        2 * listAt8 Generated.moduli2 row.val < ifmaRadix := by
      fin_cases row <;> norm_num [Generated.moduli2, listAt8, ifmaRadix]
    exact (residuesFit row).trans modulusFit
  have words := quotient_words_fit residues
    (fun row => listAt8 Generated.r2_shifted_quotient_estimations row.val)
    scalarsFit generated_r2_estimate_lt_word
  simpa [generatedR2EstimatedQuotient, quotientNumerator] using words.2.1

def r2CorrectionShift (column : Fin 8) : Nat :=
  listAt8 Generated.r2_correction column.val * ifmaRadix %
    listAt8 Generated.moduli1 column.val

def generatedR2ShiftedWideValue
    (residues : Fin 8 → Nat) (column : Fin 8) : Nat :=
  let quotient := generatedR2EstimatedQuotient residues
  matrixProductValue residues Generated.r2_rns_mat column +
    listAt8 Generated.r2_correction column.val * (quotient % ifmaRadix) +
    r2CorrectionShift column * (quotient / ifmaRadix)

def generatedR2FinalLow
    (residues : Fin 8 → Nat) (column : Fin 8) : Nat :=
  let quotient := generatedR2EstimatedQuotient residues
  lowMatrixAccumulator residues Generated.r2_rns_mat column +
    lowDigitProduct
      (listAt8 Generated.r2_correction column.val) quotient +
    lowDigitProduct (r2CorrectionShift column) (quotient / ifmaRadix)

def generatedR2FinalHigh
    (residues : Fin 8 → Nat) (column : Fin 8) : Nat :=
  let quotient := generatedR2EstimatedQuotient residues
  highMatrixAccumulator residues Generated.r2_rns_mat column +
    highDigitProduct
      (listAt8 Generated.r2_correction column.val) quotient +
    highDigitProduct (r2CorrectionShift column) (quotient / ifmaRadix)

/-- The exact split accumulators computed by `reduce_raw_shifted` reconstruct
the mathematical shifted R2 wide value, with no implicit lane truncation. -/
theorem generated_r2_accumulators_reconstruct_wide
    (residues : Fin 8 → Nat)
    (residuesFit : ∀ row,
      residues row < 2 * listAt8 Generated.moduli2 row.val)
    (column : Fin 8) :
    generatedR2FinalLow residues column + ifmaRadix *
        generatedR2FinalHigh residues column =
      generatedR2ShiftedWideValue residues column := by
  have residuesDigitFit : ∀ row, residues row < ifmaRadix := by
    intro row
    exact (residuesFit row).trans (by
      fin_cases row <;>
        norm_num [Generated.moduli2, listAt8, ifmaRadix])
  have matrixReconstruction := matrix_accumulators_reconstruct_product
    residues Generated.r2_rns_mat column residuesDigitFit
    (fun row => generated_r2_matrix_entry_lt_radix row column)
  have correctionFits :
      listAt8 Generated.r2_correction column.val < ifmaRadix := by
    fin_cases column <;>
      norm_num [Generated.r2_correction, listAt8, ifmaRadix]
  have shiftFits : r2CorrectionShift column < ifmaRadix := by
    fin_cases column <;>
      norm_num [r2CorrectionShift, Generated.r2_correction,
        Generated.moduli1, listAt8, ifmaRadix]
  have quotientBound := generated_r2_quotient_lt_2pow55 residues residuesFit
  have quotientHighFits :
      generatedR2EstimatedQuotient residues / ifmaRadix < ifmaRadix := by
    apply (Nat.div_lt_iff_lt_mul (by norm_num [ifmaRadix])).2
    exact quotientBound.trans (by norm_num [ifmaRadix])
  have firstProduct := digit_product_reconstruction
    (listAt8 Generated.r2_correction column.val)
    (generatedR2EstimatedQuotient residues)
  rw [Nat.mod_eq_of_lt correctionFits] at firstProduct
  have secondProduct := digit_product_reconstruction
    (r2CorrectionShift column)
    (generatedR2EstimatedQuotient residues / ifmaRadix)
  rw [Nat.mod_eq_of_lt shiftFits,
    Nat.mod_eq_of_lt quotientHighFits] at secondProduct
  unfold generatedR2FinalLow generatedR2FinalHigh
    generatedR2ShiftedWideValue
  dsimp only
  nlinarith

/-- The exact shifted R2 high/low IFMA schedules cannot wrap 64-bit lanes. -/
theorem r2_shifted_ifma_accumulators_fit_word
    (residues : Fin 8 → Nat) (column : Fin 8) :
    let quotient := generatedR2EstimatedQuotient residues
    lowMatrixAccumulator residues Generated.r2_rns_mat column +
          lowDigitProduct
            (listAt8 Generated.r2_correction column.val) quotient +
          lowDigitProduct (r2CorrectionShift column)
            (quotient / ifmaRadix) < uint64Modulus ∧
      highMatrixAccumulator residues Generated.r2_rns_mat column +
          highDigitProduct
            (listAt8 Generated.r2_correction column.val) quotient +
          highDigitProduct (r2CorrectionShift column)
            (quotient / ifmaRadix) < uint64Modulus := by
  dsimp only
  apply shifted_matrix_correction_fits_word
  · unfold highMatrixAccumulator
    exact fin8_sum_lt (fun index =>
      (digit_products_fit (residues index)
        (matrixEntry Generated.r2_rns_mat index.val column.val)).2)
  · unfold lowMatrixAccumulator
    exact fin8_sum_lt (fun index =>
      (digit_products_fit (residues index)
        (matrixEntry Generated.r2_rns_mat index.val column.val)).1)

theorem r2_correction_and_shift_lt_modulus (column : Fin 8) :
    listAt8 Generated.r2_correction column.val <
        listAt8 Generated.moduli1 column.val ∧
      r2CorrectionShift column < listAt8 Generated.moduli1 column.val := by
  constructor
  · fin_cases column <;>
      norm_num [Generated.r2_correction, Generated.moduli1, listAt8]
  · unfold r2CorrectionShift
    exact Nat.mod_lt _ (by
      fin_cases column <;>
        norm_num [Generated.moduli1, listAt8])

theorem r2_shifted_correction_lt_5_square
    (residues : Fin 8 → Nat)
    (residuesFit : ∀ row,
      residues row < 2 * listAt8 Generated.moduli2 row.val)
    (column : Fin 8) :
    let quotient := generatedR2EstimatedQuotient residues
    listAt8 Generated.r2_correction column.val * (quotient % ifmaRadix) +
        r2CorrectionShift column * (quotient / ifmaRadix) <
      5 * (listAt8 Generated.moduli1 column.val) ^ 2 := by
  dsimp only
  have quotientBound := generated_r2_quotient_lt_2pow55 residues residuesFit
  have quotientLow :
      generatedR2EstimatedQuotient residues % ifmaRadix < ifmaRadix :=
    Nat.mod_lt _ (by norm_num [ifmaRadix])
  have quotientHigh :
      generatedR2EstimatedQuotient residues / ifmaRadix < 8 := by
    apply (Nat.div_lt_iff_lt_mul (by norm_num [ifmaRadix])).2
    norm_num [ifmaRadix]
    exact quotientBound
  have correctionBounds := r2_correction_and_shift_lt_modulus column
  have modulusWindow :
      ifmaRadix + 8 < 5 * listAt8 Generated.moduli1 column.val := by
    fin_cases column <;>
      norm_num [Generated.moduli1, listAt8, ifmaRadix]
  nlinarith

theorem generated_r2_shifted_wide_lt_32_square
    (residues : Fin 8 → Nat)
    (residuesFit : ∀ row,
      residues row < 2 * listAt8 Generated.moduli2 row.val)
    (column : Fin 8) :
    generatedR2ShiftedWideValue residues column <
      32 * (listAt8 Generated.moduli1 column.val) ^ 2 := by
  have matrixBound := r2_matrix_product_lt_16_square
    residues residuesFit column
  have correctionBound := r2_shifted_correction_lt_5_square
    residues residuesFit column
  unfold generatedR2ShiftedWideValue
  nlinarith

theorem generated_r2_shifted_wide_represents_lane
    (residues : Fin 8 → Nat) (column : Fin 8) :
    (generatedR2ShiftedWideValue residues column :
        ZMod (listAt8 Generated.moduli1 column.val)) =
      generatedR2LaneValue residues column := by
  unfold generatedR2ShiftedWideValue generatedR2LaneValue
  unfold matrixProductValue
  push_cast
  have shiftCast :
      (r2CorrectionShift column :
          ZMod (listAt8 Generated.moduli1 column.val)) =
        listAt8 Generated.r2_correction column.val * ifmaRadix := by
    unfold r2CorrectionShift
    rw [ZMod.natCast_mod]
    push_cast
    rfl
  have quotientSplit := Nat.mod_add_div
    (generatedR2EstimatedQuotient residues) ifmaRadix
  have quotientSplitCast :
      (generatedR2EstimatedQuotient residues :
          ZMod (listAt8 Generated.moduli1 column.val)) =
        (generatedR2EstimatedQuotient residues % ifmaRadix :
          ZMod (listAt8 Generated.moduli1 column.val)) +
        ifmaRadix *
          (generatedR2EstimatedQuotient residues / ifmaRadix :
            ZMod (listAt8 Generated.moduli1 column.val)) := by
    have castEquality := congrArg
      (fun value : Nat =>
        (value : ZMod (listAt8 Generated.moduli1 column.val)))
      quotientSplit.symm
    push_cast at castEquality
    exact castEquality
  rw [shiftCast, quotientSplitCast]
  ring

end CatBLST.AVX.MatrixNoK
