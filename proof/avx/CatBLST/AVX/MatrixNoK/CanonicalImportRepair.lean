/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.CanonicalImport
import Mathlib.Tactic

/-!
# Canonical import repair witness

This is the exact canonical value that exposed `AUD-221`. It selected a
negative representative with VROOM's ceiling coefficients. The repaired
floor coefficients select the positive representative `p - t`, preserving
the field value while satisfying the following nonnegative R2 contract.
-/

namespace CatBLST.AVX.MatrixNoK

set_option exponentiation.threshold 512

def repairedImportMagnitude : Nat :=
  235613689341516701328257298793030861480281319249132620885720573137412674019990281389671108000420691382528538236

def repairedImportCanonicalValue : Nat :=
  1267253245981839142264863935160608057326608172767180188451362181929717619810821363803276830181705476116338914753674

def repairedImportDigits (index : Fin 8) : Nat :=
  listAt8 [
    500616628141194,
    1029601075259562,
    686812189624822,
    769198589660543,
    581522751939018,
    26173626664305,
    89497429892085,
    552541535
  ] index.val

theorem repaired_import_digits_fit (index : Fin 8) :
    repairedImportDigits index < 2 ^ 50 := by
  fin_cases index <;>
    norm_num [repairedImportDigits, listAt8]

theorem repaired_import_digits_encode_canonical_value :
    canonicalDigitValue repairedImportDigits =
      repairedImportCanonicalValue := by
  norm_num [canonicalDigitValue, repairedImportDigits, listAt8,
    repairedImportCanonicalValue, Fin.sum_univ_succ]

theorem repaired_import_canonical_value_is_valid :
    repairedImportCanonicalValue < Generated.target_modulus := by
  norm_num [repairedImportCanonicalValue, Generated.target_modulus]

theorem repaired_import_selects_positive_representative :
    convertToCorrectedRepresentative repairedImportDigits =
      Generated.target_modulus - repairedImportMagnitude := by
  norm_num [convertToCorrectedRepresentative,
    floorCorrectedCrtRepresentative, exactCrtNumerator,
    floorEstimatedCrtQuotient, approximateCrtNumeratorFloor,
    repairedImportDigits, repairedImportMagnitude, convertToShiftedFactor,
    uint64Modulus, Generated.convert_to_premult, Generated.target_modulus,
    listAt8, Fin.sum_univ_succ]

theorem repaired_import_representative_is_in_r2_range :
    convertToCorrectedRepresentative repairedImportDigits <
      2 * Generated.target_modulus :=
  convert_to_corrected_range repairedImportDigits repaired_import_digits_fit

theorem repaired_import_represents_canonical_value :
    (convertToCorrectedRepresentative repairedImportDigits :
        ZMod Generated.target_modulus) =
      Generated.convert_to_premult * repairedImportCanonicalValue := by
  rw [convert_to_corrected_represents_digits]
  rw [repaired_import_digits_encode_canonical_value]

end CatBLST.AVX.MatrixNoK
