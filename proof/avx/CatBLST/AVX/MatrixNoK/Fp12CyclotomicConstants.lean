/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Fp12Cyclotomic
import CatBLST.AVX.MatrixNoK.Pipeline

/-!
# Generated cyclotomic-square constants

This file binds both checked-in M1/M2 lane arrays used by the selected split
schedule to the exact standard-element encoding produced by the MatrixNoK
generator.  The semantic values are the canonical `1/3` and `-2/3` constants
used by the functional proof.
-/

namespace CatBLST.AVX.MatrixNoK

def cyclotomicCanonicalValue : Fin 2 → Nat := ![
  0x11560bf17baa99bc32126fced787c88f984f87adf7ae0c7f9a208c6b4f20a4181472aaa9cb8d555526a9ffffffffc71d,
  0x11560bf17baa99bc32126fced787c88f984f87adf7ae0c7f9a208c6b4f20a4181472aaa9cb8d555526a9ffffffffc71c]

def cyclotomicSemanticConstant : Fin 2 → SelectedBaseField := ![
  cyclotomicOneThird, cyclotomicNegTwoThird]

def generatedCyclotomicM1 : Fin 2 → List Nat := ![
  Generated.fp12_one_third_m1, Generated.fp12_neg_two_thirds_m1]

def generatedCyclotomicM2 : Fin 2 → List Nat := ![
  Generated.fp12_one_third_m2, Generated.fp12_neg_two_thirds_m2]

def generatedCyclotomicM1Lane (constant : Fin 2) (lane : Fin 8) : Nat :=
  listAt8 (generatedCyclotomicM1 constant) lane.val

def generatedCyclotomicM2Lane (constant : Fin 2) (lane : Fin 8) : Nat :=
  listAt8 (generatedCyclotomicM2 constant) lane.val

def cyclotomicStandardRepresentative (constant : Fin 2) : Nat :=
  let value := cyclotomicCanonicalValue constant
  value * Generated.M1 % Generated.target_modulus +
    value / Generated.target_modulus * Generated.target_modulus

def cyclotomicExactStandardEncoding (constant : Fin 2) : Nat :=
  cyclotomicStandardRepresentative constant * Generated.total_rotation %
    Generated.M1M2

theorem cyclotomic_semantic_constant_matches : ∀ constant : Fin 2,
    (cyclotomicCanonicalValue constant : SelectedBaseField) =
      cyclotomicSemanticConstant constant := by
  native_decide

theorem generated_cyclotomic_exact_encoding :
    ∀ (constant : Fin 2) (lane : Fin 8),
      generatedCyclotomicM1Lane constant lane =
          cyclotomicExactStandardEncoding constant % modulus1At lane ∧
        generatedCyclotomicM2Lane constant lane =
          cyclotomicExactStandardEncoding constant % modulus2At lane := by
  native_decide

theorem generated_cyclotomic_standard_contract :
    ∀ (constant : Fin 2) (lane : Fin 8),
      generatedCyclotomicM1Lane constant lane < modulus1At lane ∧
        Nat.ModEq (modulus1At lane)
          (generatedCyclotomicM1Lane constant lane)
          (cyclotomicStandardRepresentative constant * ifmaRadix *
            Generated.qr_premult) ∧
        generatedCyclotomicM2Lane constant lane < modulus2At lane ∧
        Nat.ModEq (modulus2At lane)
          (generatedCyclotomicM2Lane constant lane)
          (cyclotomicStandardRepresentative constant * ifmaRadix *
            crossBasisInverse) := by
  native_decide

theorem cyclotomic_standard_representative_semantics :
    ∀ constant : Fin 2,
      cyclotomicStandardRepresentative constant < Generated.target_modulus ∧
        (cyclotomicStandardRepresentative constant : SelectedBaseField) =
          (Generated.M1 : SelectedBaseField) *
            cyclotomicSemanticConstant constant := by
  native_decide

end CatBLST.AVX.MatrixNoK
