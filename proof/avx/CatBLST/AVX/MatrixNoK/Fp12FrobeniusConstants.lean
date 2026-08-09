/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Fp12Frobenius3
import CatBLST.AVX.MatrixNoK.Pipeline

/-!
# Generated Fp12 Frobenius constants

This file binds each checked-in M1/M2 lane array to the exact standard-element
encoding used by `Fp12ConstantsCandidate::standard`.  The intermediate
representative is the canonical coefficient imported into MatrixNoK's
Montgomery domain; the two lane bases then satisfy the accepted P2
multiplication preconditions.
-/

namespace CatBLST.AVX.MatrixNoK

open Fp12Concrete

def frobeniusCanonicalValue : Fin 10 → Nat := ![
  0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaac,
  0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaad,
  0x5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a00022e01fffffffefffe,
  0x1904d3bf02bb0667c231beb4202c0d1f0fd603fd3cbd5f4f7b2443d784bab9c4f67ea53d63e7813d8d0775ed92235fb8,
  0xfc3e2b36c4e03288e9e902231f9fb854a14787b6c7b36fec0c8ec971f63c5f282d5ac14d6c7ec22cf78a126ddc4af3,
  0x6af0e0437ff400b6831e36d6bd17ffe48395dabc2d3435e77f76e17009241c5ee67992f72ec05f4c81084fbede3cc09,
  0x135203e60180a68ee2e9c448d77a2cd91c3dedd930b1cf60ef396489f61eb45e304466cf3e67fa0af1ee7b04121bdea2,
  0x5b2cfd9013a5fd8df47fa6b48b1e045f39816240c0b8fee8beadf4d8e9c0566c63a3e6e257f87329b18fae980078116,
  0x144e4211384586c16bd3ad4afa99cc9170df3560e77982d0db45f3536814f0bd5871c1908bd478cd1ee605167ff82995,
  0x5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a00022e01fffffffeffff]

def frobeniusSemanticConstant : Fin 10 → SelectedBaseField := ![
  frob1_2, frob1_3, frob1_4, frob1_5, frob1_6,
  frob1_7, frob1_8, frob1_9, frob1_10, frob2_3]

def generatedFrobeniusM1 : Fin 10 → List Nat := ![
  Generated.fp12_frob1_2_m1, Generated.fp12_frob1_3_m1,
  Generated.fp12_frob1_4_m1, Generated.fp12_frob1_5_m1,
  Generated.fp12_frob1_6_m1, Generated.fp12_frob1_7_m1,
  Generated.fp12_frob1_8_m1, Generated.fp12_frob1_9_m1,
  Generated.fp12_frob1_10_m1, Generated.fp12_frob2_3_m1]

def generatedFrobeniusM2 : Fin 10 → List Nat := ![
  Generated.fp12_frob1_2_m2, Generated.fp12_frob1_3_m2,
  Generated.fp12_frob1_4_m2, Generated.fp12_frob1_5_m2,
  Generated.fp12_frob1_6_m2, Generated.fp12_frob1_7_m2,
  Generated.fp12_frob1_8_m2, Generated.fp12_frob1_9_m2,
  Generated.fp12_frob1_10_m2, Generated.fp12_frob2_3_m2]

def generatedFrobeniusM1Lane (constant : Fin 10) (lane : Fin 8) : Nat :=
  listAt8 (generatedFrobeniusM1 constant) lane.val

def generatedFrobeniusM2Lane (constant : Fin 10) (lane : Fin 8) : Nat :=
  listAt8 (generatedFrobeniusM2 constant) lane.val

/-- The intermediate representative constructed by `to_mont_exact`. -/
def frobeniusStandardRepresentative (constant : Fin 10) : Nat :=
  let value := frobeniusCanonicalValue constant
  value * Generated.M1 % Generated.target_modulus +
    value / Generated.target_modulus * Generated.target_modulus

/-- The common rotated RNS integer constructed by `to_mont_exact`. -/
def frobeniusExactStandardEncoding (constant : Fin 10) : Nat :=
  frobeniusStandardRepresentative constant * Generated.total_rotation %
    Generated.M1M2

theorem frobenius_semantic_constant_matches : ∀ constant : Fin 10,
    (frobeniusCanonicalValue constant : SelectedBaseField) =
      frobeniusSemanticConstant constant := by
  native_decide

theorem generated_frobenius_exact_encoding :
    ∀ (constant : Fin 10) (lane : Fin 8),
      generatedFrobeniusM1Lane constant lane =
          frobeniusExactStandardEncoding constant % modulus1At lane ∧
        generatedFrobeniusM2Lane constant lane =
          frobeniusExactStandardEncoding constant % modulus2At lane := by
  native_decide

theorem generated_frobenius_standard_contract :
    ∀ (constant : Fin 10) (lane : Fin 8),
      generatedFrobeniusM1Lane constant lane < modulus1At lane ∧
        Nat.ModEq (modulus1At lane)
          (generatedFrobeniusM1Lane constant lane)
          (frobeniusStandardRepresentative constant * ifmaRadix *
            Generated.qr_premult) ∧
        generatedFrobeniusM2Lane constant lane < modulus2At lane ∧
        Nat.ModEq (modulus2At lane)
          (generatedFrobeniusM2Lane constant lane)
          (frobeniusStandardRepresentative constant * ifmaRadix *
            crossBasisInverse) := by
  native_decide

theorem frobenius_standard_representative_semantics :
    ∀ constant : Fin 10,
      frobeniusStandardRepresentative constant < Generated.target_modulus ∧
        (frobeniusStandardRepresentative constant : SelectedBaseField) =
          (Generated.M1 : SelectedBaseField) *
            frobeniusSemanticConstant constant := by
  native_decide

end CatBLST.AVX.MatrixNoK
