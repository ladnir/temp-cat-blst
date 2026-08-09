/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Constants

/-!
# Selected GT generator constant

This file binds the twelve Fp coefficients selected by production `GT.h` to
the flattened M1/M2 MatrixNoK literal. The coefficient values are in CatBLST's
internal Fp12 order. Establishing that this canonical Fp12 value is CatBLST's
chosen GT generator remains a reviewed scalar-constant boundary and is also
covered by the scalar/AVX differential test.
-/

namespace CatBLST.AVX.MatrixNoK

def gtGeneratorCanonicalCoefficient : Fin 12 -> Nat := ![
  0x1250ebd871fc0a92a7b2d83168d0d727272d441befa15c503dd8e90ce98db3e7b6d194f60839c508a84305aaca1789b6,
  0x089a1c5b46e5110b86750ec6a532348868a84045483c92b7af5af689452eafabf1a8943e50439f1d59882a98eaa0170f,
  0x1368bb445c7c2d209703f239689ce34c0378a68e72a6b3b216da0e22a5031b54ddff57309396b38c881c4c849ec23e87,
  0x193502b86edb8857c273fa075a50512937e0794e1e65a7617c90d8bd66065b1fffe51d7a579973b1315021ec3c19934f,
  0x01b2f522473d171391125ba84dc4007cfbf2f8da752f7c74185203fcca589ac719c34dffbbaad8431dad1c1fb597aaa5,
  0x018107154f25a764bd3c79937a45b84546da634b8f6be14a8061e55cceba478b23f7dacaa35c8ca78beae9624045b4b6,
  0x19f26337d205fb469cd6bd15c3d5a04dc88784fbb3d0b2dbdea54d43b2b73f2cbb12d58386a8703e0f948226e47ee89d,
  0x06fba23eb7c5af0d9f80940ca771b6ffd5857baaf222eb95a7d2809d61bfe02e1bfd1b68ff02f0b8102ae1c2d5d5ab1a,
  0x11b8b424cd48bf38fcef68083b0b0ec5c81a93b330ee1a677d0d15ff7b984e8978ef48881e32fac91b93b47333e2ba57,
  0x03350f55a7aefcd3c31b4fcb6ce5771cc6a0e9786ab5973320c806ad360829107ba810c5a09ffdd9be2291a0c25a99a2,
  0x04c581234d086a9902249b64728ffd21a189e87935a954051c7cdba7b3872629a4fafc05066245cb9108f0242d0fe3ef,
  0x0f41e58663bf08cf068672cbd01a7ec73baca4d72ca93544deff686bfd6df543d48eaa24afe47e1efde449383b676631]

def generatedGtM1Lane (coefficient : Fin 12) (lane : Fin 8) : Nat :=
  Generated.gt_generator.getD (coefficient.val * 16 + lane.val) 0

def generatedGtM2Lane (coefficient : Fin 12) (lane : Fin 8) : Nat :=
  Generated.gt_generator.getD (coefficient.val * 16 + 8 + lane.val) 0

def gtGeneratorModulus1At (lane : Fin 8) : Nat :=
  Generated.moduli1.getD lane.val 0

def gtGeneratorModulus2At (lane : Fin 8) : Nat :=
  Generated.moduli2.getD lane.val 0

def gtGeneratorStandardRepresentative (coefficient : Fin 12) : Nat :=
  let value := gtGeneratorCanonicalCoefficient coefficient
  value * Generated.M1 % Generated.target_modulus +
    value / Generated.target_modulus * Generated.target_modulus

def gtGeneratorExactStandardEncoding (coefficient : Fin 12) : Nat :=
  gtGeneratorStandardRepresentative coefficient * Generated.total_rotation %
    Generated.M1M2

theorem generated_gt_generator_shape :
    Generated.gt_generator.length = 192 := by
  native_decide

theorem gt_generator_coefficients_are_canonical :
    forall coefficient : Fin 12,
      gtGeneratorCanonicalCoefficient coefficient < Generated.target_modulus := by
  native_decide

theorem generated_gt_generator_exact_encoding :
    forall (coefficient : Fin 12) (lane : Fin 8),
      generatedGtM1Lane coefficient lane =
          gtGeneratorExactStandardEncoding coefficient %
            gtGeneratorModulus1At lane /\
        generatedGtM2Lane coefficient lane =
          gtGeneratorExactStandardEncoding coefficient %
            gtGeneratorModulus2At lane := by
  native_decide

theorem generated_gt_generator_standard_contract :
    forall (coefficient : Fin 12) (lane : Fin 8),
      generatedGtM1Lane coefficient lane < gtGeneratorModulus1At lane /\
        Nat.ModEq (gtGeneratorModulus1At lane)
          (generatedGtM1Lane coefficient lane)
          (gtGeneratorStandardRepresentative coefficient * montgomeryRadix *
            Generated.qr_premult) /\
        generatedGtM2Lane coefficient lane < gtGeneratorModulus2At lane /\
        Nat.ModEq (gtGeneratorModulus2At lane)
          (generatedGtM2Lane coefficient lane)
          (gtGeneratorStandardRepresentative coefficient * montgomeryRadix *
            crossBasisInverse) := by
  native_decide

end CatBLST.AVX.MatrixNoK
