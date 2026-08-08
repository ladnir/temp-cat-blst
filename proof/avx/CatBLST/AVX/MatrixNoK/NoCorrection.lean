/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.RNS
import Mathlib.Tactic

/-!
# MatrixNoK correction-free R1 semantics

The R1 matrix is deliberately evaluated without the usual approximate-CRT
quotient correction.  This file isolates the reason that is valid: the R1
postmultiplier contains the target field modulus.  Changing the chosen R1
lift by any multiple of the input-basis modulus therefore changes the next
field representative only by a multiple of the target modulus.
-/

namespace CatBLST.AVX.MatrixNoK

theorem r1_lift_relation_generic
    (modulus target radix qr premult left right : Nat)
    (modulusAboveOne : 1 < modulus)
    (inverse :
      (qr * qr * premult * (modulus - target) * radix) % modulus = 1) :
    ((premult * (left * right * radix * qr * qr % modulus)) % modulus *
        (modulus - target)) % modulus =
      left * right % modulus := by
  change Nat.ModEq modulus
    ((premult * (left * right * radix * qr * qr % modulus)) % modulus *
      (modulus - target))
    (left * right)
  have inverseModEq : Nat.ModEq modulus
      (qr * qr * premult * (modulus - target) * radix) 1 := by
    change
      (qr * qr * premult * (modulus - target) * radix) % modulus =
        1 % modulus
    rw [Nat.mod_eq_of_lt modulusAboveOne]
    exact inverse
  have reducedProduct : Nat.ModEq modulus
      (left * right * radix * qr * qr % modulus)
      (left * right * radix * qr * qr) := by
    change
      (left * right * radix * qr * qr % modulus) % modulus =
      left * right * radix * qr * qr % modulus
    rw [Nat.mod_mod]
  have liftBase :=
    (Nat.ModEq.refl (n := modulus) premult).mul reducedProduct
  have liftReduced : Nat.ModEq modulus
      ((premult * (left * right * radix * qr * qr % modulus)) % modulus)
      (premult * (left * right * radix * qr * qr % modulus)) := by
    change
      (premult * (left * right * radix * qr * qr % modulus) % modulus) %
          modulus =
        (premult * (left * right * radix * qr * qr % modulus)) % modulus
    rw [Nat.mod_mod]
  have liftRelation :=
    (liftReduced.trans liftBase).mul
      (Nat.ModEq.refl (n := modulus) (modulus - target))
  have product :=
    (Nat.ModEq.refl (n := modulus) (left * right)).mul inverseModEq
  have reordered :
      premult * (left * right * radix * qr * qr) * (modulus - target) =
      left * right * (qr * qr * premult * (modulus - target) * radix) := by
    ring
  rw [reordered] at liftRelation
  rw [Nat.mul_one] at product
  exact liftRelation.trans product

/-- The canonical M1 residue selected by the reduced product and R1 premultiplier. -/
def r1CanonicalLift (left right : Nat) : Nat :=
  (Generated.r1_premult *
      (left * right * montgomeryRadix * Generated.qr_premult *
        Generated.qr_premult % Generated.M1)) % Generated.M1

theorem r1_numerator_divisible_generic
    (modulus target left right lift : Nat)
    (targetAtMost : target ≤ modulus)
    (liftRelation :
      (lift * (modulus - target)) % modulus = left * right % modulus) :
    modulus ∣ left * right + target * lift := by
  rw [Nat.dvd_iff_mod_eq_zero]
  rw [Nat.add_mod, ← liftRelation]
  rw [← Nat.add_mod]
  have identity :
      lift * (modulus - target) + target * lift = lift * modulus := by
    rw [Nat.mul_sub_left_distrib]
    rw [Nat.mul_comm target lift]
    exact Nat.sub_add_cancel (Nat.mul_le_mul_left lift targetAtMost)
  rw [identity]
  simp

theorem r1_canonical_lift_relation (left right : Nat) :
    (r1CanonicalLift left right *
      (Generated.M1 - Generated.target_modulus)) % Generated.M1 =
        left * right % Generated.M1 := by
  unfold r1CanonicalLift
  exact r1_lift_relation_generic
    Generated.M1 Generated.target_modulus montgomeryRadix
    Generated.qr_premult Generated.r1_premult left right
    (by norm_num [Generated.M1])
    r1_premult_relation

/-- The correction-free R1 lift makes the product numerator divisible by M1. -/
theorem r1_canonical_numerator_divisible (left right : Nat) :
    Generated.M1 ∣
      left * right + Generated.target_modulus * r1CanonicalLift left right := by
  have targetBelow : Generated.target_modulus ≤ Generated.M1 := by
    norm_num [Generated.target_modulus, Generated.M1]
  exact r1_numerator_divisible_generic
    Generated.M1 Generated.target_modulus left right
    (r1CanonicalLift left right) targetBelow
    (r1_canonical_lift_relation left right)

theorem r1_raw_numerator_divisible
    (left right rawLift : Nat)
    (rawRelation :
      rawLift % Generated.M1 = r1CanonicalLift left right) :
    Generated.M1 ∣
      left * right + Generated.target_modulus * rawLift := by
  have targetBelow : Generated.target_modulus ≤ Generated.M1 := by
    norm_num [Generated.target_modulus, Generated.M1]
  apply r1_numerator_divisible_generic
    Generated.M1 Generated.target_modulus left right rawLift targetBelow
  rw [Nat.mul_mod]
  rw [rawRelation]
  exact r1_canonical_lift_relation left right

/-- The next internal representative produced by the canonical R1 lift. -/
def r1ProductRepresentative (left right : Nat) : Nat :=
  (left * right + Generated.target_modulus * r1CanonicalLift left right) /
    Generated.M1

theorem r1_product_reconstruction (left right : Nat) :
    r1ProductRepresentative left right * Generated.M1 =
      left * right + Generated.target_modulus * r1CanonicalLift left right := by
  unfold r1ProductRepresentative
  exact Nat.div_mul_cancel (r1_canonical_numerator_divisible left right)

theorem r1_product_is_field_multiplication (left right : Nat) :
    (r1ProductRepresentative left right * Generated.M1) %
        Generated.target_modulus =
      left * right % Generated.target_modulus := by
  rw [r1_product_reconstruction]
  simp [Nat.add_mod]

/-- Omitting a quotient correction changes only the redundant multiple of p. -/
theorem skipped_r1_correction_preserves_field_value
    (representative quotient : Nat) :
    (representative + quotient * Generated.target_modulus) %
        Generated.target_modulus =
      representative % Generated.target_modulus := by
  simp [Nat.add_mod]

theorem cross_basis_output_generic
    (inputModulus outputModulus cross radix numerator representative : Nat)
    (outputAboveOne : 1 < outputModulus)
    (reconstruction : representative * inputModulus = numerator)
    (crossInverse : inputModulus * cross % outputModulus = 1) :
    (numerator * radix * cross * cross) % outputModulus =
      (representative * radix * cross) % outputModulus := by
  rw [← reconstruction]
  change Nat.ModEq outputModulus
    (representative * inputModulus * radix * cross * cross)
    (representative * radix * cross)
  have inverse : Nat.ModEq outputModulus (inputModulus * cross) 1 := by
    change (inputModulus * cross) % outputModulus = 1 % outputModulus
    rw [Nat.mod_eq_of_lt outputAboveOne]
    exact crossInverse
  have multiplied :=
    (Nat.ModEq.refl (n := outputModulus)
      (representative * radix * cross)).mul inverse
  have reordered :
      representative * inputModulus * radix * cross * cross =
        representative * radix * cross * (inputModulus * cross) := by
    ring
  rw [reordered]
  simpa only [Nat.mul_one] using multiplied

theorem r1_product_cross_basis_output (left right : Nat) :
    ((left * right +
        Generated.target_modulus * r1CanonicalLift left right) *
      montgomeryRadix * crossBasisInverse * crossBasisInverse) %
        Generated.M2 =
      (r1ProductRepresentative left right * montgomeryRadix *
        crossBasisInverse) % Generated.M2 := by
  exact cross_basis_output_generic
    Generated.M1 Generated.M2 crossBasisInverse montgomeryRadix
    (left * right + Generated.target_modulus * r1CanonicalLift left right)
    (r1ProductRepresentative left right)
    (by norm_num [Generated.M2])
    (r1_product_reconstruction left right)
    cross_basis_inverse_relation

end CatBLST.AVX.MatrixNoK
