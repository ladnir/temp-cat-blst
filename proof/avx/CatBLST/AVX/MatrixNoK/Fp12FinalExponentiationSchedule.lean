/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Fp12Cyclotomic
import Mathlib.Tactic

set_option maxHeartbeats 800000

/-!
# MatrixNoK final-exponentiation schedule

This file transcribes the exponent flow of
`blst/avx/fields/final_exponentiation.h`.  Exponents are signed because every
conjugation after the easy part is inversion on the cyclotomic subgroup.  The
first conjugation remains the general Fp12 Frobenius-six map and is represented
by the `p^6` term in `finalEasyExponent`.

BLST's production pairing convention is intentionally cubed: the complete
schedule computes `f ^ (3 * ((p^12 - 1) / r))`, not the uncubed textbook final
exponent.  Since three is invertible modulo the prime subgroup order, this is
still a nondegenerate pairing into the same GT subgroup.
-/

namespace CatBLST.AVX.MatrixNoK

/-- Absolute value of the negative BLS12-381 parameter used by the source. -/
def blsZMagnitude : Int := 0xd201000000010000

/-- Signed BLS12-381 parameter. -/
def blsZ : Int := -blsZMagnitude

/-- Prime-order subgroup modulus. -/
def blsScalarOrder : Nat :=
  0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001

/-- Exponent update performed by `mul_n_sqr(a, b, n)`. -/
def mulNSqrExponent (a b : Int) (n : Nat) : Int :=
  (a + b) * (2 ^ n : Nat)

/-- Literal exponent trace of `raise_to_z_div_by_2`. -/
def raiseToZDivBy2Exponent (a : Int) : Int :=
  let e0 := 2 * a
  let e1 := mulNSqrExponent e0 a 2
  let e2 := mulNSqrExponent e1 a 3
  let e3 := mulNSqrExponent e2 a 9
  let e4 := mulNSqrExponent e3 a 32
  let e5 := mulNSqrExponent e4 a (16 - 1)
  0 - e5

/-- Literal exponent trace of `raise_to_z`. -/
def raiseToZExponent (a : Int) : Int :=
  2 * raiseToZDivBy2Exponent a

theorem raise_to_z_div_by_2_exponent_correct (a : Int) :
    raiseToZDivBy2Exponent a = (blsZ / 2) * a := by
  norm_num [raiseToZDivBy2Exponent, mulNSqrExponent, blsZ, blsZMagnitude]
  ring

theorem raise_to_z_exponent_correct (a : Int) :
    raiseToZExponent a = blsZ * a := by
  rw [raiseToZExponent, raise_to_z_div_by_2_exponent_correct]
  norm_num [blsZ, blsZMagnitude]
  ring

/-- Easy-part exponent after inverse, Frobenius two, and the second multiply. -/
def finalEasyExponent : Int :=
  ((Generated.target_modulus : Int) ^ 6 - 1) *
    ((Generated.target_modulus : Int) ^ 2 + 1)

/--
Literal hard-part trace beginning with `ret2`.  Every local name and operation
matches `final_exp_after_inverse`; multiplication adds exponents, cyclotomic
squaring doubles them, Frobenius `k` multiplies them by `p^k`, and conjugation
negates them on the cyclotomic domain.
-/
def finalHardScheduleExponent (ret2 : Int) : Int :=
  let ret2Sqr := 2 * ret2
  let y0ToZ := raiseToZExponent ret2Sqr
  let y0ToZDiv2 := raiseToZDivBy2Exponent y0ToZ
  let y1MulY3 := -ret2 + y0ToZ
  let y1Final := -y1MulY3 + y0ToZDiv2
  let y1ToZ := raiseToZExponent y1Final
  let y2ToZ := raiseToZExponent y1ToZ
  let y3MulY1 := y2ToZ - y1Final
  let y1MulY2 :=
    ((Generated.target_modulus : Int) ^ 3) * y1Final +
      ((Generated.target_modulus : Int) ^ 2) * y1ToZ
  let y3ToZ := raiseToZExponent y3MulY1
  let y2MulY0 := y3ToZ + ret2Sqr
  let y2MulRet := y2MulY0 + ret2
  let y1MulY2Final := y1MulY2 + y2MulRet
  let y3Frob1 := (Generated.target_modulus : Int) * y3MulY1
  y1MulY2Final + y3Frob1

/-- Full literal exponent trace, including the easy part. -/
def cppFinalExponentiationExponent : Int :=
  finalHardScheduleExponent finalEasyExponent

/-- BLST's documented cubed final exponent. -/
def cubedFinalExponent : Int :=
  3 * (((Generated.target_modulus ^ 12 - 1) / blsScalarOrder : Nat) : Int)

theorem bls_z_is_even : blsZ % 2 = 0 := by
  norm_num [blsZ, blsZMagnitude]

theorem bls_scalar_order_parameter_identity :
    blsScalarOrder =
      Int.natAbs (blsZ ^ 4 - blsZ ^ 2 + 1) := by
  native_decide

theorem target_modulus_parameter_identity :
    (Generated.target_modulus : Int) =
      ((blsZ - 1) ^ 2 * (blsZ ^ 4 - blsZ ^ 2 + 1)) / 3 + blsZ := by
  norm_num [Generated.target_modulus, blsZ, blsZMagnitude]

theorem scalar_order_divides_target_multiplicative_order :
    Dvd.dvd blsScalarOrder (Generated.target_modulus ^ 12 - 1) := by
  native_decide

theorem final_easy_exponent_factorization :
    finalEasyExponent =
      (((Generated.target_modulus : Int) ^ 6 - 1) *
        ((Generated.target_modulus : Int) ^ 2 + 1)) := rfl

theorem cpp_final_exponentiation_exponent_correct :
    cppFinalExponentiationExponent = cubedFinalExponent := by
  native_decide

end CatBLST.AVX.MatrixNoK
