/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Fp12Conjugate
import CatBLST.AVX.MatrixNoK.Fp12FinalExponentiationSchedule
import Mathlib.Tactic

set_option maxHeartbeats 800000

/-!
# Algebra of the literal final-exponentiation schedule

The definitions below preserve the source loop and local-value order while
abstracting only the already accepted field operations.  A commutative group
is the precise algebraic domain after the easy part: the nonzero Miller result
and its supplied inverse produce a unit, and every later value remains a unit
in the cyclotomic subgroup.
-/

namespace CatBLST.AVX.MatrixNoK

variable {G : Type*} [CommGroup G]

/-- Literal `for` loop of repeated cyclotomic squaring. -/
def cppRepeatedSquare : Nat -> G -> G
  | 0, value => value
  | n + 1, value => cppRepeatedSquare n (value * value)

theorem cpp_repeated_square_eq_pow (n : Nat) (value : G) :
    cppRepeatedSquare n value = value ^ (2 ^ n) := by
  induction n generalizing value with
  | zero => simp [cppRepeatedSquare]
  | succ n ih =>
      rw [cppRepeatedSquare, ih, mul_pow]
      rw [(pow_add value (2 ^ n) (2 ^ n)).symm]
      congr 1
      omega

theorem cpp_repeated_square_zpow (base : G) (exponent : Int) (n : Nat) :
    cppRepeatedSquare n (base ^ exponent) =
      base ^ (exponent * (2 ^ n : Nat)) := by
  rw [cpp_repeated_square_eq_pow]
  rw [show (base ^ exponent) ^ (2 ^ n) =
      (base ^ exponent) ^ ((2 ^ n : Nat) : Int) from
        (zpow_natCast (base ^ exponent) (2 ^ n)).symm]
  rw [(zpow_mul base exponent ((2 ^ n : Nat) : Int)).symm]

/-- Literal multiply followed by `n` cyclotomic squares. -/
def cppGroupMulNSqr (a b : G) (n : Nat) : G :=
  cppRepeatedSquare n (a * b)

theorem cpp_group_mul_n_sqr_zpow
    (base : G) (a b : Int) (n : Nat) :
    cppGroupMulNSqr (base ^ a) (base ^ b) n =
      base ^ (mulNSqrExponent a b n) := by
  rw [cppGroupMulNSqr, (zpow_add base a b).symm,
    cpp_repeated_square_zpow]
  rfl

theorem inverse_zpow_base (base : G) (exponent : Int) :
    Inv.inv (base ^ exponent) = base ^ (-exponent) :=
  (zpow_neg base exponent).symm

theorem multiply_zpow_base (base : G) (a b : Int) :
    base ^ a * base ^ b = base ^ (a + b) :=
  (zpow_add base a b).symm

/-- Literal group semantics of `raise_to_z_div_by_2`. -/
def cppGroupRaiseToZDivBy2 (a : G) : G :=
  let e0 := cppRepeatedSquare 1 a
  let e1 := cppGroupMulNSqr e0 a 2
  let e2 := cppGroupMulNSqr e1 a 3
  let e3 := cppGroupMulNSqr e2 a 9
  let e4 := cppGroupMulNSqr e3 a 32
  let e5 := cppGroupMulNSqr e4 a (16 - 1)
  Inv.inv e5

/-- Literal group semantics of `raise_to_z`. -/
def cppGroupRaiseToZ (a : G) : G :=
  cppRepeatedSquare 1 (cppGroupRaiseToZDivBy2 a)

theorem cpp_group_raise_to_z_div_by_2_zpow (base : G) (a : Int) :
    cppGroupRaiseToZDivBy2 (base ^ a) =
      base ^ (raiseToZDivBy2Exponent a) := by
  simp only [cppGroupRaiseToZDivBy2]
  rw [cpp_repeated_square_zpow, cpp_group_mul_n_sqr_zpow,
    cpp_group_mul_n_sqr_zpow, cpp_group_mul_n_sqr_zpow,
    cpp_group_mul_n_sqr_zpow, cpp_group_mul_n_sqr_zpow,
    inverse_zpow_base]
  congr 1
  simp [raiseToZDivBy2Exponent, mulNSqrExponent]
  ring

theorem cpp_group_raise_to_z_zpow (base : G) (a : Int) :
    cppGroupRaiseToZ (base ^ a) = base ^ (raiseToZExponent a) := by
  rw [cppGroupRaiseToZ, cpp_group_raise_to_z_div_by_2_zpow,
    cpp_repeated_square_zpow]
  congr 1
  simp [raiseToZExponent]
  ring

/-- Literal hard part of `final_exp_after_inverse`, starting at `ret2`. -/
def cppGroupFinalHardSchedule
    (frob1 frob2 frob3 : G -> G) (ret2 : G) : G :=
  let ret2Sqr := cppRepeatedSquare 1 ret2
  let y0ToZ := cppGroupRaiseToZ ret2Sqr
  let y0ToZDiv2 := cppGroupRaiseToZDivBy2 y0ToZ
  let ret2Conj := Inv.inv ret2
  let y1MulY3 := ret2Conj * y0ToZ
  let y1MulY3Conj := Inv.inv y1MulY3
  let y1Final := y1MulY3Conj * y0ToZDiv2
  let y1ToZ := cppGroupRaiseToZ y1Final
  let y2ToZ := cppGroupRaiseToZ y1ToZ
  let y1FinalConj := Inv.inv y1Final
  let y3MulY1 := y2ToZ * y1FinalConj
  let y1Frob3 := frob3 y1Final
  let y1ToZFrob2 := frob2 y1ToZ
  let y1MulY2 := y1Frob3 * y1ToZFrob2
  let y3ToZ := cppGroupRaiseToZ y3MulY1
  let y2MulY0 := y3ToZ * ret2Sqr
  let y2MulRet := y2MulY0 * ret2
  let y1MulY2Final := y1MulY2 * y2MulRet
  let y3Frob1 := frob1 y3MulY1
  y1MulY2Final * y3Frob1

def exponentFrobenius (power : Nat) (value : G) : G :=
  value ^ (Generated.target_modulus ^ power)

theorem exponent_frobenius_zpow
    (base : G) (a : Int) (power : Nat) :
    exponentFrobenius power (base ^ a) =
      base ^ (a * (Generated.target_modulus ^ power : Nat)) := by
  rw [exponentFrobenius]
  rw [show (base ^ a) ^ (Generated.target_modulus ^ power) =
      (base ^ a) ^ ((Generated.target_modulus ^ power : Nat) : Int) from
        (zpow_natCast (base ^ a)
          (Generated.target_modulus ^ power)).symm]
  rw [(zpow_mul base a
    ((Generated.target_modulus ^ power : Nat) : Int)).symm]

/--
Literal easy part of `final_exp_after_inverse`.  The inversion provider's
result is represented by group inversion; proving the provider itself is a
separate base-field/inversion boundary.
-/
def cppGroupFinalEasySchedule (f : G) : G :=
  let fConj := exponentFrobenius 6 f
  let fInverse := Inv.inv f
  let ret1 := fConj * fInverse
  let frob2Ret1 := exponentFrobenius 2 ret1
  frob2Ret1 * ret1

theorem cpp_group_final_easy_schedule_zpow (base : G) (f : Int) :
    cppGroupFinalEasySchedule (base ^ f) =
      base ^ (f * finalEasyExponent) := by
  simp only [cppGroupFinalEasySchedule, exponent_frobenius_zpow,
    inverse_zpow_base, multiply_zpow_base]
  congr 1
  simp [finalEasyExponent]
  ring

/-- Complete literal group-level dataflow of `final_exp_after_inverse`. -/
def cppGroupFinalExponentiation (f : G) : G :=
  cppGroupFinalHardSchedule
    (exponentFrobenius 1) (exponentFrobenius 2)
    (exponentFrobenius 3) (cppGroupFinalEasySchedule f)

theorem cpp_group_final_hard_schedule_zpow (base : G) (ret2 : Int) :
    cppGroupFinalHardSchedule
        (exponentFrobenius 1) (exponentFrobenius 2)
        (exponentFrobenius 3) (base ^ ret2) =
      base ^ (finalHardScheduleExponent ret2) := by
  simp only [cppGroupFinalHardSchedule, cpp_repeated_square_zpow,
    cpp_group_raise_to_z_zpow, cpp_group_raise_to_z_div_by_2_zpow,
    inverse_zpow_base, multiply_zpow_base, exponent_frobenius_zpow]
  congr 1
  simp [finalHardScheduleExponent]
  simp_rw [raise_to_z_exponent_correct,
    raise_to_z_div_by_2_exponent_correct]
  ring

theorem cpp_group_final_exponentiation_is_cubed
    (base : G) :
    cppGroupFinalExponentiation base =
      base ^ cubedFinalExponent := by
  rw [cppGroupFinalExponentiation,
    show base = base ^ (1 : Int) by simp,
    cpp_group_final_easy_schedule_zpow,
    cpp_group_final_hard_schedule_zpow]
  simpa [cppFinalExponentiationExponent] using
    congrArg (fun exponent : Int => base ^ exponent)
      cpp_final_exponentiation_exponent_correct

end CatBLST.AVX.MatrixNoK
