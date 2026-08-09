/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.GtScalarSplit
import Mathlib.Tactic

set_option maxHeartbeats 800000

/-!
# Literal GT exponentiation schedules

This file models `avx_cyclotomic_exp_glv4_ct`, its unsigned variable-time
counterpart, and the binary fallback in the exponent group `ZMod r`.  Group
multiplication is addition of exponents, cyclotomic squaring doubles an
exponent, Frobenius supplies the four proved eigenvalue bases, and conjugation
negates an exponent.  This domain is GT; none of the statements extend the
contract to arbitrary Fp12 values.
-/

namespace CatBLST.AVX.MatrixNoK

/-- `(digit >> bit) & 1`, expressed with natural-number division. -/
def gtBitValue (digit bit : Nat) : Nat := digit / (2 ^ bit) % 2

theorem gt_bit_value_lt_two (digit bit : Nat) : gtBitValue digit bit < 2 := by
  exact Nat.mod_lt _ (by norm_num)

theorem gt_mod_two_pow_succ (digit bit : Nat) :
    digit % (2 ^ (bit + 1)) =
      digit % (2 ^ bit) + 2 ^ bit * gtBitValue digit bit := by
  let width := 2 ^ bit
  have widthPositive : 0 < width := pow_pos (by norm_num) bit
  have lowPart : digit % (width * 2) % width = digit % width := by
    exact Nat.mod_mul_right_mod digit width 2
  have highBit : digit % (width * 2) / width = digit / width % 2 := by
    exact Nat.mod_mul_right_div_self digit width 2
  have split := Nat.mod_add_div (digit % (width * 2)) width
  rw [lowPart, highBit] at split
  simpa [width, gtBitValue, pow_succ, Nat.mul_comm] using split.symm

/-- The disjoint-bit OR used by the source, written as its equivalent sum. -/
def gtJointIndex (digits : Fin 4 -> Nat) (bit : Nat) : Fin 16 :=
  let b0 := gtBitValue (digits 0) bit
  let b1 := gtBitValue (digits 1) bit
  let b2 := gtBitValue (digits 2) bit
  let b3 := gtBitValue (digits 3) bit
  have h0 : b0 < 2 := gt_bit_value_lt_two _ _
  have h1 : b1 < 2 := gt_bit_value_lt_two _ _
  have h2 : b2 < 2 := gt_bit_value_lt_two _ _
  have h3 : b3 < 2 := gt_bit_value_lt_two _ _
  ⟨b0 + 2 * b1 + 4 * b2 + 8 * b3, by omega⟩

/-- Exact source construction order for the sixteen-entry GLV4 table. -/
def gtSourceTableExponent (basis : Fin 4 -> GtExponent) : Fin 16 -> GtExponent :=
  ![0,
    basis 0,
    basis 1,
    basis 0 + basis 1,
    basis 2,
    basis 0 + basis 2,
    basis 1 + basis 2,
    basis 0 + basis 1 + basis 2,
    basis 3,
    basis 0 + basis 3,
    basis 1 + basis 3,
    basis 0 + basis 1 + basis 3,
    basis 2 + basis 3,
    basis 0 + basis 2 + basis 3,
    basis 1 + basis 2 + basis 3,
    basis 0 + basis 1 + basis 2 + basis 3]

theorem gt_source_table_zero (basis : Fin 4 -> GtExponent) :
    gtSourceTableExponent basis 0 = 0 := by
  rfl

private theorem gt_source_table_bits
    (basis : Fin 4 -> GtExponent) (b0 b1 b2 b3 : Nat)
    (h0 : b0 < 2) (h1 : b1 < 2) (h2 : b2 < 2) (h3 : b3 < 2) :
    gtSourceTableExponent basis
        ⟨b0 + 2 * b1 + 4 * b2 + 8 * b3, by omega⟩ =
      (b0 : GtExponent) * basis 0 +
      (b1 : GtExponent) * basis 1 +
      (b2 : GtExponent) * basis 2 +
      (b3 : GtExponent) * basis 3 := by
  have hb0 : b0 = 0 ∨ b0 = 1 := by omega
  have hb1 : b1 = 0 ∨ b1 = 1 := by omega
  have hb2 : b2 = 0 ∨ b2 = 1 := by omega
  have hb3 : b3 = 0 ∨ b3 = 1 := by omega
  rcases hb0 with rfl | rfl <;>
    rcases hb1 with rfl | rfl <;>
    rcases hb2 with rfl | rfl <;>
    rcases hb3 with rfl | rfl <;>
    simp [gtSourceTableExponent]

theorem gt_source_table_joint_lookup
    (basis : Fin 4 -> GtExponent) (digits : Fin 4 -> Nat) (bit : Nat) :
    gtSourceTableExponent basis (gtJointIndex digits bit) =
      (gtBitValue (digits 0) bit : GtExponent) * basis 0 +
      (gtBitValue (digits 1) bit : GtExponent) * basis 1 +
      (gtBitValue (digits 2) bit : GtExponent) * basis 2 +
      (gtBitValue (digits 3) bit : GtExponent) * basis 3 := by
  apply gt_source_table_bits
  all_goals exact gt_bit_value_lt_two _ _

/-- Exact descending CT loop: square, full-table select, multiply. -/
def gtJointCtLoop
    (table : Fin 16 -> GtExponent) (digits : Fin 4 -> Nat) :
    Nat -> GtExponent -> GtExponent
  | 0, out => out
  | remaining + 1, out =>
      gtJointCtLoop table digits remaining
        (out + out + table (gtJointIndex digits remaining))

def gtJointExpected
    (basis : Fin 4 -> GtExponent) (digits : Fin 4 -> Nat) (bits : Nat) :
    GtExponent :=
  Finset.univ.sum (fun index : Fin 4 =>
    ((digits index % (2 ^ bits) : Nat) : GtExponent) * basis index)

theorem gt_joint_expected_succ
    (basis : Fin 4 -> GtExponent) (digits : Fin 4 -> Nat) (bits : Nat) :
    gtJointExpected basis digits (bits + 1) =
      (2 ^ bits : GtExponent) *
        gtSourceTableExponent basis (gtJointIndex digits bits) +
      gtJointExpected basis digits bits := by
  have h0 := gt_mod_two_pow_succ (digits 0) bits
  have h1 := gt_mod_two_pow_succ (digits 1) bits
  have h2 := gt_mod_two_pow_succ (digits 2) bits
  have h3 := gt_mod_two_pow_succ (digits 3) bits
  have c0 := congrArg (fun value : Nat => (value : GtExponent)) h0
  have c1 := congrArg (fun value : Nat => (value : GtExponent)) h1
  have c2 := congrArg (fun value : Nat => (value : GtExponent)) h2
  have c3 := congrArg (fun value : Nat => (value : GtExponent)) h3
  push_cast at c0 c1 c2 c3
  rw [gt_source_table_joint_lookup]
  simp only [gtJointExpected, Fin.sum_univ_four]
  rw [c0, c1, c2, c3]
  ring

theorem gt_joint_ct_loop_correct
    (basis : Fin 4 -> GtExponent) (digits : Fin 4 -> Nat)
    (bits : Nat) (out : GtExponent) :
    gtJointCtLoop (gtSourceTableExponent basis) digits bits out =
      (2 ^ bits : GtExponent) * out + gtJointExpected basis digits bits := by
  induction bits generalizing out with
  | zero => simp [gtJointCtLoop, gtJointExpected, Nat.mod_one]
  | succ bits ih =>
      rw [gtJointCtLoop, ih, gt_joint_expected_succ]
      ring

theorem gt_radix_lt_two_pow_64 : gtRadix < 2 ^ 64 := by
  native_decide

theorem gt_canonical_digits_fit_64 (scalar : Nat) (index : Fin 4) :
    gtRadixZSplit (gtCanonicalScalar scalar) index < 2 ^ 64 :=
  (gt_canonical_split_digit_bounds scalar index).trans gt_radix_lt_two_pow_64

/-- The exact 64-round CT GLV4 executor returns the reduced scalar exponent. -/
theorem gt_unsigned_glv4_ct_correct (scalar : Nat) :
    gtJointCtLoop
        (gtSourceTableExponent gtUnsignedSourceBasisExponent)
        (gtRadixZSplit (gtCanonicalScalar scalar)) 64 0 =
      (gtCanonicalScalar scalar : GtExponent) := by
  rw [gt_joint_ct_loop_correct]
  simp only [mul_zero, zero_add]
  have fit := gt_canonical_digits_fit_64 scalar
  simp only [gtJointExpected, Fin.sum_univ_four]
  rw [Nat.mod_eq_of_lt (fit 0), Nat.mod_eq_of_lt (fit 1),
    Nat.mod_eq_of_lt (fit 2), Nat.mod_eq_of_lt (fit 3)]
  simpa [gtUnsignedSplitExponent, Fin.sum_univ_four] using
    gt_unsigned_split_exponent_correct scalar

/-- Exact source state machine for the unsigned variable-time executor. -/
def gtJointVtLoop
    (table : Fin 16 -> GtExponent) (digits : Fin 4 -> Nat) :
    Nat -> Bool -> GtExponent -> GtExponent
  | 0, _, out => out
  | remaining + 1, started, out =>
      let index := gtJointIndex digits remaining
      let squared := if started then out + out else out
      if index.val = 0 then
        gtJointVtLoop table digits remaining started squared
      else if started then
        gtJointVtLoop table digits remaining true (squared + table index)
      else
        gtJointVtLoop table digits remaining true (table index)

theorem gt_joint_vt_loop_eq_ct
    (basis : Fin 4 -> GtExponent) (digits : Fin 4 -> Nat)
    (bits : Nat) (started : Bool) (out : GtExponent)
    (notStarted : started = false -> out = 0) :
    gtJointVtLoop (gtSourceTableExponent basis) digits bits started out =
      gtJointCtLoop (gtSourceTableExponent basis) digits bits out := by
  induction bits generalizing started out with
  | zero => rfl
  | succ bits ih =>
      rw [gtJointVtLoop, gtJointCtLoop]
      by_cases hs : started = true
      · simp only [hs, if_true]
        by_cases hi : (gtJointIndex digits bits).val = 0
        · simp only [hi, if_true]
          have tableZero :
              gtSourceTableExponent basis (gtJointIndex digits bits) = 0 := by
            have indexZero : gtJointIndex digits bits = 0 := Fin.ext hi
            rw [indexZero, gt_source_table_zero]
          simpa [tableZero] using ih true (out + out) (by simp)
        · simp only [hi, if_false]
          exact ih true
            (out + out + gtSourceTableExponent basis (gtJointIndex digits bits))
            (by simp)
      · have hsFalse : started = false := Bool.eq_false_of_not_eq_true hs
        have outZero := notStarted hsFalse
        subst started
        subst out
        simp only [Bool.false_eq_true, if_false, zero_add]
        by_cases hi : (gtJointIndex digits bits).val = 0
        · simp only [hi, if_true]
          have tableZero :
              gtSourceTableExponent basis (gtJointIndex digits bits) = 0 := by
            have indexZero : gtJointIndex digits bits = 0 := Fin.ext hi
            rw [indexZero, gt_source_table_zero]
          simpa [tableZero] using ih false 0 (by simp)
        · simp only [hi, if_false]
          exact ih true (gtSourceTableExponent basis (gtJointIndex digits bits))
            (by simp)

/-- Skipping leading zero and zero-window work preserves the CT exponent. -/
theorem gt_unsigned_glv4_vartime_correct (scalar : Nat) :
    gtJointVtLoop
        (gtSourceTableExponent gtUnsignedSourceBasisExponent)
        (gtRadixZSplit (gtCanonicalScalar scalar)) 64 false 0 =
      (gtCanonicalScalar scalar : GtExponent) := by
  rw [gt_joint_vt_loop_eq_ct _ _ _ false 0 (by simp)]
  exact gt_unsigned_glv4_ct_correct scalar

/-- Exact binary CT fallback: square and select either squared or product. -/
def gtBinaryCtLoop (scalar : Nat) : Nat -> GtExponent -> GtExponent
  | 0, out => out
  | remaining + 1, out =>
      gtBinaryCtLoop scalar remaining
        (out + out + (gtBitValue scalar remaining : GtExponent))

theorem gt_binary_ct_loop_correct (scalar bits : Nat) (out : GtExponent) :
    gtBinaryCtLoop scalar bits out =
      (2 ^ bits : GtExponent) * out +
        ((scalar % (2 ^ bits) : Nat) : GtExponent) := by
  induction bits generalizing out with
  | zero => simp [gtBinaryCtLoop, Nat.mod_one]
  | succ bits ih =>
      rw [gtBinaryCtLoop, ih]
      have recurrence := gt_mod_two_pow_succ scalar bits
      have castRecurrence := congrArg
        (fun value : Nat => (value : GtExponent)) recurrence
      push_cast at castRecurrence
      rw [castRecurrence]
      ring

/-- For public `bits`, the binary fallback computes the masked scalar. -/
theorem gt_binary_ct_correct (scalar bits : Nat) :
    gtBinaryCtLoop scalar bits 0 =
      ((scalar % (2 ^ bits) : Nat) : GtExponent) := by
  rw [gt_binary_ct_loop_correct]
  ring

/-- Exact basic variable-time source loop, including leading-zero skipping. -/
def gtBinaryVtLoop (scalar : Nat) : Nat -> Bool -> GtExponent -> GtExponent
  | 0, _, out => out
  | remaining + 1, started, out =>
      let bit := gtBitValue scalar remaining
      let squared := if started then out + out else out
      if bit = 0 then
        gtBinaryVtLoop scalar remaining started squared
      else if started then
        gtBinaryVtLoop scalar remaining true (squared + 1)
      else
        gtBinaryVtLoop scalar remaining true 1

theorem gt_binary_vt_loop_eq_ct
    (scalar bits : Nat) (started : Bool) (out : GtExponent)
    (notStarted : started = false -> out = 0) :
    gtBinaryVtLoop scalar bits started out = gtBinaryCtLoop scalar bits out := by
  induction bits generalizing started out with
  | zero => rfl
  | succ bits ih =>
      rw [gtBinaryVtLoop, gtBinaryCtLoop]
      have bitBound := gt_bit_value_lt_two scalar bits
      have bitCases : gtBitValue scalar bits = 0 ∨ gtBitValue scalar bits = 1 := by
        omega
      rcases bitCases with bitZero | bitOne
      · simp only [bitZero, if_true, Nat.cast_zero]
        by_cases hs : started = true
        · simp only [hs, if_true]
          simpa using ih true (out + out) (by simp)
        · have hsFalse : started = false := Bool.eq_false_of_not_eq_true hs
          have outZero := notStarted hsFalse
          subst started
          subst out
          simp only [Bool.false_eq_true, if_false, zero_add]
          simpa using ih false 0 (by simp)
      · simp only [bitOne, Nat.cast_one]
        by_cases hs : started = true
        · simp only [hs, if_true]
          simpa using ih true (out + out + 1) (by simp)
        · have hsFalse : started = false := Bool.eq_false_of_not_eq_true hs
          have outZero := notStarted hsFalse
          subst started
          subst out
          simp only [Bool.false_eq_true, if_false, zero_add]
          simpa using ih true 1 (by simp)

/-- The source's fixed 256-bit basic VT scan computes the supplied limbs. -/
theorem gt_binary_vartime_256_correct (scalar : Nat)
    (scalarBound : scalar < 2 ^ 256) :
    gtBinaryVtLoop scalar 256 false 0 = (scalar : GtExponent) := by
  rw [gt_binary_vt_loop_eq_ct _ _ false 0 (by simp), gt_binary_ct_correct,
    Nat.mod_eq_of_lt scalarBound]

/-- Public GT CT dispatch after the C++ byte-to-limb masking step. -/
def gtPublicCtExponent (scalar bits : Nat) : GtExponent :=
  let masked := scalar % (2 ^ bits)
  if bits > 128 ∧ bits <= 256 then
    gtJointCtLoop
      (gtSourceTableExponent gtUnsignedSourceBasisExponent)
      (gtRadixZSplit (gtCanonicalScalar masked)) 64 0
  else
    gtBinaryCtLoop masked bits 0

theorem gt_public_ct_exponent_correct
    (scalar bits : Nat) (_bitsBound : bits <= 256) :
    gtPublicCtExponent scalar bits =
      (gtCanonicalScalar (scalar % (2 ^ bits)) : GtExponent) := by
  simp only [gtPublicCtExponent]
  by_cases useGlv : bits > 128 ∧ bits <= 256
  · simp only [useGlv]
    exact gt_unsigned_glv4_ct_correct (scalar % (2 ^ bits))
  · simp only [useGlv, if_false]
    rw [gt_binary_ct_correct, Nat.mod_eq_of_lt (Nat.mod_lt scalar (pow_pos (by norm_num) bits))]
    change ((scalar % 2 ^ bits : Nat) : GtExponent) =
      ((gtCanonicalScalar (scalar % 2 ^ bits) : Nat) : GtExponent)
    simp only [gtCanonicalScalar]
    exact (ZMod.natCast_mod (scalar % 2 ^ bits) blsScalarOrder).symm

end CatBLST.AVX.MatrixNoK
