/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.Lane
import Mathlib.Tactic

/-!
# AVX-512 IFMA semantics

`madd52lo` and `madd52hi` multiply the low 52 bits of two unsigned lanes,
select the low or high radix-`2^52` digit of the 104-bit product, and add that
digit to a 64-bit accumulator with machine-word wraparound.

The definitions below are the P0 mathematical specification. The theorems
prove the exact product split, digit bounds, wrapping semantics, and the
no-wrap interpretation required by later MatrixNoK bound proofs.
-/

namespace CatBLST.AVX

/-- The low 52 bits of a 64-bit lane. -/
def low52 (value : Word64) : Word52 :=
  ⟨value.val % digitModulus,
    Nat.mod_lt _ (by norm_num [digitModulus, digitBits])⟩

/-- Product of the two low-52-bit input digits, before selecting a half. -/
def product52 (a b : Word64) : Nat := (low52 a).val * (low52 b).val

def productLow52 (a b : Word64) : Nat := product52 a b % digitModulus
def productHigh52 (a b : Word64) : Nat := product52 a b / digitModulus

/-- Lane semantics of `vpmadd52luq`. -/
def madd52loLane (acc a b : Word64) : Word64 :=
  Word64.ofNat (acc.val + productLow52 a b)

/-- Lane semantics of `vpmadd52huq`. -/
def madd52hiLane (acc a b : Word64) : Word64 :=
  Word64.ofNat (acc.val + productHigh52 a b)

theorem product52_split (a b : Word64) :
    productLow52 a b + digitModulus * productHigh52 a b = product52 a b := by
  exact Nat.mod_add_div (product52 a b) digitModulus

theorem productLow52_lt (a b : Word64) :
    productLow52 a b < digitModulus := by
  exact Nat.mod_lt _ (by norm_num [digitModulus, digitBits])

theorem product52_lt_square (a b : Word64) :
    product52 a b < digitModulus * digitModulus := by
  have ha := (low52 a).isLt
  have hb := (low52 b).isLt
  simp only [product52]
  nlinarith

theorem productHigh52_lt (a b : Word64) :
    productHigh52 a b < digitModulus := by
  apply (Nat.div_lt_iff_lt_mul (by norm_num [digitModulus, digitBits])).2
  simpa [Nat.mul_comm] using product52_lt_square a b

@[simp] theorem madd52loLane_val (acc a b : Word64) :
    (madd52loLane acc a b).val =
      (acc.val + productLow52 a b) % wordModulus := rfl

@[simp] theorem madd52hiLane_val (acc a b : Word64) :
    (madd52hiLane acc a b).val =
      (acc.val + productHigh52 a b) % wordModulus := rfl

theorem madd52loLane_val_of_no_wrap {acc a b : Word64}
    (h : acc.val + productLow52 a b < wordModulus) :
    (madd52loLane acc a b).val = acc.val + productLow52 a b := by
  simp [madd52loLane, Word64.val_ofNat, Nat.mod_eq_of_lt h]

theorem madd52hiLane_val_of_no_wrap {acc a b : Word64}
    (h : acc.val + productHigh52 a b < wordModulus) :
    (madd52hiLane acc a b).val = acc.val + productHigh52 a b := by
  simp [madd52hiLane, Word64.val_ofNat, Nat.mod_eq_of_lt h]

abbrev ProductPair := Word64 × Word64

def lowContribution (term : ProductPair) : Nat :=
  productLow52 term.1 term.2

def highContribution (term : ProductPair) : Nat :=
  productHigh52 term.1 term.2

/-- Repeated low-half IFMA accumulation in one lane. -/
def accumulateLow : Word64 → List ProductPair → Word64
  | acc, [] => acc
  | acc, term :: rest =>
      accumulateLow (madd52loLane acc term.1 term.2) rest

/-- Repeated high-half IFMA accumulation in one lane. -/
def accumulateHigh : Word64 → List ProductPair → Word64
  | acc, [] => acc
  | acc, term :: rest =>
      accumulateHigh (madd52hiLane acc term.1 term.2) rest

/-
These two theorems are the reusable no-unintended-wrap rule. Later kernels
must prove the stated bound from their concrete residue and matrix bounds; P0
then removes all machine-word modular reductions from the accumulated result.
-/
theorem accumulateLow_val_of_no_wrap (acc : Word64) (terms : List ProductPair)
    (h : acc.val + (terms.map lowContribution).sum < wordModulus) :
    (accumulateLow acc terms).val =
      acc.val + (terms.map lowContribution).sum := by
  induction terms generalizing acc with
  | nil => simp [accumulateLow]
  | cons term rest ih =>
      have hstep : acc.val + lowContribution term < wordModulus := by
        simp only [List.map_cons, List.sum_cons] at h
        omega
      have hval :
          (madd52loLane acc term.1 term.2).val =
            acc.val + lowContribution term := by
        exact madd52loLane_val_of_no_wrap hstep
      have hrest :
          (madd52loLane acc term.1 term.2).val +
              (rest.map lowContribution).sum < wordModulus := by
        rw [hval]
        simp only [List.map_cons, List.sum_cons] at h
        omega
      have hind := ih (madd52loLane acc term.1 term.2) hrest
      rw [accumulateLow, hind, hval]
      simp only [List.map_cons, List.sum_cons]
      omega

theorem accumulateHigh_val_of_no_wrap (acc : Word64) (terms : List ProductPair)
    (h : acc.val + (terms.map highContribution).sum < wordModulus) :
    (accumulateHigh acc terms).val =
      acc.val + (terms.map highContribution).sum := by
  induction terms generalizing acc with
  | nil => simp [accumulateHigh]
  | cons term rest ih =>
      have hstep : acc.val + highContribution term < wordModulus := by
        simp only [List.map_cons, List.sum_cons] at h
        omega
      have hval :
          (madd52hiLane acc term.1 term.2).val =
            acc.val + highContribution term := by
        exact madd52hiLane_val_of_no_wrap hstep
      have hrest :
          (madd52hiLane acc term.1 term.2).val +
              (rest.map highContribution).sum < wordModulus := by
        rw [hval]
        simp only [List.map_cons, List.sum_cons] at h
        omega
      have hind := ih (madd52hiLane acc term.1 term.2) hrest
      rw [accumulateHigh, hind, hval]
      simp only [List.map_cons, List.sum_cons]
      omega

namespace Vec8

def madd52lo (acc a b : Vec8) : Vec8 :=
  fun lane => madd52loLane (acc lane) (a lane) (b lane)

def madd52hi (acc a b : Vec8) : Vec8 :=
  fun lane => madd52hiLane (acc lane) (a lane) (b lane)

@[simp] theorem madd52lo_apply (acc a b : Vec8) (lane : LaneIndex) :
    madd52lo acc a b lane = madd52loLane (acc lane) (a lane) (b lane) := rfl

@[simp] theorem madd52hi_apply (acc a b : Vec8) (lane : LaneIndex) :
    madd52hi acc a b lane = madd52hiLane (acc lane) (a lane) (b lane) := rfl

theorem madd52lo_val_of_no_wrap {acc a b : Vec8}
    (h : ∀ lane, (acc lane).val + productLow52 (a lane) (b lane) < wordModulus) :
    ∀ lane, (madd52lo acc a b lane).val =
      (acc lane).val + productLow52 (a lane) (b lane) := by
  intro lane
  exact madd52loLane_val_of_no_wrap (h lane)

theorem madd52hi_val_of_no_wrap {acc a b : Vec8}
    (h : ∀ lane, (acc lane).val + productHigh52 (a lane) (b lane) < wordModulus) :
    ∀ lane, (madd52hi acc a b lane).val =
      (acc lane).val + productHigh52 (a lane) (b lane) := by
  intro lane
  exact madd52hiLane_val_of_no_wrap (h lane)

end Vec8

end CatBLST.AVX
