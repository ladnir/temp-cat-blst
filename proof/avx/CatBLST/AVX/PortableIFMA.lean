/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.IFMA

/-!
# Portable IFMA decomposition

The portable CatBLST implementation computes a 52-by-52-bit product from
26-bit halves so that every C++ multiplication fits in `uint64_t`. These
definitions mirror that source schedule. The central theorem proves that its
low and unmasked-high outputs reconstruct the original product exactly.
-/

namespace CatBLST.AVX.PortableIFMA

def halfDigitBits : Nat := 26
def halfDigitModulus : Nat := 2 ^ halfDigitBits
def halfDigitMask : Nat := halfDigitModulus - 1

def halfLow (value : Nat) : Nat := value % halfDigitModulus
def halfHigh (value : Nat) : Nat := value / halfDigitModulus

def p0 (x y : Nat) : Nat := halfLow x * halfLow y
def p1 (x y : Nat) : Nat := halfLow x * halfHigh y + halfHigh x * halfLow y
def p2 (x y : Nat) : Nat := halfHigh x * halfHigh y

def carry0 (x y : Nat) : Nat := p0 x y / halfDigitModulus
def coefficient0 (x y : Nat) : Nat := p0 x y % halfDigitModulus

def coefficient1Wide (x y : Nat) : Nat :=
  carry0 x y + p1 x y % halfDigitModulus

def coefficient1 (x y : Nat) : Nat :=
  coefficient1Wide x y % halfDigitModulus

def carry1 (x y : Nat) : Nat :=
  coefficient1Wide x y / halfDigitModulus

def coefficient2Wide (x y : Nat) : Nat :=
  p1 x y / halfDigitModulus + p2 x y % halfDigitModulus + carry1 x y

def coefficient2 (x y : Nat) : Nat :=
  coefficient2Wide x y % halfDigitModulus

def carry2 (x y : Nat) : Nat :=
  coefficient2Wide x y / halfDigitModulus

def coefficient3 (x y : Nat) : Nat :=
  p2 x y / halfDigitModulus + carry2 x y

def low (x y : Nat) : Nat :=
  coefficient0 x y + halfDigitModulus * coefficient1 x y

/-- High half before the final defensive 52-bit source mask. -/
def highRaw (x y : Nat) : Nat :=
  coefficient2 x y + halfDigitModulus * coefficient3 x y

/-- Exact source expression, including the final mask on coefficient 3. -/
def high (x y : Nat) : Nat :=
  coefficient2 x y +
    halfDigitModulus * (coefficient3 x y % halfDigitModulus)

/-
The following two expressions retain the source's shifts, masks, and bitwise
OR. They are kept separate from `low`/`high` so the algebraic proof is readable
while the correspondence boundary remains exact.
-/
def sourceLow (x y : Nat) : Nat :=
  coefficient0 x y ||| (coefficient1 x y <<< halfDigitBits)

def sourceHigh (x y : Nat) : Nat :=
  coefficient2 x y |||
    ((coefficient3 x y &&& halfDigitMask) <<< halfDigitBits)

theorem digitModulus_eq :
    digitModulus = halfDigitModulus * halfDigitModulus := by
  norm_num [digitModulus, digitBits, halfDigitModulus, halfDigitBits]

theorem half_split (value : Nat) :
    halfLow value + halfDigitModulus * halfHigh value = value := by
  exact Nat.mod_add_div value halfDigitModulus

theorem halfLow_eq_source_mask (value : Nat) :
    halfLow value = value &&& halfDigitMask := by
  simp [halfLow, halfDigitMask, halfDigitModulus]

theorem halfHigh_eq_source_shift (value : Nat) :
    halfHigh value = value >>> halfDigitBits := by
  simp [halfHigh, halfDigitModulus, Nat.shiftRight_eq_div_pow]

theorem sourceLow_eq_low (x y : Nat) : sourceLow x y = low x y := by
  have hc0Base : coefficient0 x y < halfDigitModulus :=
    Nat.mod_lt _ (by norm_num [halfDigitModulus, halfDigitBits])
  have hc0 : coefficient0 x y < 2 ^ halfDigitBits := by
    simpa [halfDigitModulus] using hc0Base
  have hor := Nat.shiftLeft_add_eq_or_of_lt hc0 (coefficient1 x y)
  calc
    sourceLow x y =
        (coefficient1 x y <<< halfDigitBits) ||| coefficient0 x y := by
      simp only [sourceLow, Nat.or_comm]
    _ = (coefficient1 x y <<< halfDigitBits) + coefficient0 x y := hor.symm
    _ = low x y := by
      rw [Nat.shiftLeft_eq]
      simp only [low, halfDigitModulus]
      ring

theorem sourceHigh_eq_high (x y : Nat) : sourceHigh x y = high x y := by
  have hc2Base : coefficient2 x y < halfDigitModulus :=
    Nat.mod_lt _ (by norm_num [halfDigitModulus, halfDigitBits])
  have hc2 : coefficient2 x y < 2 ^ halfDigitBits := by
    simpa [halfDigitModulus] using hc2Base
  have hor := Nat.shiftLeft_add_eq_or_of_lt hc2
    (coefficient3 x y &&& halfDigitMask)
  calc
    sourceHigh x y =
        ((coefficient3 x y &&& halfDigitMask) <<< halfDigitBits) |||
          coefficient2 x y := by
      simp only [sourceHigh, Nat.or_comm]
    _ = ((coefficient3 x y &&& halfDigitMask) <<< halfDigitBits) +
          coefficient2 x y := hor.symm
    _ = high x y := by
      rw [Nat.shiftLeft_eq]
      simp only [high, halfDigitMask, halfDigitModulus,
        Nat.and_two_pow_sub_one_eq_mod]
      ring

theorem decomposition (x y : Nat) :
    low x y + digitModulus * highRaw x y = x * y := by
  have hp0 := Nat.mod_add_div (p0 x y) halfDigitModulus
  have hp1 := Nat.mod_add_div (p1 x y) halfDigitModulus
  have hp2 := Nat.mod_add_div (p2 x y) halfDigitModulus
  have hc1 := Nat.mod_add_div (coefficient1Wide x y) halfDigitModulus
  have hc2 := Nat.mod_add_div (coefficient2Wide x y) halfDigitModulus
  have hcoeff :
      low x y + digitModulus * highRaw x y =
        p0 x y + halfDigitModulus * p1 x y + digitModulus * p2 x y := by
    rw [digitModulus_eq]
    simp only [low, highRaw, coefficient0, coefficient1, coefficient1Wide,
      carry0, coefficient2, coefficient2Wide, carry1, coefficient3, carry2]
    norm_num [halfDigitModulus, halfDigitBits] at hp0 hp1 hp2 hc1 hc2 ⊢
    omega
  calc
    low x y + digitModulus * highRaw x y =
        p0 x y + halfDigitModulus * p1 x y + digitModulus * p2 x y := hcoeff
    _ = (halfLow x + halfDigitModulus * halfHigh x) *
          (halfLow y + halfDigitModulus * halfHigh y) := by
        rw [digitModulus_eq]
        simp only [p0, p1, p2]
        ring
    _ = x * y := by rw [half_split, half_split]

theorem low_lt_digitModulus (x y : Nat) : low x y < digitModulus := by
  have hc0 : coefficient0 x y < halfDigitModulus :=
    Nat.mod_lt _ (by norm_num [halfDigitModulus, halfDigitBits])
  have hc1 : coefficient1 x y < halfDigitModulus :=
    Nat.mod_lt _ (by norm_num [halfDigitModulus, halfDigitBits])
  rw [digitModulus_eq]
  unfold low
  nlinarith

theorem low_eq_product_mod (x y : Nat) :
    low x y = (x * y) % digitModulus := by
  rw [← decomposition x y]
  simp [Nat.mod_eq_of_lt (low_lt_digitModulus x y)]

theorem highRaw_eq_product_div (x y : Nat) :
    highRaw x y = (x * y) / digitModulus := by
  rw [← decomposition x y]
  rw [Nat.add_mul_div_left _ _ (by norm_num [digitModulus, digitBits])]
  simp [Nat.div_eq_of_lt (low_lt_digitModulus x y)]

theorem highRaw_lt_digitModulus {x y : Nat}
    (hx : x < digitModulus) (hy : y < digitModulus) :
    highRaw x y < digitModulus := by
  rw [highRaw_eq_product_div]
  apply (Nat.div_lt_iff_lt_mul (by norm_num [digitModulus, digitBits])).2
  nlinarith

theorem coefficient3_lt {x y : Nat}
    (hx : x < digitModulus) (hy : y < digitModulus) :
    coefficient3 x y < halfDigitModulus := by
  have hraw := highRaw_lt_digitModulus hx hy
  have hc2 : coefficient2 x y < halfDigitModulus :=
    Nat.mod_lt _ (by norm_num [halfDigitModulus, halfDigitBits])
  rw [digitModulus_eq] at hraw
  unfold highRaw at hraw
  nlinarith

theorem high_eq_highRaw {x y : Nat}
    (hx : x < digitModulus) (hy : y < digitModulus) :
    high x y = highRaw x y := by
  simp [high, highRaw, Nat.mod_eq_of_lt (coefficient3_lt hx hy)]

theorem halfLow_lt (value : Nat) : halfLow value < halfDigitModulus :=
  Nat.mod_lt _ (by norm_num [halfDigitModulus, halfDigitBits])

theorem halfHigh_lt {value : Nat} (hvalue : value < digitModulus) :
    halfHigh value < halfDigitModulus := by
  apply (Nat.div_lt_iff_lt_mul
    (by norm_num [halfDigitModulus, halfDigitBits])).2
  rw [← digitModulus_eq]
  exact hvalue

/--
Every non-accumulator temporary in the portable C++ schedule fits in its
`uint64_t` owner. The final accumulator addition is intentionally excluded:
its specified behavior is wrapping addition and is modeled in `IFMA.lean`.
-/
theorem schedule_temporaries_fit_word {x y : Nat}
    (hx : x < digitModulus) (hy : y < digitModulus) :
    p0 x y < wordModulus ∧
    p1 x y < wordModulus ∧
    p2 x y < wordModulus ∧
    coefficient1Wide x y < wordModulus ∧
    coefficient2Wide x y < wordModulus ∧
    sourceLow x y < wordModulus ∧
    sourceHigh x y < wordModulus := by
  have hx0 := halfLow_lt x
  have hx1 := halfHigh_lt hx
  have hy0 := halfLow_lt y
  have hy1 := halfHigh_lt hy
  have hp0Digit : p0 x y < digitModulus := by
    rw [digitModulus_eq]
    unfold p0
    nlinarith
  have hp1TwoDigit : p1 x y < 2 * digitModulus := by
    rw [digitModulus_eq]
    unfold p1
    nlinarith
  have hp2Digit : p2 x y < digitModulus := by
    rw [digitModulus_eq]
    unfold p2
    nlinarith
  have hcarry0 : carry0 x y < halfDigitModulus := by
    unfold carry0
    apply (Nat.div_lt_iff_lt_mul
      (by norm_num [halfDigitModulus, halfDigitBits])).2
    rw [← digitModulus_eq]
    exact hp0Digit
  have hp1Mod : p1 x y % halfDigitModulus < halfDigitModulus :=
    Nat.mod_lt _ (by norm_num [halfDigitModulus, halfDigitBits])
  have hc1Wide : coefficient1Wide x y < 2 * halfDigitModulus := by
    unfold coefficient1Wide
    omega
  have hp1Div : p1 x y / halfDigitModulus < 2 * halfDigitModulus := by
    apply (Nat.div_lt_iff_lt_mul
      (by norm_num [halfDigitModulus, halfDigitBits])).2
    rw [digitModulus_eq] at hp1TwoDigit
    simpa [Nat.mul_assoc] using hp1TwoDigit
  have hp2Mod : p2 x y % halfDigitModulus < halfDigitModulus :=
    Nat.mod_lt _ (by norm_num [halfDigitModulus, halfDigitBits])
  have hcarry1 : carry1 x y < 2 := by
    unfold carry1
    apply (Nat.div_lt_iff_lt_mul
      (by norm_num [halfDigitModulus, halfDigitBits])).2
    simpa [Nat.mul_comm] using hc1Wide
  have hc2Wide : coefficient2Wide x y <
      3 * halfDigitModulus + 2 := by
    unfold coefficient2Wide
    omega
  have hsourceLowDigit : sourceLow x y < digitModulus := by
    rw [sourceLow_eq_low]
    exact low_lt_digitModulus x y
  have hsourceHighDigit : sourceHigh x y < digitModulus := by
    rw [sourceHigh_eq_high, high_eq_highRaw hx hy]
    exact highRaw_lt_digitModulus hx hy
  have hdigitWord : digitModulus < wordModulus := by
    norm_num [digitModulus, digitBits, wordModulus, wordBits]
  have htwoDigitWord : 2 * digitModulus < wordModulus := by
    norm_num [digitModulus, digitBits, wordModulus, wordBits]
  have htwoHalfWord : 2 * halfDigitModulus < wordModulus := by
    norm_num [halfDigitModulus, halfDigitBits, wordModulus, wordBits]
  have hthreeHalfWord : 3 * halfDigitModulus + 2 < wordModulus := by
    norm_num [halfDigitModulus, halfDigitBits, wordModulus, wordBits]
  exact ⟨hp0Digit.trans hdigitWord,
    hp1TwoDigit.trans htwoDigitWord,
    hp2Digit.trans hdigitWord,
    hc1Wide.trans htwoHalfWord,
    hc2Wide.trans hthreeHalfWord,
    hsourceLowDigit.trans hdigitWord,
    hsourceHighDigit.trans hdigitWord⟩

theorem low_matches_ifma (a b : Word64) :
    sourceLow (low52 a).val (low52 b).val = productLow52 a b := by
  rw [sourceLow_eq_low]
  exact low_eq_product_mod _ _

theorem high_matches_ifma (a b : Word64) :
    sourceHigh (low52 a).val (low52 b).val = productHigh52 a b := by
  rw [sourceHigh_eq_high]
  rw [high_eq_highRaw (low52 a).isLt (low52 b).isLt]
  exact highRaw_eq_product_div _ _

end CatBLST.AVX.PortableIFMA
