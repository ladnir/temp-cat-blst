/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import Mathlib.Data.Fin.Basic
import Mathlib.Data.Nat.Bitwise

/-!
# AVX-512 lane semantics

This file gives the mathematical model used by the CatBLST AVX proof. A
machine lane is an unsigned integer modulo `2^64`; a vector is exactly eight
such lanes. The definitions are specifications, not executable claims about a
particular compiler or processor.
-/

namespace CatBLST.AVX

def wordBits : Nat := 64
def digitBits : Nat := 52
def laneCount : Nat := 8

def wordModulus : Nat := 2 ^ wordBits
def digitModulus : Nat := 2 ^ digitBits

abbrev Word64 := Fin wordModulus
abbrev Word52 := Fin digitModulus
abbrev LaneIndex := Fin laneCount
abbrev Vec8 := LaneIndex → Word64

namespace Word64

/-- Reduce a natural number to the unsigned 64-bit machine-word domain. -/
def ofNat (value : Nat) : Word64 :=
  ⟨value % wordModulus, Nat.mod_lt _ (by decide)⟩

@[simp] theorem val_ofNat (value : Nat) :
    (ofNat value).val = value % wordModulus := rfl

/-- Unsigned wrapping addition. -/
def add (a b : Word64) : Word64 := ofNat (a.val + b.val)

/-- Unsigned wrapping subtraction. -/
def sub (a b : Word64) : Word64 := ofNat (a.val + wordModulus - b.val)

/-- Logical right shift of an unsigned lane. -/
def shiftRight (a : Word64) (bits : Nat) : Word64 := ofNat (a.val >>> bits)

/-- Logical left shift with the high bits discarded. -/
def shiftLeft (a : Word64) (bits : Nat) : Word64 := ofNat (a.val <<< bits)

/-- Bitwise conjunction in the 64-bit lane domain. -/
def bitAnd (a b : Word64) : Word64 := ofNat (a.val &&& b.val)

theorem add_val_of_no_wrap {a b : Word64}
    (h : a.val + b.val < wordModulus) :
    (add a b).val = a.val + b.val := by
  simp [add, val_ofNat, Nat.mod_eq_of_lt h]

end Word64

namespace Vec8

def set1 (value : Word64) : Vec8 := fun _ => value

def add (a b : Vec8) : Vec8 := fun lane => Word64.add (a lane) (b lane)

def sub (a b : Vec8) : Vec8 := fun lane => Word64.sub (a lane) (b lane)

def shiftRight (a : Vec8) (bits : Nat) : Vec8 :=
  fun lane => Word64.shiftRight (a lane) bits

def shiftLeft (a : Vec8) (bits : Nat) : Vec8 :=
  fun lane => Word64.shiftLeft (a lane) bits

def bitAnd (a b : Vec8) : Vec8 := fun lane => Word64.bitAnd (a lane) (b lane)

/--
AVX-512 variable lane permutation. Only the low three bits of each index lane
select the source, exactly as `_mm512_permutexvar_epi64` does.
-/
def permute (indices values : Vec8) : Vec8 := fun lane =>
  values ⟨(indices lane).val % laneCount,
    Nat.mod_lt _ (by decide)⟩

def lowLane (values : Vec8) : Word64 := values ⟨0, by decide⟩

def broadcastLow (values : Vec8) : Vec8 := set1 (lowLane values)

@[simp] theorem set1_apply (value : Word64) (lane : LaneIndex) :
    set1 value lane = value := rfl

@[simp] theorem broadcastLow_apply (values : Vec8) (lane : LaneIndex) :
    broadcastLow values lane = lowLane values := rfl

theorem lowLaneExtractionAndBroadcast (values : Vec8) :
    broadcastLow values = fun _ => values ⟨0, by decide⟩ := rfl

end Vec8

end CatBLST.AVX
