/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Fp2
import Mathlib.Tactic

/-!
# Concrete selected Fp2 algebra

This packages the already accepted `Fp2Value` formulas as a commutative ring.
It lets the P4 curve and Miller proofs use ordinary polynomial reasoning while
retaining the exact `u^2 = -1` representation proved in P3-001.
-/

namespace CatBLST.AVX.MatrixNoK

namespace Fp2Concrete

variable {R : Type*}

def zero [Zero R] : Fp2Value R := ⟨0, 0⟩
def one [Zero R] [One R] : Fp2Value R := ⟨1, 0⟩
def add [Add R] (x y : Fp2Value R) : Fp2Value R := Fp2Value.add x y
def neg [Neg R] (x : Fp2Value R) : Fp2Value R := Fp2Value.neg x

instance fp2Zero [Zero R] : Zero (Fp2Value R) := ⟨zero⟩
instance fp2One [Zero R] [One R] : One (Fp2Value R) := ⟨one⟩
instance fp2Add [Add R] : Add (Fp2Value R) := ⟨add⟩
instance fp2Neg [Neg R] : Neg (Fp2Value R) := ⟨neg⟩
instance fp2Mul [Ring R] : Mul (Fp2Value R) := ⟨Fp2Value.mul⟩

private theorem add_assoc_fn [CommRing R] (x y z : Fp2Value R) :
    add (add x y) z = add x (add y z) := by
  ext <;> simp [add, Fp2Value.add] <;> ring

private theorem zero_add_fn [CommRing R] (x : Fp2Value R) :
    add zero x = x := by
  ext <;> simp [zero, add, Fp2Value.add]

private theorem add_zero_fn [CommRing R] (x : Fp2Value R) :
    add x zero = x := by
  ext <;> simp [zero, add, Fp2Value.add]

private theorem add_comm_fn [CommRing R] (x y : Fp2Value R) :
    add x y = add y x := by
  ext <;> simp [add, Fp2Value.add] <;> ring

private theorem neg_add_cancel_fn [CommRing R] (x : Fp2Value R) :
    add (neg x) x = zero := by
  ext <;> simp [neg, add, zero, Fp2Value.neg, Fp2Value.add]

private theorem mul_assoc_fn [CommRing R] (x y z : Fp2Value R) :
    Fp2Value.mul (Fp2Value.mul x y) z =
      Fp2Value.mul x (Fp2Value.mul y z) := by
  ext <;> simp [Fp2Value.mul] <;> ring

private theorem one_mul_fn [CommRing R] (x : Fp2Value R) :
    Fp2Value.mul one x = x := by
  ext <;> simp [one, Fp2Value.mul]

private theorem mul_one_fn [CommRing R] (x : Fp2Value R) :
    Fp2Value.mul x one = x := by
  ext <;> simp [one, Fp2Value.mul]

private theorem zero_mul_fn [CommRing R] (x : Fp2Value R) :
    Fp2Value.mul zero x = zero := by
  ext <;> simp [zero, Fp2Value.mul]

private theorem mul_zero_fn [CommRing R] (x : Fp2Value R) :
    Fp2Value.mul x zero = zero := by
  ext <;> simp [zero, Fp2Value.mul]

private theorem left_distrib_fn [CommRing R] (x y z : Fp2Value R) :
    Fp2Value.mul x (add y z) =
      add (Fp2Value.mul x y) (Fp2Value.mul x z) := by
  ext <;> simp [add, Fp2Value.add, Fp2Value.mul] <;> ring

private theorem right_distrib_fn [CommRing R] (x y z : Fp2Value R) :
    Fp2Value.mul (add x y) z =
      add (Fp2Value.mul x z) (Fp2Value.mul y z) := by
  ext <;> simp [add, Fp2Value.add, Fp2Value.mul] <;> ring

private theorem mul_comm_fn [CommRing R] (x y : Fp2Value R) :
    Fp2Value.mul x y = Fp2Value.mul y x := by
  ext <;> simp [Fp2Value.mul] <;> ring

def scalar [Zero R] (x : R) : Fp2Value R := ⟨x, 0⟩

private theorem natCast_succ_fn [CommRing R] (n : Nat) :
    scalar ((n + 1 : Nat) : R) = add (scalar (n : R)) one := by
  ext <;> simp [scalar, add, one, Fp2Value.add]

instance fp2AddCommGroup [CommRing R] : AddCommGroup (Fp2Value R) where
  add_assoc := add_assoc_fn
  zero_add := zero_add_fn
  add_zero := add_zero_fn
  add_comm := add_comm_fn
  neg_add_cancel := neg_add_cancel_fn
  nsmul := nsmulRec
  zsmul := zsmulRec

instance fp2CommMonoid [CommRing R] : CommMonoid (Fp2Value R) where
  mul_assoc := mul_assoc_fn
  one_mul := one_mul_fn
  mul_one := mul_one_fn
  mul_comm := mul_comm_fn
  npow := npowRec

instance fp2CommRing [CommRing R] : CommRing (Fp2Value R) where
  __ := fp2AddCommGroup
  __ := fp2CommMonoid
  natCast := fun n => scalar (n : R)
  natCast_zero := by
    rw [Nat.cast_zero]
    rfl
  natCast_succ := natCast_succ_fn
  zero_mul := zero_mul_fn
  mul_zero := mul_zero_fn
  left_distrib := left_distrib_fn
  right_distrib := right_distrib_fn

theorem scalar_eq_zero_iff [CommRing R] (x : R) :
    scalar x = (0 : Fp2Value R) ↔ x = 0 := by
  constructor
  · intro h
    exact congrArg Fp2Value.x h
  · intro h
    subst x
    rfl

theorem natCast_eq_scalar [CommRing R] (n : Nat) :
    (n : Fp2Value R) = scalar (n : R) := rfl

end Fp2Concrete

end CatBLST.AVX.MatrixNoK
