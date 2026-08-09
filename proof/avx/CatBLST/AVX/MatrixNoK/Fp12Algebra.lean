/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Fp12Raw
import Mathlib.Algebra.CharP.Frobenius
import Mathlib.Tactic

set_option maxHeartbeats 800000

/-!
# Concrete selected MatrixNoK Fp12 algebra

This file gives the tower from `Fp12.lean` its concrete commutative-ring
structure. Keeping the algebra laws separate prevents later correspondence
packets from repeatedly elaborating the large associativity proof.
-/

namespace CatBLST.AVX.MatrixNoK

namespace Fp12Concrete

variable {R : Type*}

def zero [Zero R] : Fp12Value R :=
  ⟨⟨⟨0, 0⟩, ⟨0, 0⟩, ⟨0, 0⟩⟩, ⟨⟨0, 0⟩, ⟨0, 0⟩, ⟨0, 0⟩⟩⟩

def one [Zero R] [One R] : Fp12Value R :=
  ⟨⟨⟨1, 0⟩, ⟨0, 0⟩, ⟨0, 0⟩⟩, ⟨⟨0, 0⟩, ⟨0, 0⟩, ⟨0, 0⟩⟩⟩

def add [Add R] (x y : Fp12Value R) : Fp12Value R :=
  Fp6Value.add x.a y.a |> fun a =>
    ⟨a, Fp6Value.add x.b y.b⟩

def neg [Neg R] (x : Fp12Value R) : Fp12Value R :=
  ⟨⟨⟨-x.a.c0.x, -x.a.c0.y⟩, ⟨-x.a.c1.x, -x.a.c1.y⟩,
      ⟨-x.a.c2.x, -x.a.c2.y⟩⟩,
   ⟨⟨-x.b.c0.x, -x.b.c0.y⟩, ⟨-x.b.c1.x, -x.b.c1.y⟩,
      ⟨-x.b.c2.x, -x.b.c2.y⟩⟩⟩

instance fp12Zero [Zero R] : Zero (Fp12Value R) := ⟨zero⟩
instance fp12One [Zero R] [One R] : One (Fp12Value R) := ⟨one⟩
instance fp12Add [Add R] : Add (Fp12Value R) := ⟨add⟩
instance fp12Neg [Neg R] : Neg (Fp12Value R) := ⟨neg⟩
instance fp12Mul [CommRing R] : Mul (Fp12Value R) := ⟨Fp12Value.mul⟩

private theorem add_assoc_fn [CommRing R] (x y z : Fp12Value R) :
    add (add x y) z = add x (add y z) := by
  ext <;> simp [add, Fp6Value.add, Fp2Value.add] <;> ring

private theorem zero_add_fn [CommRing R] (x : Fp12Value R) :
    add zero x = x := by
  ext <;> simp [zero, add, Fp6Value.add, Fp2Value.add]

private theorem add_zero_fn [CommRing R] (x : Fp12Value R) :
    add x zero = x := by
  ext <;> simp [zero, add, Fp6Value.add, Fp2Value.add]

private theorem add_comm_fn [CommRing R] (x y : Fp12Value R) :
    add x y = add y x := by
  ext <;> simp [add, Fp6Value.add, Fp2Value.add] <;> ring

private theorem neg_add_cancel_fn [CommRing R] (x : Fp12Value R) :
    add (neg x) x = zero := by
  ext <;> simp [neg, add, zero, Fp6Value.add, Fp2Value.add]

private theorem mul_assoc_fn [CommRing R] (x y z : Fp12Value R) :
    Fp12Value.mul (Fp12Value.mul x y) z =
      Fp12Value.mul x (Fp12Value.mul y z) := by
  ext <;>
    simp [Fp12Value.mul, Fp6Value.add, Fp6Value.mul,
      Fp6Value.mulByV, Fp2Value.add, Fp2Value.mul,
      Fp2Value.mulByNonresidue] <;>
    ring

private theorem one_mul_fn [CommRing R] (x : Fp12Value R) :
    Fp12Value.mul one x = x := by
  ext <;>
    simp [one, Fp12Value.mul, Fp6Value.add, Fp6Value.mul,
      Fp6Value.mulByV, Fp2Value.add, Fp2Value.mul,
      Fp2Value.mulByNonresidue]

private theorem mul_one_fn [CommRing R] (x : Fp12Value R) :
    Fp12Value.mul x one = x := by
  ext <;>
    simp [one, Fp12Value.mul, Fp6Value.add, Fp6Value.mul,
      Fp6Value.mulByV, Fp2Value.add, Fp2Value.mul,
      Fp2Value.mulByNonresidue]

private theorem zero_mul_fn [CommRing R] (x : Fp12Value R) :
    Fp12Value.mul zero x = zero := by
  ext <;>
    simp [zero, Fp12Value.mul, Fp6Value.add, Fp6Value.mul,
      Fp6Value.mulByV, Fp2Value.add, Fp2Value.mul,
      Fp2Value.mulByNonresidue]

private theorem mul_zero_fn [CommRing R] (x : Fp12Value R) :
    Fp12Value.mul x zero = zero := by
  ext <;>
    simp [zero, Fp12Value.mul, Fp6Value.add, Fp6Value.mul,
      Fp6Value.mulByV, Fp2Value.add, Fp2Value.mul,
      Fp2Value.mulByNonresidue]

private theorem left_distrib_fn [CommRing R] (x y z : Fp12Value R) :
    Fp12Value.mul x (add y z) =
      add (Fp12Value.mul x y) (Fp12Value.mul x z) := by
  ext <;>
    simp [add, Fp12Value.mul, Fp6Value.add, Fp6Value.mul,
      Fp6Value.mulByV, Fp2Value.add, Fp2Value.mul,
      Fp2Value.mulByNonresidue] <;>
    ring

private theorem right_distrib_fn [CommRing R] (x y z : Fp12Value R) :
    Fp12Value.mul (add x y) z =
      add (Fp12Value.mul x z) (Fp12Value.mul y z) := by
  ext <;>
    simp [add, Fp12Value.mul, Fp6Value.add, Fp6Value.mul,
      Fp6Value.mulByV, Fp2Value.add, Fp2Value.mul,
      Fp2Value.mulByNonresidue] <;>
    ring

private theorem mul_comm_fn [CommRing R] (x y : Fp12Value R) :
    Fp12Value.mul x y = Fp12Value.mul y x := by
  ext <;>
    simp [Fp12Value.mul, Fp6Value.add, Fp6Value.mul,
      Fp6Value.mulByV, Fp2Value.add, Fp2Value.mul,
      Fp2Value.mulByNonresidue] <;>
    ring

def scalar [Zero R] (x : R) : Fp12Value R :=
  ⟨⟨⟨x, 0⟩, ⟨0, 0⟩, ⟨0, 0⟩⟩, ⟨⟨0, 0⟩, ⟨0, 0⟩, ⟨0, 0⟩⟩⟩

private theorem natCast_succ_fn [CommRing R] (n : Nat) :
    scalar ((n + 1 : Nat) : R) =
      add (scalar (n : R)) one := by
  ext <;> simp [scalar, add, one, Fp6Value.add, Fp2Value.add]

instance fp12AddCommGroup [CommRing R] : AddCommGroup (Fp12Value R) where
  add_assoc := add_assoc_fn
  zero_add := zero_add_fn
  add_zero := add_zero_fn
  add_comm := add_comm_fn
  neg_add_cancel := neg_add_cancel_fn
  nsmul := nsmulRec
  zsmul := zsmulRec

instance fp12CommMonoid [CommRing R] : CommMonoid (Fp12Value R) where
  mul_assoc := mul_assoc_fn
  one_mul := one_mul_fn
  mul_one := mul_one_fn
  mul_comm := mul_comm_fn
  npow := npowRec

instance [CommRing R] : CommRing (Fp12Value R) where
  __ := fp12AddCommGroup
  __ := fp12CommMonoid
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
    scalar x = (0 : Fp12Value R) ↔ x = 0 := by
  constructor
  · intro h
    exact congrArg (fun value => value.a.c0.x) h
  · intro h
    subst x
    rfl

theorem natCast_eq_scalar [CommRing R] (n : Nat) :
    (n : Fp12Value R) = scalar (n : R) := rfl

instance [CommRing R] {p : Nat} [CharP R p] :
    CharP (Fp12Value R) p where
  cast_eq_zero_iff n := by
    rw [natCast_eq_scalar, scalar_eq_zero_iff,
      CharP.cast_eq_zero_iff R p]

end Fp12Concrete

end CatBLST.AVX.MatrixNoK
