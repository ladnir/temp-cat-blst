/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Fp12CyclotomicSchedule

set_option maxHeartbeats 800000
set_option maxRecDepth 100000

/-!
Exact source/direct correspondence for the cyclotomic predicate.

Composition is proved through the literal flat twelve-coordinate source
schedule.  This avoids asking the kernel to reduce a doubly nested packed Fp12
term merely to expose one coordinate.
-/

namespace CatBLST.AVX.MatrixNoK

theorem cpp_frobenius2_flat_twice_eq_frobenius4_flat
    (x : FlatFp12 SelectedBaseField) :
    cppFp12Frobenius2Flat (cppFp12Frobenius2Flat x) =
      cppFp12Frobenius4Flat x := by
  funext index
  fin_cases index <;>
    simp [cppFp12Frobenius2Flat, cppFp12Frobenius4Flat,
      mul_frob1_4_twice, mul_frob1_2_twice, mul_frob2_3_twice,
      mul_frob1_3_twice]

theorem cpp_frobenius2_twice_eq_frobenius4 (x : SelectedFp12) :
    cppFp12Frobenius2 (cppFp12Frobenius2 x) =
      cppFp12Frobenius4Direct x := by
  apply unpack_flat_fp12_injective
  funext index
  rw [unpack_cpp_fp12_frobenius2]
  have inner : unpackFlatFp12 (cppFp12Frobenius2 x) =
      cppFp12Frobenius2Flat (unpackFlatFp12 x) := by
    funext coordinate
    exact unpack_cpp_fp12_frobenius2 x coordinate
  rw [inner, congrFun (cpp_frobenius2_flat_twice_eq_frobenius4_flat
    (unpackFlatFp12 x)) index, unpack_frobenius4_direct]

def cppFp12IsCyclotomicSource (x : SelectedFp12) : Prop :=
  cppFp12Frobenius2 x =
    cppFp12Frobenius2 (cppFp12Frobenius2 x) * x

def cppFp12IsCyclotomic (x : SelectedFp12) : Prop :=
  cppFp12Frobenius2Direct x = cppFp12Frobenius4Direct x * x

theorem cpp_fp12_is_cyclotomic_source_iff_direct (x : SelectedFp12) :
    cppFp12IsCyclotomicSource x ↔ cppFp12IsCyclotomic x := by
  unfold cppFp12IsCyclotomicSource cppFp12IsCyclotomic
  rw [cpp_frobenius2_twice_eq_frobenius4, cpp_frobenius2_eq_direct]

def cppFp12IsCyclotomicFlat
    (x : FlatFp12 SelectedBaseField) : Prop :=
  cppFp12Frobenius2Flat x =
    cppFp12LiteralMul 0 (cppFp12Frobenius4Flat x) x

theorem pack_flat_fp12_injective :
    Function.Injective
      (packFlatFp12 : FlatFp12 SelectedBaseField → SelectedFp12) := by
  intro x y equal
  have unpacked := congrArg unpackFlatFp12 equal
  simpa only [unpack_pack_fp12] using unpacked

theorem pack_frobenius2_flat (x : SelectedFp12) :
    packFlatFp12 (cppFp12Frobenius2Flat (unpackFlatFp12 x)) =
      cppFp12Frobenius2Direct x := by
  exact (cpp_frobenius2_eq_direct x).symm.trans rfl

theorem pack_frobenius4_flat (x : SelectedFp12) :
    packFlatFp12 (cppFp12Frobenius4Flat (unpackFlatFp12 x)) =
      cppFp12Frobenius4Direct x := by
  ext <;> rfl

theorem pack_cyclotomic_relation_product (x : SelectedFp12) :
    packFlatFp12
        (cppFp12LiteralMul 0
          (cppFp12Frobenius4Flat (unpackFlatFp12 x))
          (unpackFlatFp12 x)) =
      cppFp12Frobenius4Direct x * x := by
  calc
    packFlatFp12
        (cppFp12LiteralMul 0
          (cppFp12Frobenius4Flat (unpackFlatFp12 x))
          (unpackFlatFp12 x)) =
        cppFp12FlatMul
          (packFlatFp12 (cppFp12Frobenius4Flat (unpackFlatFp12 x)))
          (packFlatFp12 (unpackFlatFp12 x)) :=
      cpp_fp12_literal_mul_packs _ _
    _ = Fp12Value.mul
          (packFlatFp12 (cppFp12Frobenius4Flat (unpackFlatFp12 x)))
          (packFlatFp12 (unpackFlatFp12 x)) :=
      cpp_fp12_flat_mul_correct _ _
    _ = cppFp12Frobenius4Direct x * x := by
      rw [pack_frobenius4_flat, pack_unpack_fp12, fp12_mul_def]

theorem cpp_fp12_is_cyclotomic_iff_flat (x : SelectedFp12) :
    cppFp12IsCyclotomic x ↔
      cppFp12IsCyclotomicFlat (unpackFlatFp12 x) := by
  unfold cppFp12IsCyclotomic cppFp12IsCyclotomicFlat
  constructor
  · intro relation
    apply pack_flat_fp12_injective
    rw [pack_frobenius2_flat, pack_cyclotomic_relation_product]
    exact relation
  · intro relation
    have packed := congrArg packFlatFp12 relation
    rw [pack_frobenius2_flat, pack_cyclotomic_relation_product] at packed
    exact packed

end CatBLST.AVX.MatrixNoK
