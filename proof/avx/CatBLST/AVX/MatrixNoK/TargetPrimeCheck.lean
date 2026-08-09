/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Generated
import Mathlib.Tactic

namespace CatBLST.AVX.MatrixNoK

structure LucasCheck where
  modulus : Nat
  witness : Nat
  factors : List Nat
deriving DecidableEq

def lucasCheckValid (check : LucasCheck) : Bool :=
  decide ((check.witness : ZMod check.modulus) ^
      (check.modulus - 1) = 1) &&
    check.factors.all (fun factor =>
      decide ((check.witness : ZMod check.modulus) ^
        ((check.modulus - 1) / factor) ≠ 1))

theorem lucas_check_full_order {check : LucasCheck}
    (valid : lucasCheckValid check = true) :
    (check.witness : ZMod check.modulus) ^
      (check.modulus - 1) = 1 := by
  exact of_decide_eq_true (Bool.and_eq_true_iff.mp valid).1

theorem lucas_check_factor_order {check : LucasCheck}
    (valid : lucasCheckValid check = true) {factor : Nat}
    (member : factor ∈ check.factors) :
    (check.witness : ZMod check.modulus) ^
      ((check.modulus - 1) / factor) ≠ 1 := by
  have checked :=
    (List.all_eq_true.mp (Bool.and_eq_true_iff.mp valid).2) factor member
  exact of_decide_eq_true checked

end CatBLST.AVX.MatrixNoK
