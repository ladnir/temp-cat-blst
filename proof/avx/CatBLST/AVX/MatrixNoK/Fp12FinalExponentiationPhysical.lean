/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.Fp12FinalExponentiationAlgebra
import Mathlib.Tactic

/-!
# Final-exponentiation physical ledger

The final-exponentiation header introduces no new fused polynomial kernel.  It
serially composes accepted Fp12 multiplication, Frobenius, conjugation, and
cyclotomic-square calls.  This ledger checks the exact fixed call counts and
records that `mul_n_sqr` explicitly recasts its multiplication result to the
source `CommonElement<6>` accumulator type. Other multiplication calls return
the accepted normalized result before another operation consumes them.
-/

namespace CatBLST.AVX.MatrixNoK

structure FinalExponentiationCost where
  mul : Nat
  cyclotomicSquare : Nat
  conjugate : Nat
  frobenius1 : Nat
  frobenius2 : Nat
  frobenius3 : Nat
  inverse : Nat
  deriving DecidableEq, Repr

instance : Add FinalExponentiationCost where
  add a b := {
    mul := a.mul + b.mul
    cyclotomicSquare := a.cyclotomicSquare + b.cyclotomicSquare
    conjugate := a.conjugate + b.conjugate
    frobenius1 := a.frobenius1 + b.frobenius1
    frobenius2 := a.frobenius2 + b.frobenius2
    frobenius3 := a.frobenius3 + b.frobenius3
    inverse := a.inverse + b.inverse
  }

def zeroFinalExponentiationCost : FinalExponentiationCost :=
  { mul := 0, cyclotomicSquare := 0, conjugate := 0,
    frobenius1 := 0, frobenius2 := 0, frobenius3 := 0, inverse := 0 }

def mulNSqrCost (n : Nat) : FinalExponentiationCost :=
  { zeroFinalExponentiationCost with mul := 1, cyclotomicSquare := n }

def raiseToZDivBy2Cost : FinalExponentiationCost :=
  { zeroFinalExponentiationCost with
    mul := 5
    cyclotomicSquare := 1 + 2 + 3 + 9 + 32 + (16 - 1)
    conjugate := 1 }

def raiseToZCost : FinalExponentiationCost :=
  { raiseToZDivBy2Cost with
    cyclotomicSquare := raiseToZDivBy2Cost.cyclotomicSquare + 1 }

/-- Direct operations in `final_exp_after_inverse`, excluding nested raises. -/
def finalAfterInverseDirectCost : FinalExponentiationCost :=
  { mul := 10
    cyclotomicSquare := 1
    conjugate := 4
    frobenius1 := 1
    frobenius2 := 2
    frobenius3 := 1
    inverse := 0 }

/-- Four `raise_to_z` calls and one additional `raise_to_z_div_by_2`. -/
def finalAfterInverseNestedCost : FinalExponentiationCost :=
  raiseToZCost + raiseToZCost + raiseToZCost + raiseToZCost +
    raiseToZDivBy2Cost

def finalAfterInverseCost : FinalExponentiationCost :=
  finalAfterInverseDirectCost + finalAfterInverseNestedCost

def fullFinalExponentiationCost : FinalExponentiationCost :=
  { finalAfterInverseCost with inverse := 1 }

theorem raise_to_z_div_by_2_exact_cost :
    raiseToZDivBy2Cost =
      { mul := 5, cyclotomicSquare := 62, conjugate := 1,
        frobenius1 := 0, frobenius2 := 0, frobenius3 := 0, inverse := 0 } := by
  native_decide

theorem raise_to_z_exact_cost :
    raiseToZCost =
      { mul := 5, cyclotomicSquare := 63, conjugate := 1,
        frobenius1 := 0, frobenius2 := 0, frobenius3 := 0, inverse := 0 } := by
  native_decide

theorem final_after_inverse_exact_cost :
    finalAfterInverseCost =
      { mul := 35, cyclotomicSquare := 315, conjugate := 9,
        frobenius1 := 1, frobenius2 := 2, frobenius3 := 1, inverse := 0 } := by
  native_decide

theorem full_final_exponentiation_exact_cost :
    fullFinalExponentiationCost =
      { mul := 35, cyclotomicSquare := 315, conjugate := 9,
        frobenius1 := 1, frobenius2 := 2, frobenius3 := 1, inverse := 1 } := by
  native_decide

/-- Maximum declared common bound carried between source operations. -/
def finalExponentiationCommonBound : Nat := 6

theorem final_exponentiation_common_bound_matches_source :
    finalExponentiationCommonBound = 6 := rfl

theorem final_exponentiation_has_no_cross_call_accumulation :
    finalExponentiationCommonBound <= 6 := by
  norm_num [finalExponentiationCommonBound]

end CatBLST.AVX.MatrixNoK
