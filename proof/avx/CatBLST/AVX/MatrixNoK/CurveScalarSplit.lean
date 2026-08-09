/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.GtScalarSplit
import Mathlib.Tactic

set_option maxHeartbeats 800000

/-!
# G1 GLV and G2 GLS scalar decomposition

This file models the exact unsigned decompositions used by the selected AVX
`GLVMult::multiply` and `GLSMult::multiply` owners after public scalar
preparation.  The G1 source divides once by `z^2`; G2 divides by `z^2` and
then divides both halves by `z`.  Pair ordering is remainder followed by
quotient, matching `div_by_zz_into` and `div_by_z_into`.

The semantic endpoint is the prime-order subgroup exponent ring.  The concrete
coordinate endomorphisms are connected to these basis exponents by the
subgroup precondition in the curve-operation packet.
-/

namespace CatBLST.AVX.MatrixNoK

abbrev CurveExponent := ZMod blsScalarOrder

/-- Positive magnitude of the negative BLS12-381 parameter. -/
def curveRadix : Nat := gtRadix

/-- `z^2`, the divisor and G1 endomorphism eigenvalue used by the source. -/
def curveRadixSquared : Nat := curveRadix ^ 2

theorem curve_radix_positive : 0 < curveRadix := by
  exact gt_radix_positive

theorem curve_radix_squared_positive : 0 < curveRadixSquared := by
  exact pow_pos curve_radix_positive 2

theorem curve_scalar_order_lt_radix_fourth :
    blsScalarOrder < curveRadix ^ 4 := by
  exact gt_scalar_order_lt_radix_fourth

/-- Exact two-component output of the G1 `div_by_zz` split. -/
def g1GlvSplit (scalar : Nat) : Fin 2 -> Nat :=
  ![scalar % curveRadixSquared, scalar / curveRadixSquared]

theorem g1_glv_split_recompose (scalar : Nat) :
    scalar =
      g1GlvSplit scalar 0 + curveRadixSquared * g1GlvSplit scalar 1 := by
  simpa [g1GlvSplit, curveRadixSquared] using
    (Nat.mod_add_div scalar curveRadixSquared).symm

theorem g1_canonical_split_bounds (scalar : Nat)
    (canonical : scalar < blsScalarOrder) :
    forall index : Fin 2, g1GlvSplit scalar index < curveRadixSquared := by
  intro index
  fin_cases index
  · exact Nat.mod_lt _ curve_radix_squared_positive
  · simp only [g1GlvSplit]
    apply (Nat.div_lt_iff_lt_mul curve_radix_squared_positive).2
    calc
      scalar < blsScalarOrder := canonical
      _ < curveRadix ^ 4 := curve_scalar_order_lt_radix_fourth
      _ = curveRadixSquared * curveRadixSquared := by
        simp [curveRadixSquared]
        ring

/-- The selected G1 table bases: `P` and the `z^2` eigenpoint. -/
def g1GlvBasisExponent : Fin 2 -> CurveExponent :=
  ![1, curveRadixSquared]

def g1GlvSplitExponent (scalar : Nat) : CurveExponent :=
  Finset.univ.sum fun index : Fin 2 =>
    (g1GlvSplit scalar index : CurveExponent) * g1GlvBasisExponent index

theorem g1_glv_split_exponent_correct (scalar : Nat) :
    g1GlvSplitExponent scalar = (scalar : CurveExponent) := by
  have recomposed := g1_glv_split_recompose scalar
  have castRecomposed := congrArg
    (fun value : Nat => (value : CurveExponent)) recomposed
  push_cast at castRecomposed
  simp [g1GlvSplitExponent, g1GlvBasisExponent, Fin.sum_univ_two]
  simpa [mul_comm] using castRecomposed.symm

/-- Exact four-component output of the selected G2 split. -/
def g2GlsSplit (scalar : Nat) : Fin 4 -> Nat :=
  gtRadixZSplit scalar

theorem g2_gls_split_recompose (scalar : Nat) :
    scalar =
      g2GlsSplit scalar 0 +
      curveRadix * g2GlsSplit scalar 1 +
      curveRadix ^ 2 * g2GlsSplit scalar 2 +
      curveRadix ^ 3 * g2GlsSplit scalar 3 := by
  exact gt_radix_z_split_recompose scalar

theorem g2_canonical_split_bounds (scalar : Nat)
    (canonical : scalar < blsScalarOrder) :
    forall index : Fin 4, g2GlsSplit scalar index < curveRadix := by
  apply gt_radix_z_split_digit_bounds
  exact canonical.trans curve_scalar_order_lt_radix_fourth

/-- The selected G2 bases after each application of `minus_psi`. -/
def g2GlsBasisExponent : Fin 4 -> CurveExponent :=
  ![1, curveRadix, curveRadix ^ 2, curveRadix ^ 3]

def g2GlsSplitExponent (scalar : Nat) : CurveExponent :=
  Finset.univ.sum fun index : Fin 4 =>
    (g2GlsSplit scalar index : CurveExponent) * g2GlsBasisExponent index

theorem g2_gls_split_exponent_correct (scalar : Nat) :
    g2GlsSplitExponent scalar = (scalar : CurveExponent) := by
  have recomposed := g2_gls_split_recompose scalar
  have castRecomposed := congrArg
    (fun value : Nat => (value : CurveExponent)) recomposed
  push_cast at castRecomposed
  simp [g2GlsSplitExponent, g2GlsBasisExponent, Fin.sum_univ_four]
  simpa [mul_comm, mul_left_comm, mul_assoc] using castRecomposed.symm

end CatBLST.AVX.MatrixNoK
