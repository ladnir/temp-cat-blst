/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.CurveBoothExecutor
import CatBLST.AVX.MatrixNoK.CurveWnafExecutor
import CatBLST.AVX.MatrixNoK.CurveGenerators
import Mathlib.Algebra.Module.Torsion.Basic
import Mathlib.Algebra.Group.Subgroup.ZPowers.Basic
import Mathlib.Tactic

/-!
# Selected G1/G2 scalar-multiplication dispatch

This file composes the source-shaped GLV/GLS scalar splits and executor proofs
with an actual additive group.  The group layer is intentionally a
`ZMod blsScalarOrder` module: that is precisely the algebra available for a
point in the prime-order subgroup, and it makes reduction modulo `r` explicit
rather than hiding it in integer arithmetic.

The four point executors below are the mathematical meanings of the selected
production owners:

* G1 constant-time: two radix-`z^2` components and the fixed Booth schedule;
* G1 variable-time: the same components and the selected width-5 wNAF;
* G2 constant-time: four radix-`z` components and the fixed Booth schedule;
* G2 variable-time: the same components and the selected width-4 wNAF.

Each component acts on the input point by the exponent of the corresponding
source table.  The preceding files prove the literal recoders and table
schedules produce those exponents.  Consequently, the theorems here are the
missing composition from those source schedules to scalar multiplication in
the established projective group law.
-/

namespace CatBLST.AVX.MatrixNoK

/-- The public boundary first ignores bits above `nbits`. -/
def curveMaskedScalar (scalar nbits : Nat) : Nat :=
  scalar % (2 ^ nbits)

/-- The selected GLV/GLS owners consume the canonical representative modulo `r`. -/
def curvePreparedScalar (scalar nbits : Nat) : Nat :=
  curveMaskedScalar scalar nbits % blsScalarOrder

theorem curve_prepared_scalar_canonical (scalar nbits : Nat) :
    curvePreparedScalar scalar nbits < blsScalarOrder := by
  exact Nat.mod_lt _ (by native_decide)

theorem curve_prepared_scalar_cast (scalar nbits : Nat) :
    (curvePreparedScalar scalar nbits : CurveExponent) =
      (curveMaskedScalar scalar nbits : CurveExponent) := by
  simp [curvePreparedScalar]

/-- Exact symbolic result of the selected G1 constant-time component schedule. -/
def g1BoothDispatchExponent (scalar : Nat) : CurveExponent :=
  Finset.univ.sum fun index : Fin 2 =>
    boothComponentExponent 26 curveRadixSquared
      (g1GlvBasisExponent index) (g1GlvSplit scalar index)

/-- Exact symbolic result of the selected G2 constant-time component schedule. -/
def g2BoothDispatchExponent (scalar : Nat) : CurveExponent :=
  Finset.univ.sum fun index : Fin 4 =>
    boothComponentExponent 13 curveRadix
      (g2GlsBasisExponent index) (g2GlsSplit scalar index)

/-- Exact symbolic result of the selected G1 variable-time component schedule. -/
def g1WnafDispatchExponent (scalar : Nat) : CurveExponent :=
  Finset.univ.sum fun index : Fin 2 =>
    wnafComponentExponent 5 129
      (g1GlvBasisExponent index) (g1GlvSplit scalar index)

/-- Exact symbolic result of the selected G2 variable-time component schedule. -/
def g2WnafDispatchExponent (scalar : Nat) : CurveExponent :=
  Finset.univ.sum fun index : Fin 4 =>
    wnafComponentExponent 4 65
      (g2GlsBasisExponent index) (g2GlsSplit scalar index)

theorem g1_booth_dispatch_exponent_correct (scalar : Nat)
    (canonical : scalar < blsScalarOrder) :
    g1BoothDispatchExponent scalar = (scalar : CurveExponent) := by
  exact g1_booth_executor_correct scalar canonical

theorem g2_booth_dispatch_exponent_correct (scalar : Nat)
    (canonical : scalar < blsScalarOrder) :
    g2BoothDispatchExponent scalar = (scalar : CurveExponent) := by
  exact g2_booth_executor_correct scalar canonical

theorem g1_wnaf_dispatch_exponent_correct (scalar : Nat)
    (canonical : scalar < blsScalarOrder) :
    g1WnafDispatchExponent scalar = (scalar : CurveExponent) := by
  exact g1_wnaf_executor_correct scalar canonical

theorem g2_wnaf_dispatch_exponent_correct (scalar : Nat)
    (canonical : scalar < blsScalarOrder) :
    g2WnafDispatchExponent scalar = (scalar : CurveExponent) := by
  exact g2_wnaf_executor_correct scalar canonical

section GroupExecution

variable {G : Type*} [AddCommGroup G] [Module CurveExponent G]

/-- Point-level meaning of the selected G1 constant-time table executor. -/
def g1BoothGroupExecute (point : G) (scalar : Nat) : G :=
  Finset.univ.sum fun index : Fin 2 =>
    boothComponentExponent 26 curveRadixSquared
      (g1GlvBasisExponent index) (g1GlvSplit scalar index) • point

/-- Point-level meaning of the selected G2 constant-time table executor. -/
def g2BoothGroupExecute (point : G) (scalar : Nat) : G :=
  Finset.univ.sum fun index : Fin 4 =>
    boothComponentExponent 13 curveRadix
      (g2GlsBasisExponent index) (g2GlsSplit scalar index) • point

/-- Point-level meaning of the selected G1 variable-time table executor. -/
def g1WnafGroupExecute (point : G) (scalar : Nat) : G :=
  Finset.univ.sum fun index : Fin 2 =>
    wnafComponentExponent 5 129
      (g1GlvBasisExponent index) (g1GlvSplit scalar index) • point

/-- Point-level meaning of the selected G2 variable-time table executor. -/
def g2WnafGroupExecute (point : G) (scalar : Nat) : G :=
  Finset.univ.sum fun index : Fin 4 =>
    wnafComponentExponent 4 65
      (g2GlsBasisExponent index) (g2GlsSplit scalar index) • point

theorem g1_booth_group_execute_correct (point : G) (scalar : Nat)
    (canonical : scalar < blsScalarOrder) :
    g1BoothGroupExecute point scalar = (scalar : CurveExponent) • point := by
  rw [g1BoothGroupExecute, ← Finset.sum_smul,
    g1_booth_executor_correct scalar canonical]

theorem g2_booth_group_execute_correct (point : G) (scalar : Nat)
    (canonical : scalar < blsScalarOrder) :
    g2BoothGroupExecute point scalar = (scalar : CurveExponent) • point := by
  rw [g2BoothGroupExecute, ← Finset.sum_smul,
    g2_booth_executor_correct scalar canonical]

theorem g1_wnaf_group_execute_correct (point : G) (scalar : Nat)
    (canonical : scalar < blsScalarOrder) :
    g1WnafGroupExecute point scalar = (scalar : CurveExponent) • point := by
  rw [g1WnafGroupExecute, ← Finset.sum_smul,
    g1_wnaf_executor_correct scalar canonical]

theorem g2_wnaf_group_execute_correct (point : G) (scalar : Nat)
    (canonical : scalar < blsScalarOrder) :
    g2WnafGroupExecute point scalar = (scalar : CurveExponent) • point := by
  rw [g2WnafGroupExecute, ← Finset.sum_smul,
    g2_wnaf_executor_correct scalar canonical]

/-- Selected public G1 `mult` after masking and canonical preparation. -/
def selectedG1Mult (point : G) (scalar nbits : Nat) : G :=
  g1BoothGroupExecute point (curvePreparedScalar scalar nbits)

/-- Selected public G1 `mult_vartime` after masking and canonical preparation. -/
def selectedG1MultVartime (point : G) (scalar nbits : Nat) : G :=
  g1WnafGroupExecute point (curvePreparedScalar scalar nbits)

/-- Selected public G2 `mult` after masking and canonical preparation. -/
def selectedG2Mult (point : G) (scalar nbits : Nat) : G :=
  g2BoothGroupExecute point (curvePreparedScalar scalar nbits)

/-- Selected public G2 `mult_vartime` after masking and canonical preparation. -/
def selectedG2MultVartime (point : G) (scalar nbits : Nat) : G :=
  g2WnafGroupExecute point (curvePreparedScalar scalar nbits)

theorem selected_g1_mult_correct (point : G) (scalar nbits : Nat) :
    selectedG1Mult point scalar nbits =
      (curveMaskedScalar scalar nbits : CurveExponent) • point := by
  rw [selectedG1Mult,
    g1_booth_group_execute_correct _ _
      (curve_prepared_scalar_canonical scalar nbits),
    curve_prepared_scalar_cast]

theorem selected_g1_mult_vartime_correct (point : G) (scalar nbits : Nat) :
    selectedG1MultVartime point scalar nbits =
      (curveMaskedScalar scalar nbits : CurveExponent) • point := by
  rw [selectedG1MultVartime,
    g1_wnaf_group_execute_correct _ _
      (curve_prepared_scalar_canonical scalar nbits),
    curve_prepared_scalar_cast]

theorem selected_g2_mult_correct (point : G) (scalar nbits : Nat) :
    selectedG2Mult point scalar nbits =
      (curveMaskedScalar scalar nbits : CurveExponent) • point := by
  rw [selectedG2Mult,
    g2_booth_group_execute_correct _ _
      (curve_prepared_scalar_canonical scalar nbits),
    curve_prepared_scalar_cast]

theorem selected_g2_mult_vartime_correct (point : G) (scalar nbits : Nat) :
    selectedG2MultVartime point scalar nbits =
      (curveMaskedScalar scalar nbits : CurveExponent) • point := by
  rw [selectedG2MultVartime,
    g2_wnaf_group_execute_correct _ _
      (curve_prepared_scalar_canonical scalar nbits),
    curve_prepared_scalar_cast]

end GroupExecution

/-!
The concrete carriers below are the `r`-torsion subgroups of the selected
Mathlib projective curves.  Their addition is therefore the projective group
law proved equivalent to `cppCompleteAdd` in `CurveGroupAdd.lean`, while their
scalar action is canonically over `ZMod r`.
-/

noncomputable abbrev SelectedG1PrimeSubgroup :=
  Submodule.torsionBy ℤ selectedG1Weierstrass.toProjective.Point
    (blsScalarOrder : ℤ)

noncomputable abbrev SelectedG2PrimeSubgroup :=
  Submodule.torsionBy ℤ selectedG2Weierstrass.toProjective.Point
    (blsScalarOrder : ℤ)

noncomputable instance selectedG1PrimeSubgroupModule :
    Module CurveExponent SelectedG1PrimeSubgroup :=
  AddSubgroup.torsionBy.zmodModule

noncomputable instance selectedG2PrimeSubgroupModule :
    Module CurveExponent SelectedG2PrimeSubgroup :=
  AddSubgroup.torsionBy.zmodModule

/-!
CatBLST's semantic G1/G2 domains are the cyclic subgroups generated by the
standard generators, not every abstract `r`-torsion point over an algebraic
closure.  The concrete generator proofs above show these subgroups are killed
by `r`, so they inherit the same canonical `ZMod r` action without a
point-count assumption.
-/

noncomputable abbrev SelectedG1CyclicSubgroup :=
  AddSubgroup.zmultiples selectedG1GeneratorPoint

noncomputable abbrev SelectedG2CyclicSubgroup :=
  AddSubgroup.zmultiples selectedG2GeneratorPoint

noncomputable instance selectedG1CyclicSubgroupModule :
    Module CurveExponent SelectedG1CyclicSubgroup :=
  AddCommGroup.zmodModule (n := blsScalarOrder) (by
    intro point
    apply Subtype.ext
    change blsScalarOrder •
      (point : selectedG1Weierstrass.toProjective.Point) = 0
    rcases AddSubgroup.mem_zmultiples_iff.mp point.property with ⟨scalar, span⟩
    rw [← span]
    cases scalar with
    | ofNat scalar =>
        change blsScalarOrder • (scalar • selectedG1GeneratorPoint) = 0
        rw [nsmul_left_comm, selected_g1_generator_point_order_divides_r]
        simp
    | negSucc scalar =>
        change blsScalarOrder • (-((scalar + 1) • selectedG1GeneratorPoint)) = 0
        have positive : blsScalarOrder •
            ((scalar + 1) • selectedG1GeneratorPoint) = 0 := by
          rw [nsmul_left_comm,
            selected_g1_generator_point_order_divides_r]
          simp
        simpa using congrArg Neg.neg positive)

noncomputable instance selectedG2CyclicSubgroupModule :
    Module CurveExponent SelectedG2CyclicSubgroup :=
  AddCommGroup.zmodModule (n := blsScalarOrder) (by
    intro point
    apply Subtype.ext
    change blsScalarOrder •
      (point : selectedG2Weierstrass.toProjective.Point) = 0
    rcases AddSubgroup.mem_zmultiples_iff.mp point.property with ⟨scalar, span⟩
    rw [← span]
    cases scalar with
    | ofNat scalar =>
        change blsScalarOrder • (scalar • selectedG2GeneratorPoint) = 0
        rw [nsmul_left_comm, selected_g2_generator_point_order_divides_r]
        simp
    | negSucc scalar =>
        change blsScalarOrder • (-((scalar + 1) • selectedG2GeneratorPoint)) = 0
        have positive : blsScalarOrder •
            ((scalar + 1) • selectedG2GeneratorPoint) = 0 := by
          rw [nsmul_left_comm,
            selected_g2_generator_point_order_divides_r]
          simp
        simpa using congrArg Neg.neg positive)

theorem selected_g1_projective_mult_correct
    (point : SelectedG1PrimeSubgroup) (scalar nbits : Nat) :
    selectedG1Mult point scalar nbits =
      (curveMaskedScalar scalar nbits : CurveExponent) • point :=
  selected_g1_mult_correct point scalar nbits

theorem selected_g1_projective_mult_vartime_correct
    (point : SelectedG1PrimeSubgroup) (scalar nbits : Nat) :
    selectedG1MultVartime point scalar nbits =
      (curveMaskedScalar scalar nbits : CurveExponent) • point :=
  selected_g1_mult_vartime_correct point scalar nbits

theorem selected_g2_projective_mult_correct
    (point : SelectedG2PrimeSubgroup) (scalar nbits : Nat) :
    selectedG2Mult point scalar nbits =
      (curveMaskedScalar scalar nbits : CurveExponent) • point :=
  selected_g2_mult_correct point scalar nbits

theorem selected_g2_projective_mult_vartime_correct
    (point : SelectedG2PrimeSubgroup) (scalar nbits : Nat) :
    selectedG2MultVartime point scalar nbits =
      (curveMaskedScalar scalar nbits : CurveExponent) • point :=
  selected_g2_mult_vartime_correct point scalar nbits

theorem selected_g1_cyclic_mult_correct
    (point : SelectedG1CyclicSubgroup) (scalar nbits : Nat) :
    selectedG1Mult point scalar nbits =
      (curveMaskedScalar scalar nbits : CurveExponent) • point :=
  selected_g1_mult_correct point scalar nbits

theorem selected_g1_cyclic_mult_vartime_correct
    (point : SelectedG1CyclicSubgroup) (scalar nbits : Nat) :
    selectedG1MultVartime point scalar nbits =
      (curveMaskedScalar scalar nbits : CurveExponent) • point :=
  selected_g1_mult_vartime_correct point scalar nbits

theorem selected_g2_cyclic_mult_correct
    (point : SelectedG2CyclicSubgroup) (scalar nbits : Nat) :
    selectedG2Mult point scalar nbits =
      (curveMaskedScalar scalar nbits : CurveExponent) • point :=
  selected_g2_mult_correct point scalar nbits

theorem selected_g2_cyclic_mult_vartime_correct
    (point : SelectedG2CyclicSubgroup) (scalar nbits : Nat) :
    selectedG2MultVartime point scalar nbits =
      (curveMaskedScalar scalar nbits : CurveExponent) • point :=
  selected_g2_mult_vartime_correct point scalar nbits

end CatBLST.AVX.MatrixNoK
