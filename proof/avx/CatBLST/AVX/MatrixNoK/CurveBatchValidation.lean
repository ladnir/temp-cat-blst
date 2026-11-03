/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.CurveGroupAdd
import CatBLST.AVX.MatrixNoK.CanonicalExport
import CatBLST.AVX.MatrixNoK.CurveSubgroup
import CatBLST.AVX.MatrixNoK.CurveGeneratorOrder
import Mathlib.GroupTheory.Coset.Card

/-!
# Curve validity through AVX batch addition

The scalar boundary treats every field-zero `Z` as infinity. MatrixNoK `dadd`
selects an operand only when its *RNS representation* of `Z` is exactly zero.
These predicates differ. This model keeps both observations and proves that
batch addition preserves a valid projective point, or an exact-zero infinity
that `dadd` will select away on the next addition. Exact RNS zero is proved to
decode to field zero using the existing ConvertFrom model. Identifying each
physical add output with the reduced field formula remains a separate proof
boundary.
-/

namespace CatBLST.AVX.MatrixNoK

/-- The literal field-level operation order of `matrixnok_curve.h::add_point`.
The six products are the first `batch_reduce`; the six cross/linear terms
are `batch_expand`; the final coordinates are `batch_reduce_expand`. -/
def cppMatrixNoKAddSchedule [CommRing K] (b : K)
    (left right : HomogeneousPoint K) : HomogeneousPoint K :=
  let leftXY := left.x + left.y
  let leftXZ := left.x + left.z
  let leftYZ := left.y + left.z
  let rightXY := right.x + right.y
  let rightYZ := right.y + right.z
  let rightXZ := right.x + right.z
  let xx := left.x * right.x
  let yy := left.y * right.y
  let zz := left.z * right.z
  let xy := leftXY * rightXY
  let yz := leftYZ * rightYZ
  let xz := leftXZ * rightXZ
  let xxPlusYY := xx + yy
  let yyPlusZZ := yy + zz
  let xxPlusZZ := xx + zz
  let xzCross := xz - xxPlusZZ
  let twoXX := xx + xx
  let b3ZZ := 3 * b * zz
  let yzCross := yz - yyPlusZZ
  let xyCross := xy - xxPlusYY
  let b3XZCross := 3 * b * xzCross
  let yyMinusB3ZZ := yy - b3ZZ
  let threeXX := xx + twoXX
  let yyPlusB3ZZ := yy + b3ZZ
  let xNegative := yzCross * -b3XZCross
  let xPositive := yyMinusB3ZZ * xyCross
  let yLeft := threeXX * b3XZCross
  let yRight := yyMinusB3ZZ * yyPlusB3ZZ
  let zLeft := yzCross * yyPlusB3ZZ
  let zRight := threeXX * xyCross
  ⟨xPositive + xNegative, yRight + yLeft, zLeft + zRight⟩

theorem cpp_matrixnok_add_schedule_eq_complete_add [CommRing K] (b : K)
    (left right : HomogeneousPoint K) :
    cppMatrixNoKAddSchedule b left right = cppCompleteAdd b left right := by
  ext <;> simp only [cppMatrixNoKAddSchedule, cppCompleteAdd] <;> ring

/-- The physical two-basis coordinate inspected by `ExpandedElement::is_zero`.
`convert_from_raw` exports the second basis, so its field observation is
derived from `m2` by the existing fixed ConvertFrom model. -/
structure BatchRnsFp where
  m1 : Fin 8 → Nat
  m2 : Fin 8 → Nat

def BatchRnsFp.exactZero (value : BatchRnsFp) : Prop :=
  (∀ lane, value.m1 lane = 0) ∧ (∀ lane, value.m2 lane = 0)

def BatchRnsFp.fieldValue (value : BatchRnsFp) : SelectedCurveFp :=
  generatedConvertFromValue value.m2

theorem batch_rns_fp_exact_zero_implies_field_zero
    (value : BatchRnsFp) (exact : value.exactZero) :
    value.fieldValue = 0 := by
  have m2Zero : value.m2 = fun _ => 0 := funext exact.2
  simp [BatchRnsFp.fieldValue, generatedConvertFromValue,
    generatedConvertFromEstimatedQuotient, m2Zero]

/-- The G2 `FP2::is_zero` predicate conjoins the exact-zero predicates of
its two base-field components. -/
structure BatchRnsFp2 where
  x : BatchRnsFp
  y : BatchRnsFp

def BatchRnsFp2.exactZero (value : BatchRnsFp2) : Prop :=
  value.x.exactZero ∧ value.y.exactZero

def BatchRnsFp2.fieldValue (value : BatchRnsFp2) : SelectedMillerFp2 :=
  ⟨value.x.fieldValue, value.y.fieldValue⟩

theorem batch_rns_fp2_exact_zero_implies_field_zero
    (value : BatchRnsFp2) (exact : value.exactZero) :
    value.fieldValue = 0 := by
  have xZero := batch_rns_fp_exact_zero_implies_field_zero value.x exact.1
  have yZero := batch_rns_fp_exact_zero_implies_field_zero value.y exact.2
  change (⟨value.x.fieldValue, value.y.fieldValue⟩ :
    SelectedMillerFp2) = ⟨0, 0⟩
  simp [xZero, yZero]

structure BatchRnsPoint (K : Type*) where
  fieldPoint : HomogeneousPoint K
  exactZero : Bool

def BatchRnsPoint.exactZeroSound [Zero K] (point : BatchRnsPoint K) : Prop :=
  point.exactZero = true → point.fieldPoint.z = 0

/-- A G1 point whose semantic coordinates and exact-zero test are both
computed from the same physical RNS coordinates. -/
noncomputable def BatchRnsPoint.fromFpCoordinates (x y z : BatchRnsFp) :
    BatchRnsPoint SelectedCurveFp := by
  classical
  exact ⟨⟨x.fieldValue, y.fieldValue, z.fieldValue⟩, decide z.exactZero⟩

theorem batch_rns_g1_exact_zero_sound (x y z : BatchRnsFp) :
    (BatchRnsPoint.fromFpCoordinates x y z).exactZeroSound := by
  classical
  intro exact
  have exactPhysical : z.exactZero := by
    apply of_decide_eq_true
    simpa [BatchRnsPoint.fromFpCoordinates] using exact
  exact batch_rns_fp_exact_zero_implies_field_zero z
    exactPhysical

/-- The same concrete construction for the two-component G2 coordinate. -/
noncomputable def BatchRnsPoint.fromFp2Coordinates (x y z : BatchRnsFp2) :
    BatchRnsPoint SelectedMillerFp2 := by
  classical
  exact ⟨⟨x.fieldValue, y.fieldValue, z.fieldValue⟩, decide z.exactZero⟩

theorem batch_rns_g2_exact_zero_sound (x y z : BatchRnsFp2) :
    (BatchRnsPoint.fromFp2Coordinates x y z).exactZeroSound := by
  classical
  intro exact
  have exactPhysical : z.exactZero := by
    apply of_decide_eq_true
    simpa [BatchRnsPoint.fromFp2Coordinates] using exact
  exact batch_rns_fp2_exact_zero_implies_field_zero z
    exactPhysical

/-- The generic `rawSound` premise follows for a physical G1 result once
the raw addition's exported Z is identified with the reduced-field formula. -/
theorem batch_g1_raw_sound_of_physical (b : SelectedCurveFp)
    (left right : BatchRnsPoint SelectedCurveFp)
    (rawX rawY rawZ : BatchRnsFp)
    (resultZ : (cppCompleteAdd b left.fieldPoint right.fieldPoint).z =
      rawZ.fieldValue) :
    (BatchRnsPoint.fromFpCoordinates rawX rawY rawZ).exactZero = true →
      (cppCompleteAdd b left.fieldPoint right.fieldPoint).z = 0 := by
  intro exact
  rw [resultZ]
  exact batch_rns_g1_exact_zero_sound rawX rawY rawZ exact

/-- The corresponding physical-result bridge for G2. -/
theorem batch_g2_raw_sound_of_physical (b : SelectedMillerFp2)
    (left right : BatchRnsPoint SelectedMillerFp2)
    (rawX rawY rawZ : BatchRnsFp2)
    (resultZ : (cppCompleteAdd b left.fieldPoint right.fieldPoint).z =
      rawZ.fieldValue) :
    (BatchRnsPoint.fromFp2Coordinates rawX rawY rawZ).exactZero = true →
      (cppCompleteAdd b left.fieldPoint right.fieldPoint).z = 0 := by
  intro exact
  rw [resultZ]
  exact batch_rns_g2_exact_zero_sound rawX rawY rawZ exact

/-- The scalar validator's curve interpretation of a projective input. -/
def BatchRnsPoint.semanticallyOnCurve [CommRing K] (b : K)
    (point : BatchRnsPoint K) : Prop :=
  point.fieldPoint.z = 0 ∨ point.fieldPoint.onCurve b

/-- Scalar `on_curve` after the homogeneous-to-Jacobian conversion
`(X,Y,Z) ↦ (XZ,YZ²,Z)`, including its unconditional `Z = 0` acceptance. -/
def BatchRnsPoint.scalarConvertedOnCurve [CommRing K] (b : K)
    (point : BatchRnsPoint K) : Prop :=
  point.fieldPoint.z = 0 ∨
    (point.fieldPoint.y * point.fieldPoint.z ^ 2) ^ 2 =
      (point.fieldPoint.x * point.fieldPoint.z) ^ 3 +
        b * point.fieldPoint.z ^ 6

theorem batch_scalar_curve_iff_semantic_curve [Field K] (b : K)
    (point : BatchRnsPoint K) :
    point.scalarConvertedOnCurve b ↔ point.semanticallyOnCurve b := by
  by_cases zZero : point.fieldPoint.z = 0
  · simp [BatchRnsPoint.scalarConvertedOnCurve,
      BatchRnsPoint.semanticallyOnCurve, zZero]
  constructor
  · intro scalarValid
    right
    have jacobian := (Or.resolve_left scalarValid zZero)
    have scaledResidual : point.fieldPoint.z ^ 3 *
        (point.fieldPoint.y ^ 2 * point.fieldPoint.z -
          point.fieldPoint.x ^ 3 - b * point.fieldPoint.z ^ 3) = 0 := by
      linear_combination jacobian
    have residual := (mul_eq_zero.mp scaledResidual).resolve_left
      (pow_ne_zero 3 zZero)
    simp only [HomogeneousPoint.onCurve]
    linear_combination residual
  · intro semanticValid
    right
    have curveEquation := Or.resolve_left semanticValid zZero
    simp only [HomogeneousPoint.onCurve] at curveEquation
    linear_combination point.fieldPoint.z ^ 3 * curveEquation

/-- The invariant needed before a representation-sensitive `dadd`. -/
def BatchRnsPoint.additionSafe [CommRing K] (b : K)
    (point : BatchRnsPoint K) : Prop :=
  (point.exactZero = true ∧ point.fieldPoint.z = 0) ∨
    point.fieldPoint.projectivelyValid b

theorem batch_addition_safe_semantically_on_curve [Field K] (b : K)
    (point : BatchRnsPoint K) (safe : point.additionSafe b) :
    point.semanticallyOnCurve b := by
  rcases safe with exact | valid
  · exact Or.inl exact.2
  · exact Or.inr valid.1

/-- What prevalidation establishes when every semantic infinity has an exact
zero RNS `Z`. This is the normal randomized-batch input path. -/
theorem batch_validated_input_addition_safe [Field K] (b : K)
    (point : BatchRnsPoint K)
    (validated : point.semanticallyOnCurve b)
    (canonicalInfinity : point.fieldPoint.z = 0 → point.exactZero = true) :
    point.additionSafe b := by
  rcases validated with zZero | onCurve
  · exact Or.inl ⟨canonicalInfinity zZero, zZero⟩
  · by_cases zZero : point.fieldPoint.z = 0
    · exact Or.inl ⟨canonicalInfinity zZero, zZero⟩
    · exact Or.inr ⟨onCurve, Or.inr (Or.inr zZero)⟩

theorem batch_scalar_validated_input_addition_safe [Field K] (b : K)
    (point : BatchRnsPoint K)
    (validated : point.scalarConvertedOnCurve b)
    (canonicalInfinity : point.fieldPoint.z = 0 → point.exactZero = true) :
    point.additionSafe b := by
  exact batch_validated_input_addition_safe b point
    ((batch_scalar_curve_iff_semantic_curve b point).mp validated)
    canonicalInfinity

/-- Field-level meaning of the two exact-zero selections in `G1Curve::dadd`
and `G2Curve::dadd`. The raw formula is evaluated first in C++, but the
result is replaced according to these branches. -/
def cppBatchDadd [CommRing K] (b : K)
    (left right : BatchRnsPoint K) (rawExactZero : Bool) : BatchRnsPoint K :=
  if left.exactZero = true then right
  else if right.exactZero = true then left
  else ⟨cppCompleteAdd b left.fieldPoint right.fieldPoint, rawExactZero⟩

theorem cpp_batch_dadd_addition_safe [Field K] [DecidableEq K]
    (b : K) (noTwoTorsion : NoAffineTwoTorsion b)
    (twoNonzero : (2 : K) ≠ 0)
    (left right : BatchRnsPoint K) (rawExactZero : Bool)
    (leftSafe : left.additionSafe b) (rightSafe : right.additionSafe b) :
    (cppBatchDadd b left right rawExactZero).additionSafe b := by
  by_cases leftExact : left.exactZero = true
  · simpa [cppBatchDadd, leftExact] using rightSafe
  by_cases rightExact : right.exactZero = true
  · simpa [cppBatchDadd, leftExact, rightExact] using leftSafe
  have leftValid : left.fieldPoint.projectivelyValid b := by
    rcases leftSafe with ⟨exact, _⟩ | valid
    · exact False.elim (leftExact exact)
    · exact valid
  have rightValid : right.fieldPoint.projectivelyValid b := by
    rcases rightSafe with ⟨exact, _⟩ | valid
    · exact False.elim (rightExact exact)
    · exact valid
  simp only [cppBatchDadd, leftExact, rightExact,
    BatchRnsPoint.additionSafe]
  exact Or.inr (cpp_complete_add_projectively_valid b noTwoTorsion
    twoNonzero left.fieldPoint right.fieldPoint leftValid rightValid)

theorem cpp_batch_dadd_exact_zero_sound [Field K] (b : K)
    (left right : BatchRnsPoint K) (rawExactZero : Bool)
    (leftSound : left.exactZeroSound) (rightSound : right.exactZeroSound)
    (rawSound : rawExactZero = true →
      (cppCompleteAdd b left.fieldPoint right.fieldPoint).z = 0) :
    (cppBatchDadd b left right rawExactZero).exactZeroSound := by
  by_cases leftExact : left.exactZero = true
  · simpa [cppBatchDadd, leftExact] using rightSound
  by_cases rightExact : right.exactZero = true
  · simpa [cppBatchDadd, leftExact, rightExact] using leftSound
  simpa [cppBatchDadd, leftExact, rightExact,
    BatchRnsPoint.exactZeroSound] using rawSound

/-- A point represents zero through the exact-zero branch, or a projectively
valid field point through the ordinary branch. -/
def BatchRnsPoint.represents [Field K] (b : K)
    [elliptic : (shortWeierstrass b).IsElliptic]
    (point : BatchRnsPoint K)
    (value : (shortWeierstrass b).toProjective.Point) : Prop :=
  (point.exactZero = true ∧ point.fieldPoint.z = 0 ∧ value = 0) ∨
    ∃ valid : point.fieldPoint.projectivelyValid b,
      toProjectivePoint b point.fieldPoint valid = value

theorem batch_exact_zero_represents_zero [Field K] (b : K)
    [elliptic : (shortWeierstrass b).IsElliptic]
    (point : BatchRnsPoint K)
    (value : (shortWeierstrass b).toProjective.Point)
    (exact : point.exactZero = true) (sound : point.exactZeroSound)
    (represents : point.represents b value) : value = 0 := by
  rcases represents with ⟨_, _, valueZero⟩ | ⟨valid, valueEq⟩
  · exact valueZero
  · exact valueEq.symm.trans
      (to_projective_point_eq_zero_of_z_eq_zero b point.fieldPoint
        valid (sound exact))

theorem cpp_batch_dadd_represents_sum [Field K] [DecidableEq K]
    (b : K) [elliptic : (shortWeierstrass b).IsElliptic]
    (noTwoTorsion : NoAffineTwoTorsion b)
    (twoNonzero : (2 : K) ≠ 0)
    (left right : BatchRnsPoint K) (rawExactZero : Bool)
    (leftValue rightValue : (shortWeierstrass b).toProjective.Point)
    (leftSound : left.exactZeroSound) (rightSound : right.exactZeroSound)
    (leftRepresents : left.represents b leftValue)
    (rightRepresents : right.represents b rightValue) :
    (cppBatchDadd b left right rawExactZero).represents b
      (leftValue + rightValue) := by
  by_cases leftExact : left.exactZero = true
  · have leftZero := batch_exact_zero_represents_zero b left leftValue
      leftExact leftSound leftRepresents
    simpa [cppBatchDadd, leftExact, leftZero] using rightRepresents
  by_cases rightExact : right.exactZero = true
  · have rightZero := batch_exact_zero_represents_zero b right rightValue
      rightExact rightSound rightRepresents
    simpa [cppBatchDadd, leftExact, rightExact, rightZero] using leftRepresents
  rcases leftRepresents with ⟨leftExact', _, _⟩ | ⟨leftValid, leftEq⟩
  · exact False.elim (leftExact leftExact')
  rcases rightRepresents with ⟨rightExact', _, _⟩ | ⟨rightValid, rightEq⟩
  · exact False.elim (rightExact rightExact')
  have rawEq : cppBatchDadd b left right rawExactZero =
      ⟨cppCompleteAdd b left.fieldPoint right.fieldPoint, rawExactZero⟩ := by
    simp [cppBatchDadd, leftExact, rightExact]
  rw [rawEq]
  change (rawExactZero = true ∧
      (cppCompleteAdd b left.fieldPoint right.fieldPoint).z = 0 ∧
      leftValue + rightValue = 0) ∨
    ∃ valid : (cppCompleteAdd b left.fieldPoint right.fieldPoint).projectivelyValid b,
      toProjectivePoint b (cppCompleteAdd b left.fieldPoint right.fieldPoint)
        valid = leftValue + rightValue
  right
  refine ⟨cpp_complete_add_projectively_valid b noTwoTorsion twoNonzero
    left.fieldPoint right.fieldPoint leftValid rightValid, ?_⟩
  rw [projective_complete_add_eq_add b noTwoTorsion twoNonzero
    left.fieldPoint right.fieldPoint leftValid rightValid,
    leftEq, rightEq]

def cppBatchNegate [Neg K] (point : BatchRnsPoint K) : BatchRnsPoint K :=
  ⟨point.fieldPoint.negate, point.exactZero⟩

theorem cpp_batch_negate_addition_safe [Field K] (b : K)
    (point : BatchRnsPoint K) (safe : point.additionSafe b) :
    (cppBatchNegate point).additionSafe b := by
  rcases safe with exact | valid
  · exact Or.inl (by simpa [cppBatchNegate, HomogeneousPoint.negate] using exact)
  · right
    constructor
    · exact cpp_negate_preserves_curve b point.fieldPoint valid.1
    · rcases valid.2 with xNonzero | yNonzero | zNonzero
      · exact Or.inl (by simpa [cppBatchNegate, HomogeneousPoint.negate] using xNonzero)
      · exact Or.inr (Or.inl (by simpa [cppBatchNegate, HomogeneousPoint.negate]
          using yNonzero))
      · exact Or.inr (Or.inr (by simpa [cppBatchNegate, HomogeneousPoint.negate]
          using zNonzero))

theorem cpp_batch_negate_exact_zero_sound [Field K]
    (point : BatchRnsPoint K) (sound : point.exactZeroSound) :
    (cppBatchNegate point).exactZeroSound := by
  simpa [cppBatchNegate, BatchRnsPoint.exactZeroSound,
    HomogeneousPoint.negate] using sound

/-- Model the physical exact-zero result with an arbitrary oracle. The caller
must establish that an exact-zero RNS result represents field zero. -/
def cppBatchFold [CommRing K] (b : K)
    (rawExactZero : BatchRnsPoint K → BatchRnsPoint K → Bool)
    (start : BatchRnsPoint K) : List (BatchRnsPoint K) → BatchRnsPoint K
  | [] => start
  | point :: rest =>
      cppBatchFold b rawExactZero
        (cppBatchDadd b start point (rawExactZero start point)) rest

theorem cpp_batch_fold_invariant [Field K] [DecidableEq K]
    (b : K) (noTwoTorsion : NoAffineTwoTorsion b)
    (twoNonzero : (2 : K) ≠ 0)
    (rawExactZero : BatchRnsPoint K → BatchRnsPoint K → Bool)
    (rawSound : ∀ left right, rawExactZero left right = true →
      (cppCompleteAdd b left.fieldPoint right.fieldPoint).z = 0)
    (start : BatchRnsPoint K) (points : List (BatchRnsPoint K))
    (startSafe : start.additionSafe b) (startSound : start.exactZeroSound)
    (pointsSafe : ∀ point ∈ points,
      point.additionSafe b ∧ point.exactZeroSound) :
    (cppBatchFold b rawExactZero start points).additionSafe b ∧
      (cppBatchFold b rawExactZero start points).exactZeroSound := by
  induction points generalizing start with
  | nil => exact ⟨startSafe, startSound⟩
  | cons point rest ih =>
      have pointInvariant := pointsSafe point (by simp)
      have restInvariant : ∀ q ∈ rest, q.additionSafe b ∧ q.exactZeroSound := by
        intro q membership
        exact pointsSafe q (by simp [membership])
      simpa [cppBatchFold] using
        ih (cppBatchDadd b start point (rawExactZero start point))
          (cpp_batch_dadd_addition_safe b noTwoTorsion twoNonzero
            start point (rawExactZero start point) startSafe pointInvariant.1)
          (cpp_batch_dadd_exact_zero_sound b start point
            (rawExactZero start point) startSound pointInvariant.2
            (rawSound start point)) restInvariant

theorem cpp_batch_fold_semantically_on_curve [Field K] [DecidableEq K]
    (b : K) (noTwoTorsion : NoAffineTwoTorsion b)
    (twoNonzero : (2 : K) ≠ 0)
    (rawExactZero : BatchRnsPoint K → BatchRnsPoint K → Bool)
    (rawSound : ∀ left right, rawExactZero left right = true →
      (cppCompleteAdd b left.fieldPoint right.fieldPoint).z = 0)
    (start : BatchRnsPoint K) (points : List (BatchRnsPoint K))
    (startSafe : start.additionSafe b) (startSound : start.exactZeroSound)
    (pointsSafe : ∀ point ∈ points,
      point.additionSafe b ∧ point.exactZeroSound) :
    (cppBatchFold b rawExactZero start points).semanticallyOnCurve b := by
  exact batch_addition_safe_semantically_on_curve b _
    (cpp_batch_fold_invariant b noTwoTorsion twoNonzero rawExactZero
      rawSound start points startSafe startSound pointsSafe).1

/-- The source-shaped positive and negative accumulators of one batch round.
`positiveInputs` and `negativeInputs` are the points selected by coefficients
`+1` and `-1`; coefficient generation is outside this curve-invariant proof. -/
def cppBatchRound [CommRing K] (b : K)
    (rawExactZero : BatchRnsPoint K → BatchRnsPoint K → Bool)
    (shift : BatchRnsPoint K)
    (positiveInputs negativeInputs : List (BatchRnsPoint K))
    (hasNegative : Bool) : BatchRnsPoint K :=
  let positive := cppBatchFold b rawExactZero shift positiveInputs
  let negative := cppBatchFold b rawExactZero
    ⟨cppProjectiveZero, true⟩ negativeInputs
  if hasNegative then
    cppBatchDadd b positive (cppBatchNegate negative)
      (rawExactZero positive (cppBatchNegate negative))
  else positive

theorem cpp_batch_round_invariant [Field K] [DecidableEq K]
    (b : K) (noTwoTorsion : NoAffineTwoTorsion b)
    (twoNonzero : (2 : K) ≠ 0)
    (rawExactZero : BatchRnsPoint K → BatchRnsPoint K → Bool)
    (rawSound : ∀ left right, rawExactZero left right = true →
      (cppCompleteAdd b left.fieldPoint right.fieldPoint).z = 0)
    (shift : BatchRnsPoint K)
    (positiveInputs negativeInputs : List (BatchRnsPoint K))
    (hasNegative : Bool)
    (shiftSafe : shift.additionSafe b) (shiftSound : shift.exactZeroSound)
    (positiveSafe : ∀ point ∈ positiveInputs,
      point.additionSafe b ∧ point.exactZeroSound)
    (negativeSafe : ∀ point ∈ negativeInputs,
      point.additionSafe b ∧ point.exactZeroSound) :
    (cppBatchRound b rawExactZero shift positiveInputs negativeInputs
      hasNegative).additionSafe b ∧
    (cppBatchRound b rawExactZero shift positiveInputs negativeInputs
      hasNegative).exactZeroSound := by
  let positive := cppBatchFold b rawExactZero shift positiveInputs
  let negative := cppBatchFold b rawExactZero
    ⟨cppProjectiveZero, true⟩ negativeInputs
  have positiveInvariant : positive.additionSafe b ∧
      positive.exactZeroSound :=
    cpp_batch_fold_invariant b noTwoTorsion twoNonzero rawExactZero
      rawSound shift positiveInputs shiftSafe shiftSound positiveSafe
  have negativeStartSafe :
      (⟨cppProjectiveZero, true⟩ : BatchRnsPoint K).additionSafe b := by
    exact Or.inl ⟨rfl, rfl⟩
  have negativeStartSound :
      (⟨cppProjectiveZero, true⟩ : BatchRnsPoint K).exactZeroSound := by
    intro _
    rfl
  have negativeInvariant : negative.additionSafe b ∧
      negative.exactZeroSound :=
    cpp_batch_fold_invariant b noTwoTorsion twoNonzero rawExactZero
      rawSound _ negativeInputs negativeStartSafe negativeStartSound negativeSafe
  by_cases h : hasNegative = true
  · have negatedSafe := cpp_batch_negate_addition_safe b negative
      negativeInvariant.1
    have negatedSound := cpp_batch_negate_exact_zero_sound negative
      negativeInvariant.2
    simpa [cppBatchRound, h, positive, negative] using
      And.intro
        (cpp_batch_dadd_addition_safe b noTwoTorsion twoNonzero
          positive (cppBatchNegate negative)
          (rawExactZero positive (cppBatchNegate negative))
          positiveInvariant.1 negatedSafe)
        (cpp_batch_dadd_exact_zero_sound b positive
          (cppBatchNegate negative)
          (rawExactZero positive (cppBatchNegate negative))
          positiveInvariant.2 negatedSound
          (rawSound positive (cppBatchNegate negative)))
  · simpa [cppBatchRound, h, positive, negative] using positiveInvariant

theorem cpp_batch_round_semantically_on_curve [Field K] [DecidableEq K]
    (b : K) (noTwoTorsion : NoAffineTwoTorsion b)
    (twoNonzero : (2 : K) ≠ 0)
    (rawExactZero : BatchRnsPoint K → BatchRnsPoint K → Bool)
    (rawSound : ∀ left right, rawExactZero left right = true →
      (cppCompleteAdd b left.fieldPoint right.fieldPoint).z = 0)
    (shift : BatchRnsPoint K)
    (positiveInputs negativeInputs : List (BatchRnsPoint K))
    (hasNegative : Bool)
    (shiftSafe : shift.additionSafe b) (shiftSound : shift.exactZeroSound)
    (positiveSafe : ∀ point ∈ positiveInputs,
      point.additionSafe b ∧ point.exactZeroSound)
    (negativeSafe : ∀ point ∈ negativeInputs,
      point.additionSafe b ∧ point.exactZeroSound) :
    (cppBatchRound b rawExactZero shift positiveInputs negativeInputs
      hasNegative).semanticallyOnCurve b := by
  exact batch_addition_safe_semantically_on_curve b _
    (cpp_batch_round_invariant b noTwoTorsion twoNonzero rawExactZero
      rawSound shift positiveInputs negativeInputs hasNegative shiftSafe
      shiftSound positiveSafe negativeSafe).1

/-- Curve validation is redundant on a round result satisfying the invariant;
the remaining relation predicate is deliberately arbitrary. -/
theorem cpp_batch_round_full_check_iff_relation [Field K] [DecidableEq K]
    (b : K) (noTwoTorsion : NoAffineTwoTorsion b)
    (twoNonzero : (2 : K) ≠ 0)
    (rawExactZero : BatchRnsPoint K → BatchRnsPoint K → Bool)
    (rawSound : ∀ left right, rawExactZero left right = true →
      (cppCompleteAdd b left.fieldPoint right.fieldPoint).z = 0)
    (shift : BatchRnsPoint K)
    (positiveInputs negativeInputs : List (BatchRnsPoint K))
    (hasNegative : Bool)
    (shiftSafe : shift.additionSafe b) (shiftSound : shift.exactZeroSound)
    (positiveSafe : ∀ point ∈ positiveInputs,
      point.additionSafe b ∧ point.exactZeroSound)
    (negativeSafe : ∀ point ∈ negativeInputs,
      point.additionSafe b ∧ point.exactZeroSound)
    (relation : BatchRnsPoint K → Prop) :
    ((cppBatchRound b rawExactZero shift positiveInputs negativeInputs
      hasNegative).scalarConvertedOnCurve b ∧
      relation (cppBatchRound b rawExactZero shift positiveInputs
        negativeInputs hasNegative)) ↔
      relation (cppBatchRound b rawExactZero shift positiveInputs
        negativeInputs hasNegative) := by
  have roundCurve := cpp_batch_round_semantically_on_curve b noTwoTorsion
    twoNonzero rawExactZero rawSound shift positiveInputs negativeInputs
    hasNegative shiftSafe shiftSound positiveSafe negativeSafe
  have scalarCurve := (batch_scalar_curve_iff_semantic_curve b
    (cppBatchRound b rawExactZero shift positiveInputs negativeInputs
      hasNegative)).mpr roundCurve
  simpa using and_iff_right scalarCurve

/- The prevalidation result alone is insufficient when semantic infinity is
not an exact-zero RNS representation. The production path detects this case
and uses deterministic individual checks. -/
example : ∃ point : BatchRnsPoint (ZMod 7),
    point.scalarConvertedOnCurve 4 ∧ ¬ point.additionSafe 4 := by
  refine ⟨⟨⟨1, 2, 0⟩, false⟩, ?_, ?_⟩
  · exact Or.inl rfl
  · simp [BatchRnsPoint.additionSafe, HomogeneousPoint.projectivelyValid,
      HomogeneousPoint.onCurve]
    intro impossible
    have zeroNeOne : (0 : ZMod 7) ≠ 1 := by decide
    exact False.elim (zeroNeOne impossible)

theorem selected_g1_batch_round_full_check_iff_relation
    (rawExactZero : BatchRnsPoint SelectedCurveFp →
      BatchRnsPoint SelectedCurveFp → Bool)
    (rawSound : ∀ left right, rawExactZero left right = true →
      (cppCompleteAdd 4 left.fieldPoint right.fieldPoint).z = 0)
    (shift : BatchRnsPoint SelectedCurveFp)
    (positiveInputs negativeInputs : List (BatchRnsPoint SelectedCurveFp))
    (hasNegative : Bool)
    (shiftSafe : shift.additionSafe 4) (shiftSound : shift.exactZeroSound)
    (positiveSafe : ∀ point ∈ positiveInputs,
      point.additionSafe 4 ∧ point.exactZeroSound)
    (negativeSafe : ∀ point ∈ negativeInputs,
      point.additionSafe 4 ∧ point.exactZeroSound)
    (relation : BatchRnsPoint SelectedCurveFp → Prop) :
    ((cppBatchRound 4 rawExactZero shift positiveInputs negativeInputs
      hasNegative).scalarConvertedOnCurve 4 ∧
      relation (cppBatchRound 4 rawExactZero shift positiveInputs
        negativeInputs hasNegative)) ↔
      relation (cppBatchRound 4 rawExactZero shift positiveInputs
        negativeInputs hasNegative) := by
  exact cpp_batch_round_full_check_iff_relation 4
    selected_g1_no_affine_two_torsion selected_g1_two_ne_zero
    rawExactZero rawSound shift positiveInputs negativeInputs hasNegative
    shiftSafe shiftSound positiveSafe negativeSafe relation

theorem selected_g2_batch_round_full_check_iff_relation
    (rawExactZero : BatchRnsPoint SelectedMillerFp2 →
      BatchRnsPoint SelectedMillerFp2 → Bool)
    (rawSound : ∀ left right, rawExactZero left right = true →
      (cppCompleteAdd selectedG2CurveB left.fieldPoint right.fieldPoint).z = 0)
    (shift : BatchRnsPoint SelectedMillerFp2)
    (positiveInputs negativeInputs : List (BatchRnsPoint SelectedMillerFp2))
    (hasNegative : Bool)
    (shiftSafe : shift.additionSafe selectedG2CurveB)
    (shiftSound : shift.exactZeroSound)
    (positiveSafe : ∀ point ∈ positiveInputs,
      point.additionSafe selectedG2CurveB ∧ point.exactZeroSound)
    (negativeSafe : ∀ point ∈ negativeInputs,
      point.additionSafe selectedG2CurveB ∧ point.exactZeroSound)
    (relation : BatchRnsPoint SelectedMillerFp2 → Prop) :
    ((cppBatchRound selectedG2CurveB rawExactZero shift positiveInputs
      negativeInputs hasNegative).scalarConvertedOnCurve selectedG2CurveB ∧
      relation (cppBatchRound selectedG2CurveB rawExactZero shift
        positiveInputs negativeInputs hasNegative)) ↔
      relation (cppBatchRound selectedG2CurveB rawExactZero shift
        positiveInputs negativeInputs hasNegative) := by
  exact cpp_batch_round_full_check_iff_relation selectedG2CurveB
    selected_g2_no_affine_two_torsion selected_g2_two_ne_zero
    rawExactZero rawSound shift positiveInputs negativeInputs hasNegative
    shiftSafe shiftSound positiveSafe negativeSafe relation

/-!
The friendly direct path does not use AVX addition. Native AVX `in_subgroup`
is `on_curve && g1_in_subgroup_relation` or
`on_curve && g2_in_subgroup_relation` in `point.h`. The direct path calls the
same relation after validation. This Boolean theorem covers the entire input
list, including an early false from either predicate. The separate
semantic-infinity branch is a representation-level boundary. -/

def cppFullIndividualChecks {P : Type*}
    (onCurve relation : P → Bool) (points : List P) : Bool :=
  points.all (fun point => onCurve point && relation point)

def cppPrevalidatedRelationChecks {P : Type*}
    (onCurve relation : P → Bool) (points : List P) : Bool :=
  points.all onCurve && points.all relation

theorem cpp_prevalidated_relation_checks_eq_full {P : Type*}
    (onCurve relation : P → Bool) (points : List P) :
    cppPrevalidatedRelationChecks onCurve relation points =
      cppFullIndividualChecks onCurve relation points := by
  induction points with
  | nil => simp [cppPrevalidatedRelationChecks, cppFullIndividualChecks]
  | cons point rest ih =>
      simp only [cppPrevalidatedRelationChecks, cppFullIndividualChecks,
        List.all_cons] at ih ⊢
      rw [← ih]
      simp [Bool.and_assoc, Bool.and_left_comm]

theorem cpp_relation_only_after_prevalidation {P : Type*}
    (onCurve relation : P → Bool) (points : List P)
    (validated : points.all onCurve = true) :
    points.all relation = cppFullIndividualChecks onCurve relation points := by
  rw [← cpp_prevalidated_relation_checks_eq_full]
  simp [cppPrevalidatedRelationChecks, validated]

/-! The G2 coefficient alphabet in `avx_c_api.cpp` is 0, ±1, ..., ±6.
The quotient by the prime-order subgroup has the concrete G2 cofactor as its
order. That cofactor is coprime to every nonzero difference between two
coefficients, so at most one choice can cancel a fixed nonzero quotient
component in each row. -/

theorem g1_cofactor_divisible_by_three : 3 ∣ g1Cofactor := by
  native_decide

theorem g2_cofactor_coprime_small (d : Nat) (lower : 1 ≤ d)
    (upper : d ≤ 12) : g2Cofactor.Coprime d := by
  interval_cases d <;> native_decide

/-- The quotient-cardinality premise used below follows from the ambient
curve order and the already-proved exact order of the selected generator. -/
theorem selected_g2_quotient_card
    (fullOrder : Nat.card selectedG2Weierstrass.toProjective.Point =
      g2Cofactor * blsScalarOrder) :
    Nat.card (selectedG2Weierstrass.toProjective.Point ⧸
      AddSubgroup.zmultiples selectedG2GeneratorPoint) = g2Cofactor := by
  have product := AddSubgroup.card_eq_card_quotient_mul_card_addSubgroup
    (AddSubgroup.zmultiples selectedG2GeneratorPoint)
  rw [fullOrder, selected_g2_generator_span_natCard] at product
  exact mul_right_cancel₀ (by native_decide : blsScalarOrder ≠ 0)
    product.symm

def cppThirteenShiftedCoefficient (i : Fin 13) : Nat :=
  match i.val with
  | 0 => 6
  | 1 => 7
  | 2 => 5
  | 3 => 8
  | 4 => 4
  | 5 => 9
  | 6 => 3
  | 7 => 10
  | 8 => 2
  | 9 => 11
  | 10 => 1
  | 11 => 12
  | _ => 0

theorem cpp_thirteen_shifted_coefficient_injective :
    Function.Injective cppThirteenShiftedCoefficient := by
  decide

def cppThirteenCoefficient {G : Type*} [AddCommGroup G]
    (x : G) (i : Fin 13) : G :=
  match i.val with
  | 0 => 0
  | 1 => x
  | 2 => -x
  | 3 => 2 • x
  | 4 => -(2 • x)
  | 5 => 3 • x
  | 6 => -(3 • x)
  | 7 => 4 • x
  | 8 => -(4 • x)
  | 9 => 5 • x
  | 10 => -(5 • x)
  | 11 => 6 • x
  | _ => -(6 • x)

theorem cpp_thirteen_coefficient_shift {G : Type*} [AddCommGroup G]
    (x : G) (i : Fin 13) :
    cppThirteenCoefficient x i + 6 • x =
      cppThirteenShiftedCoefficient i • x := by
  fin_cases i <;>
    simp [cppThirteenCoefficient, cppThirteenShiftedCoefficient,
      succ_nsmul]
  all_goals abel

theorem cpp_thirteen_coefficients_distinct {G : Type*}
    [AddCommGroup G] [Finite G]
    (card : Nat.card G = g2Cofactor) (x : G) (nonzero : x ≠ 0) :
    Function.Injective (cppThirteenCoefficient x) := by
  intro i j equal
  have shiftedEqual : cppThirteenShiftedCoefficient i • x =
      cppThirteenShiftedCoefficient j • x := by
    rw [← cpp_thirteen_coefficient_shift,
      ← cpp_thirteen_coefficient_shift, equal]
  have distinctOrEq := lt_trichotomy
    (cppThirteenShiftedCoefficient i)
    (cppThirteenShiftedCoefficient j)
  rcases distinctOrEq with smaller | same | larger
  · let d := cppThirteenShiftedCoefficient j -
      cppThirteenShiftedCoefficient i
    have dLower : 1 ≤ d := by dsimp [d]; omega
    have dUpper : d ≤ 12 := by
      dsimp [d]
      have hi : cppThirteenShiftedCoefficient i ≤ 12 := by
        fin_cases i <;> decide
      have hj : cppThirteenShiftedCoefficient j ≤ 12 := by
        fin_cases j <;> decide
      omega
    have zero : d • x = 0 := by
      rw [show d = cppThirteenShiftedCoefficient j -
        cppThirteenShiftedCoefficient i from rfl,
        sub_nsmul x smaller.le, shiftedEqual, add_neg_cancel]
    have coprime : (Nat.card G).Coprime d := by
      rw [card]
      exact g2_cofactor_coprime_small d dLower dUpper
    have xZero : x = 0 :=
      (Nat.Coprime.nsmul_right_bijective coprime).injective
        (by simpa using zero)
    exact False.elim (nonzero xZero)
  · exact cpp_thirteen_shifted_coefficient_injective same
  · let d := cppThirteenShiftedCoefficient i -
      cppThirteenShiftedCoefficient j
    have dLower : 1 ≤ d := by dsimp [d]; omega
    have dUpper : d ≤ 12 := by
      dsimp [d]
      have hi : cppThirteenShiftedCoefficient i ≤ 12 := by
        fin_cases i <;> decide
      have hj : cppThirteenShiftedCoefficient j ≤ 12 := by
        fin_cases j <;> decide
      omega
    have zero : d • x = 0 := by
      rw [show d = cppThirteenShiftedCoefficient i -
        cppThirteenShiftedCoefficient j from rfl,
        sub_nsmul x larger.le, ← shiftedEqual, add_neg_cancel]
    have coprime : (Nat.card G).Coprime d := by
      rw [card]
      exact g2_cofactor_coprime_small d dLower dUpper
    have xZero : x = 0 :=
      (Nat.Coprime.nsmul_right_bijective coprime).injective
        (by simpa using zero)
    exact False.elim (nonzero xZero)

theorem cpp_thirteen_way_round_at_most_one {G : Type*}
    [AddCommGroup G] [Finite G] [DecidableEq G]
    (card : Nat.card G = g2Cofactor) (x : G) (nonzero : x ≠ 0)
    (offset : G) :
    (Finset.univ.filter (fun i : Fin 13 =>
      offset + cppThirteenCoefficient x i = 0)).card ≤ 1 := by
  classical
  apply Finset.card_le_one_iff.mpr
  intro i j iMember jMember
  have iZero := (Finset.mem_filter.mp iMember).2
  have jZero := (Finset.mem_filter.mp jMember).2
  apply cpp_thirteen_coefficients_distinct card x nonzero
  exact add_left_cancel (iZero.trans jZero.symm)

theorem cpp_thirteen_way_conditional_rows_at_most_one_any_rounds {G : Type*}
    [AddCommGroup G] [Finite G] [DecidableEq G]
    (card : Nat.card G = g2Cofactor) (x : G) (nonzero : x ≠ 0)
    (rounds : Nat) (offset : Fin rounds → G) :
    Set.Subsingleton {choices : Fin rounds → Fin 13 |
      ∀ round, offset round +
        cppThirteenCoefficient x (choices round) = 0} := by
  intro first firstMember second secondMember
  apply funext
  intro round
  have firstZero := firstMember round
  have secondZero := secondMember round
  apply cpp_thirteen_coefficients_distinct card x nonzero
  exact add_left_cancel (firstZero.trans secondZero.symm)

theorem cpp_thirteen_way_conditional_rows_at_most_one {G : Type*}
    [AddCommGroup G] [Finite G] [DecidableEq G]
    (card : Nat.card G = g2Cofactor) (x : G) (nonzero : x ≠ 0)
    (offset : Fin 35 → G) :
    Set.Subsingleton {choices : Fin 35 → Fin 13 |
      ∀ round, offset round +
        cppThirteenCoefficient x (choices round) = 0} := by
  exact cpp_thirteen_way_conditional_rows_at_most_one_any_rounds
    card x nonzero 35 offset

theorem cpp_thirteen_way_challenge_space_card :
    Fintype.card (Fin 35 → Fin 13) = 13 ^ 35 := by
  simp

theorem g2_thirteen_way_round_bound : 13 ^ 35 > 2 ^ 128 := by
  native_decide

theorem cpp_thirteen_way_22_rows_at_most_one {G : Type*}
    [AddCommGroup G] [Finite G] [DecidableEq G]
    (card : Nat.card G = g2Cofactor) (x : G) (nonzero : x ≠ 0)
    (offset : Fin 22 → G) :
    Set.Subsingleton {choices : Fin 22 → Fin 13 |
      ∀ round, offset round +
        cppThirteenCoefficient x (choices round) = 0} := by
  exact cpp_thirteen_way_conditional_rows_at_most_one_any_rounds
    card x nonzero 22 offset

/-- Instantiate the 22-row uniqueness bound for the selected G2 quotient.
The ambient curve-order identity is the remaining cardinality premise. -/
theorem selected_g2_thirteen_way_22_rows_at_most_one
    (fullOrder : Nat.card selectedG2Weierstrass.toProjective.Point =
      g2Cofactor * blsScalarOrder)
    (x : selectedG2Weierstrass.toProjective.Point ⧸
      AddSubgroup.zmultiples selectedG2GeneratorPoint)
    (nonzero : x ≠ 0)
    (offset : Fin 22 →
      selectedG2Weierstrass.toProjective.Point ⧸
        AddSubgroup.zmultiples selectedG2GeneratorPoint) :
    Set.Subsingleton {choices : Fin 22 → Fin 13 |
      ∀ round, offset round +
        cppThirteenCoefficient x (choices round) = 0} := by
  classical
  letI : Finite selectedG2Weierstrass.toProjective.Point :=
    Nat.finite_of_card_ne_zero (by rw [fullOrder]; native_decide)
  exact cpp_thirteen_way_22_rows_at_most_one
    (selected_g2_quotient_card fullOrder) x nonzero offset

theorem cpp_thirteen_way_22_challenge_space_card :
    Fintype.card (Fin 22 → Fin 13) = 13 ^ 22 := by
  simp

theorem g2_thirteen_way_22_round_bound : 13 ^ 22 > 2 ^ 80 := by
  native_decide

/-! The large-table fallback uses the first nine coefficients of the
thirteen-value alphabet, so its uniqueness proof follows by restriction. -/

def cppNineCoefficient {G : Type*} [AddCommGroup G]
    (x : G) (i : Fin 9) : G :=
  cppThirteenCoefficient x
    ⟨i.val, Nat.lt_trans i.isLt (by decide : 9 < 13)⟩

theorem cpp_nine_coefficients_distinct {G : Type*}
    [AddCommGroup G] [Finite G]
    (card : Nat.card G = g2Cofactor) (x : G) (nonzero : x ≠ 0) :
    Function.Injective (cppNineCoefficient x) := by
  intro i j equal
  have lifted := cpp_thirteen_coefficients_distinct card x nonzero equal
  exact Fin.ext (congrArg (fun k : Fin 13 => k.val) lifted)

theorem cpp_nine_way_conditional_rows_at_most_one {G : Type*}
    [AddCommGroup G] [Finite G] [DecidableEq G]
    (card : Nat.card G = g2Cofactor) (x : G) (nonzero : x ≠ 0)
    (offset : Fin 41 → G) :
    Set.Subsingleton {choices : Fin 41 → Fin 9 |
      ∀ round, offset round +
        cppNineCoefficient x (choices round) = 0} := by
  intro first firstMember second secondMember
  apply funext
  intro round
  have firstZero := firstMember round
  have secondZero := secondMember round
  apply cpp_nine_coefficients_distinct card x nonzero
  exact add_left_cancel (firstZero.trans secondZero.symm)

theorem cpp_nine_way_challenge_space_card :
    Fintype.card (Fin 41 → Fin 9) = 9 ^ 41 := by
  simp

theorem g2_nine_way_round_bound : 9 ^ 41 > 2 ^ 128 := by
  native_decide

end CatBLST.AVX.MatrixNoK
