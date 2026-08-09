/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.CurveComplete
import CatBLST.AVX.MatrixNoK.Fp12FrobeniusConstants
import Mathlib.Tactic

set_option maxHeartbeats 1600000

/-!
# Concrete G1 GLV and G2 GLS endomorphisms

This file transcribes the coordinate maps in `GLVMult::batch_endomorphism`
and `GLSMult::minus_psi`.  It proves that their literal constants have the
required field identities, that both maps preserve the selected curves, and
that their signs agree with BLST's subgroup-relation predicates.

The distinction between the signed BLS parameter and its positive magnitude
is intentional.  Scalar BLST proves `sigma^2(P) = [-|z|^2]P` for G1 and
`psi(P) = [-|z|]P` for G2.  The AVX maps negate those source endomorphisms, so
the unsigned split bases are `[|z|^2]P` and `[|z|]P`, respectively.
-/

namespace CatBLST.AVX.MatrixNoK

def scalarG1SigmaBeta : SelectedCurveFp :=
  0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaac

def avxG1Beta : SelectedCurveFp :=
  0x5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a00022e01fffffffefffe

theorem scalar_g1_sigma_beta_cube : scalarG1SigmaBeta ^ 3 = 1 := by
  native_decide

theorem avx_g1_beta_cube : avxG1Beta ^ 3 = 1 := by
  native_decide

theorem avx_g1_beta_is_scalar_sigma_beta_squared :
    avxG1Beta = scalarG1SigmaBeta ^ 2 := by
  native_decide

/-- Homogeneous form of scalar BLST's `sigma` map. -/
def cppG1Sigma (point : HomogeneousPoint SelectedCurveFp) :
    HomogeneousPoint SelectedCurveFp :=
  ⟨scalarG1SigmaBeta * point.x, point.y, point.z⟩

/-- Exact reduced-field map in `GLVMult::batch_endomorphism`. -/
def cppG1AvxEndomorphism (point : HomogeneousPoint SelectedCurveFp) :
    HomogeneousPoint SelectedCurveFp :=
  ⟨avxG1Beta * point.x, -point.y, point.z⟩

theorem cpp_g1_sigma_cube (point : HomogeneousPoint SelectedCurveFp) :
    cppG1Sigma (cppG1Sigma (cppG1Sigma point)) = point := by
  apply HomogeneousPoint.ext
  · change scalarG1SigmaBeta *
        (scalarG1SigmaBeta * (scalarG1SigmaBeta * point.x)) = point.x
    calc
      _ = scalarG1SigmaBeta ^ 3 * point.x := by ring
      _ = point.x := by rw [scalar_g1_sigma_beta_cube]; simp
  · rfl
  · rfl

theorem cpp_g1_avx_endomorphism_eq_neg_sigma_squared
    (point : HomogeneousPoint SelectedCurveFp) :
    cppG1AvxEndomorphism point = (cppG1Sigma (cppG1Sigma point)).negate := by
  ext
  all_goals
    simp [cppG1AvxEndomorphism, cppG1Sigma, HomogeneousPoint.negate,
      avx_g1_beta_is_scalar_sigma_beta_squared, pow_two] <;> ring

theorem cpp_g1_avx_endomorphism_preserves_curve
    (b : SelectedCurveFp) (point : HomogeneousPoint SelectedCurveFp)
    (onCurve : point.onCurve b) :
    (cppG1AvxEndomorphism point).onCurve b := by
  simp only [HomogeneousPoint.onCurve] at onCurve ⊢
  simp only [cppG1AvxEndomorphism]
  rw [mul_pow, avx_g1_beta_cube]
  simpa using onCurve

def scaleXPoint [Mul K] (beta : K) (point : HomogeneousPoint K) :
    HomogeneousPoint K :=
  ⟨beta * point.x, point.y, point.z⟩

theorem complete_add_scale_x [CommRing K]
    (b beta : K) (betaCube : beta ^ 3 = 1)
    (left right : HomogeneousPoint K) :
    scaleXPoint beta (cppCompleteAdd b left right) =
      cppCompleteAdd b (scaleXPoint beta left) (scaleXPoint beta right) := by
  rw [cpp_complete_add_coordinates, cpp_complete_add_coordinates]
  apply HomogeneousPoint.ext
  · simp only [scaleXPoint]
    ring
  · simp only [scaleXPoint]
    linear_combination
      -(3 * left.x * right.x *
        (3 * b * (left.x * right.z + left.z * right.x))) * betaCube
  · simp only [scaleXPoint]
    linear_combination
      -(3 * left.x * right.x *
        (left.x * right.y + left.y * right.x)) * betaCube

theorem cpp_g1_sigma_complete_add
    (b : SelectedCurveFp)
    (left right : HomogeneousPoint SelectedCurveFp) :
    cppG1Sigma (cppCompleteAdd b left right) =
      cppCompleteAdd b (cppG1Sigma left) (cppG1Sigma right) := by
  change scaleXPoint scalarG1SigmaBeta (cppCompleteAdd b left right) =
    cppCompleteAdd b (scaleXPoint scalarG1SigmaBeta left)
      (scaleXPoint scalarG1SigmaBeta right)
  exact complete_add_scale_x b scalarG1SigmaBeta scalar_g1_sigma_beta_cube
    left right

def mapHomogeneousPoint [CommRing K] (phi : K →+* K)
    (point : HomogeneousPoint K) : HomogeneousPoint K :=
  ⟨phi point.x, phi point.y, phi point.z⟩

theorem complete_add_map [CommRing K]
    (phi : K →+* K) (b : K) (left right : HomogeneousPoint K) :
    mapHomogeneousPoint phi (cppCompleteAdd b left right) =
      cppCompleteAdd (phi b) (mapHomogeneousPoint phi left)
        (mapHomogeneousPoint phi right) := by
  rw [cpp_complete_add_coordinates, cpp_complete_add_coordinates]
  ext <;>
    simp only [mapHomogeneousPoint, map_add, map_sub, map_mul, map_ofNat]

def scaleXYPoint [Mul K] (cx cy : K) (point : HomogeneousPoint K) :
    HomogeneousPoint K :=
  ⟨cx * point.x, cy * point.y, point.z⟩

def homogeneousScale [Mul K] (scale : K) (point : HomogeneousPoint K) :
    HomogeneousPoint K :=
  ⟨scale * point.x, scale * point.y, scale * point.z⟩

theorem complete_add_scale_xy [CommRing K]
    (sourceB targetB cx cy : K)
    (coefficientPowers : cy ^ 2 = cx ^ 3)
    (curveCoefficient : targetB = cy ^ 2 * sourceB)
    (left right : HomogeneousPoint K) :
    cppCompleteAdd targetB (scaleXYPoint cx cy left)
        (scaleXYPoint cx cy right) =
      homogeneousScale (cy ^ 3)
        (scaleXYPoint cx cy (cppCompleteAdd sourceB left right)) := by
  subst targetB
  rw [cpp_complete_add_coordinates, cpp_complete_add_coordinates]
  apply HomogeneousPoint.ext
  · simp only [scaleXYPoint, homogeneousScale]
    ring
  · simp only [scaleXYPoint, homogeneousScale]
    ring_nf
    rw [← coefficientPowers]
    ring
  · simp only [scaleXYPoint, homogeneousScale]
    ring_nf
    rw [← coefficientPowers]
    ring

theorem homogeneous_scale_equivalent [CommRing K]
    (factor : K) (point : HomogeneousPoint K) (_factorNonzero : factor ≠ 0) :
    (homogeneousScale factor point).equivalent point := by
  simp only [HomogeneousPoint.equivalent, homogeneousScale]
  by_cases atInfinity : point.z = 0
  · left
    exact ⟨by simp [atInfinity], atInfinity⟩
  · right
    constructor <;> ring

theorem homogeneous_equivalent_negate [CommRing K]
    (left right : HomogeneousPoint K)
    (equivalent : left.equivalent right) :
    left.negate.equivalent right.negate := by
  rcases equivalent with atInfinity | finite
  · exact Or.inl atInfinity
  · right
    constructor
    · exact finite.1
    · simpa [HomogeneousPoint.negate] using
        congrArg (fun value => -value) finite.2

theorem homogeneous_negate_negate [Ring K] (point : HomogeneousPoint K) :
    point.negate.negate = point := by
  ext <;> simp [HomogeneousPoint.negate]

/-!
This is the exact sign bridge needed by the unsigned G1 split.  The premise is
the relation tested by scalar BLST, expressed against the point `[|z|^2]P`.
No subgroup membership is inferred from curve membership alone.
-/
theorem g1_avx_endomorphism_unsigned_basis
    (point radixSquaredPoint : HomogeneousPoint SelectedCurveFp)
    (subgroupRelation :
      (cppG1Sigma (cppG1Sigma point)).equivalent
        radixSquaredPoint.negate) :
    (cppG1AvxEndomorphism point).equivalent radixSquaredPoint := by
  rw [cpp_g1_avx_endomorphism_eq_neg_sigma_squared]
  have negated := homogeneous_equivalent_negate _ _ subgroupRelation
  simpa [homogeneous_negate_negate] using negated

def fp2Conjugate (value : SelectedMillerFp2) : SelectedMillerFp2 :=
  ⟨value.x, -value.y⟩

@[simp] theorem selected_fp2_mul_x (left right : SelectedMillerFp2) :
    (left * right).x = left.x * right.x - left.y * right.y := rfl

@[simp] theorem selected_fp2_mul_y (left right : SelectedMillerFp2) :
    (left * right).y = left.x * right.y + left.y * right.x := rfl

@[simp] theorem selected_fp2_neg_x (value : SelectedMillerFp2) :
    (-value).x = -value.x := rfl

@[simp] theorem selected_fp2_neg_y (value : SelectedMillerFp2) :
    (-value).y = -value.y := rfl

@[simp] theorem selected_fp2_zero_x :
    (0 : SelectedMillerFp2).x = 0 := rfl

theorem fp2_conjugate_add (left right : SelectedMillerFp2) :
    fp2Conjugate (left + right) = fp2Conjugate left + fp2Conjugate right := by
  apply Fp2Value.ext
  · rfl
  · change -(left.y + right.y) = -left.y + -right.y
    ring

theorem fp2_conjugate_mul (left right : SelectedMillerFp2) :
    fp2Conjugate (left * right) = fp2Conjugate left * fp2Conjugate right := by
  apply Fp2Value.ext
  · change left.x * right.x - left.y * right.y =
      left.x * right.x - (-left.y) * (-right.y)
    ring
  · change -(left.x * right.y + left.y * right.x) =
      left.x * (-right.y) + (-left.y) * right.x
    ring

@[simp] theorem fp2_conjugate_conjugate (value : SelectedMillerFp2) :
    fp2Conjugate (fp2Conjugate value) = value := by
  apply Fp2Value.ext
  · rfl
  · change - -value.y = value.y
    ring

theorem fp2_conjugate_pow (value : SelectedMillerFp2) (exponent : Nat) :
    fp2Conjugate (value ^ exponent) = fp2Conjugate value ^ exponent := by
  induction exponent with
  | zero =>
      apply Fp2Value.ext
      · rfl
      · change -(0 : SelectedCurveFp) = 0
        simp
  | succ exponent ih =>
      rw [pow_succ, fp2_conjugate_mul, ih, pow_succ]

def fp2ConjugateHom : SelectedMillerFp2 →+* SelectedMillerFp2 where
  toFun := fp2Conjugate
  map_zero' := by
    apply Fp2Value.ext
    · rfl
    · change -(0 : SelectedCurveFp) = 0
      simp
  map_one' := by
    apply Fp2Value.ext
    · rfl
    · change -(0 : SelectedCurveFp) = 0
      simp
  map_add' := fp2_conjugate_add
  map_mul' := fp2_conjugate_mul

def avxG2PsiXImaginary : SelectedCurveFp :=
  0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaad

def avxG2PsiYReal : SelectedCurveFp :=
  0x135203e60180a68ee2e9c448d77a2cd91c3dedd930b1cf60ef396489f61eb45e304466cf3e67fa0af1ee7b04121bdea2

def avxG2PsiYImaginary : SelectedCurveFp :=
  0x06af0e0437ff400b6831e36d6bd17ffe48395dabc2d3435e77f76e17009241c5ee67992f72ec05f4c81084fbede3cc09

def g2PsiXCoefficient : SelectedMillerFp2 :=
  ⟨0, avxG2PsiXImaginary⟩

def g2PsiYCoefficient : SelectedMillerFp2 :=
  ⟨avxG2PsiYReal, avxG2PsiYImaginary⟩

def selectedG2CurveB : SelectedMillerFp2 := ⟨4, 4⟩

def g2PsiSquaredXCoefficient : SelectedMillerFp2 :=
  g2PsiXCoefficient * fp2Conjugate g2PsiXCoefficient

def g2PsiSquaredYCoefficient : SelectedMillerFp2 :=
  g2PsiYCoefficient * fp2Conjugate g2PsiYCoefficient

theorem avx_g2_psi_y_imaginary_is_negative_real :
    avxG2PsiYImaginary = -avxG2PsiYReal := by
  native_decide

theorem g2_psi_coefficient_square_cube :
    g2PsiYCoefficient ^ 2 = g2PsiXCoefficient ^ 3 := by
  native_decide

theorem g2_psi_curve_coefficient_compatibility :
    g2PsiYCoefficient ^ 2 * fp2Conjugate selectedG2CurveB =
      selectedG2CurveB := by
  native_decide

theorem g2_psi_x_cube_curve_coefficient_compatibility :
    g2PsiXCoefficient ^ 3 * fp2Conjugate selectedG2CurveB =
      selectedG2CurveB := by
  rw [← g2_psi_coefficient_square_cube]
  exact g2_psi_curve_coefficient_compatibility

theorem g2_psi_squared_x_coefficient_cube :
    g2PsiSquaredXCoefficient ^ 3 = 1 := by
  native_decide

theorem g2_psi_squared_y_coefficient_cube :
    g2PsiSquaredYCoefficient ^ 3 = -1 := by
  native_decide

theorem g2_psi_y_coefficient_nonzero : g2PsiYCoefficient ≠ 0 := by
  intro coefficientZero
  have realZero : avxG2PsiYReal = 0 := by
    simpa [g2PsiYCoefficient] using
      congrArg Fp2Value.x coefficientZero
  have realNonzero : avxG2PsiYReal ≠ 0 := by native_decide
  exact realNonzero realZero

theorem g2_psi_y_coefficient_cube_nonzero :
    g2PsiYCoefficient ^ 3 ≠ 0 := by
  intro cubeZero
  have cubeRealZero : (g2PsiYCoefficient ^ 3).x = 0 := by
    simpa using congrArg Fp2Value.x cubeZero
  have cubeRealNonzero : (g2PsiYCoefficient ^ 3).x ≠ 0 := by
    native_decide
  exact cubeRealNonzero cubeRealZero

/-- Scalar BLST's `psi`: conjugate each coordinate, then scale X and Y. -/
def cppG2Psi (point : HomogeneousPoint SelectedMillerFp2) :
    HomogeneousPoint SelectedMillerFp2 :=
  ⟨g2PsiXCoefficient * fp2Conjugate point.x,
    g2PsiYCoefficient * fp2Conjugate point.y,
    fp2Conjugate point.z⟩

theorem cpp_g2_psi_complete_add_scaled
    (left right : HomogeneousPoint SelectedMillerFp2) :
    cppCompleteAdd selectedG2CurveB (cppG2Psi left) (cppG2Psi right) =
      homogeneousScale (g2PsiYCoefficient ^ 3)
        (cppG2Psi (cppCompleteAdd selectedG2CurveB left right)) := by
  change cppCompleteAdd selectedG2CurveB
      (scaleXYPoint g2PsiXCoefficient g2PsiYCoefficient
        (mapHomogeneousPoint fp2ConjugateHom left))
      (scaleXYPoint g2PsiXCoefficient g2PsiYCoefficient
        (mapHomogeneousPoint fp2ConjugateHom right)) =
    homogeneousScale (g2PsiYCoefficient ^ 3)
      (scaleXYPoint g2PsiXCoefficient g2PsiYCoefficient
        (mapHomogeneousPoint fp2ConjugateHom
          (cppCompleteAdd selectedG2CurveB left right)))
  rw [complete_add_scale_xy
    (fp2ConjugateHom selectedG2CurveB) selectedG2CurveB
    g2PsiXCoefficient g2PsiYCoefficient
    g2_psi_coefficient_square_cube
    g2_psi_curve_coefficient_compatibility.symm]
  rw [complete_add_map]

theorem cpp_g2_psi_complete_add
    (left right : HomogeneousPoint SelectedMillerFp2) :
    (cppCompleteAdd selectedG2CurveB (cppG2Psi left)
      (cppG2Psi right)).equivalent
      (cppG2Psi (cppCompleteAdd selectedG2CurveB left right)) := by
  rw [cpp_g2_psi_complete_add_scaled]
  exact homogeneous_scale_equivalent _ _
    g2_psi_y_coefficient_cube_nonzero

def cppG2PsiSquaredScale (point : HomogeneousPoint SelectedMillerFp2) :
    HomogeneousPoint SelectedMillerFp2 :=
  ⟨g2PsiSquaredXCoefficient * point.x,
    g2PsiSquaredYCoefficient * point.y,
    point.z⟩

theorem cpp_g2_psi_squared (point : HomogeneousPoint SelectedMillerFp2) :
    cppG2Psi (cppG2Psi point) = cppG2PsiSquaredScale point := by
  ext <;>
    simp only [cppG2Psi, cppG2PsiSquaredScale, fp2_conjugate_mul,
      fp2_conjugate_conjugate, g2PsiSquaredXCoefficient,
      g2PsiSquaredYCoefficient] <;>
    ring_nf

theorem cpp_g2_psi_six (point : HomogeneousPoint SelectedMillerFp2) :
    cppG2Psi (cppG2Psi (cppG2Psi (cppG2Psi (cppG2Psi (cppG2Psi point))))) =
      point.negate := by
  rw [cpp_g2_psi_squared, cpp_g2_psi_squared, cpp_g2_psi_squared]
  apply HomogeneousPoint.ext
  · change g2PsiSquaredXCoefficient * (g2PsiSquaredXCoefficient *
        (g2PsiSquaredXCoefficient * point.x)) = point.x
    calc
      _ = g2PsiSquaredXCoefficient ^ 3 * point.x := by ring
      _ = point.x := by rw [g2_psi_squared_x_coefficient_cube]; simp
  · change g2PsiSquaredYCoefficient * (g2PsiSquaredYCoefficient *
        (g2PsiSquaredYCoefficient * point.y)) = -point.y
    calc
      _ = g2PsiSquaredYCoefficient ^ 3 * point.y := by ring
      _ = -point.y := by rw [g2_psi_squared_y_coefficient_cube]; simp
  · rfl

/-- Exact component schedule in `GLSMult::minus_psi`. -/
def cppG2AvxMinusPsi (point : HomogeneousPoint SelectedMillerFp2) :
    HomogeneousPoint SelectedMillerFp2 :=
  ⟨⟨point.x.y * avxG2PsiXImaginary,
      point.x.x * avxG2PsiXImaginary⟩,
    ⟨point.y.y * avxG2PsiYReal + point.y.x * avxG2PsiYImaginary,
      point.y.y * avxG2PsiYReal + point.y.x * avxG2PsiYReal⟩,
    ⟨point.z.x, -point.z.y⟩⟩

theorem cpp_g2_avx_minus_psi_eq_neg_psi
    (point : HomogeneousPoint SelectedMillerFp2) :
    cppG2AvxMinusPsi point = (cppG2Psi point).negate := by
  unfold cppG2AvxMinusPsi cppG2Psi
  rw [avx_g2_psi_y_imaginary_is_negative_real]
  simp only [HomogeneousPoint.negate]
  apply HomogeneousPoint.ext <;> apply Fp2Value.ext <;>
    simp only [g2PsiXCoefficient, g2PsiYCoefficient, fp2Conjugate,
      selected_fp2_mul_x, selected_fp2_mul_y,
      selected_fp2_neg_x, selected_fp2_neg_y,
      avx_g2_psi_y_imaginary_is_negative_real] <;>
    ring

theorem cpp_g2_psi_preserves_curve
    (point : HomogeneousPoint SelectedMillerFp2)
    (onCurve : point.onCurve selectedG2CurveB) :
    (cppG2Psi point).onCurve selectedG2CurveB := by
  simp only [HomogeneousPoint.onCurve] at onCurve ⊢
  have conjugated := congrArg fp2Conjugate onCurve
  rw [fp2_conjugate_mul, fp2_conjugate_add,
    fp2_conjugate_pow, fp2_conjugate_mul, fp2_conjugate_pow,
    fp2_conjugate_pow] at conjugated
  simp only [cppG2Psi]
  calc
    (g2PsiYCoefficient * fp2Conjugate point.y) ^ 2 *
          fp2Conjugate point.z =
        g2PsiYCoefficient ^ 2 *
          (fp2Conjugate point.y ^ 2 * fp2Conjugate point.z) := by ring
    _ = g2PsiYCoefficient ^ 2 *
          (fp2Conjugate point.x ^ 3 +
            fp2Conjugate selectedG2CurveB * fp2Conjugate point.z ^ 3) := by
        rw [conjugated]
    _ = g2PsiYCoefficient ^ 2 * fp2Conjugate point.x ^ 3 +
          (g2PsiYCoefficient ^ 2 * fp2Conjugate selectedG2CurveB) *
            fp2Conjugate point.z ^ 3 := by ring
    _ = g2PsiXCoefficient ^ 3 * fp2Conjugate point.x ^ 3 +
          selectedG2CurveB * fp2Conjugate point.z ^ 3 := by
        rw [g2_psi_coefficient_square_cube,
          g2_psi_x_cube_curve_coefficient_compatibility]
    _ = (g2PsiXCoefficient * fp2Conjugate point.x) ^ 3 +
          selectedG2CurveB * fp2Conjugate point.z ^ 3 := by ring

theorem cpp_g2_avx_minus_psi_preserves_curve
    (point : HomogeneousPoint SelectedMillerFp2)
    (onCurve : point.onCurve selectedG2CurveB) :
    (cppG2AvxMinusPsi point).onCurve selectedG2CurveB := by
  rw [cpp_g2_avx_minus_psi_eq_neg_psi]
  exact cpp_negate_preserves_curve selectedG2CurveB _
    (cpp_g2_psi_preserves_curve point onCurve)

/-!
The scalar subgroup test uses the signed BLS parameter `-|z|`, hence its
relation is `psi(P) = [-|z|]P`.  Negating the source `psi` yields the positive
radix basis selected by the AVX unsigned decomposition.
-/
theorem g2_avx_minus_psi_unsigned_basis
    (point radixPoint : HomogeneousPoint SelectedMillerFp2)
    (subgroupRelation :
      (cppG2Psi point).equivalent radixPoint.negate) :
    (cppG2AvxMinusPsi point).equivalent radixPoint := by
  rw [cpp_g2_avx_minus_psi_eq_neg_psi]
  have negated := homogeneous_equivalent_negate _ _ subgroupRelation
  simpa [homogeneous_negate_negate] using negated

end CatBLST.AVX.MatrixNoK
