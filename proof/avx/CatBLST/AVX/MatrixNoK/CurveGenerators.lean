/-
Copyright 2026 Category Labs, Inc.
SPDX-License-Identifier: Apache-2.0
-/

import CatBLST.AVX.MatrixNoK.CurveBinaryMul
import CatBLST.AVX.MatrixNoK.CurveEndomorphism

/-!
# Concrete BLS12-381 G1/G2 generators

The canonical coordinates below are the standard BLS12-381 generators used by
the scalar implementation and by the independent MatrixNoK constant
regenerator.  The proofs use the logarithmic `cppBinaryMul` evaluator, whose
correspondence to Mathlib projective scalar multiplication is established in
`CurveBinaryMul.lean`.
-/

namespace CatBLST.AVX.MatrixNoK

def selectedG1Generator : HomogeneousPoint SelectedCurveFp :=
  ⟨
    0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb,
    0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1,
    1
  ⟩

def selectedG2Generator : HomogeneousPoint SelectedMillerFp2 :=
  ⟨
    ⟨
      0x024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8,
      0x13e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e
    ⟩,
    ⟨
      0x0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160d12c923ac9cc3baca289e193548608b82801,
      0x0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e99ab3f370d275cec1da1aaa9075ff05f79be
    ⟩,
    1
  ⟩

theorem selected_g1_generator_valid :
    selectedG1Generator.projectivelyValid 4 := by
  constructor
  · change selectedG1Generator.y ^ 2 * selectedG1Generator.z =
      selectedG1Generator.x ^ 3 + 4 * selectedG1Generator.z ^ 3
    native_decide
  · right
    right
    norm_num [selectedG1Generator]

theorem selected_g2_generator_valid :
    selectedG2Generator.projectivelyValid selectedG2CurveB := by
  constructor
  · change selectedG2Generator.y ^ 2 * selectedG2Generator.z =
      selectedG2Generator.x ^ 3 +
        selectedG2CurveB * selectedG2Generator.z ^ 3
    native_decide
  · right
    right
    norm_num [selectedG2Generator]

/-- The concrete G1 generator is killed by the BLS12-381 scalar order. -/
theorem selected_g1_generator_order_divides_r :
    (cppBinaryMul 4 selectedG1Generator blsScalarOrder).equivalent
      cppProjectiveZero := by
  left
  constructor <;> native_decide

/-- The concrete G2 generator is killed by the BLS12-381 scalar order. -/
theorem selected_g2_generator_order_divides_r :
    (cppBinaryMul selectedG2CurveB selectedG2Generator blsScalarOrder).equivalent
      cppProjectiveZero := by
  left
  constructor <;> native_decide

/-- The selected AVX G1 endomorphism is `[z^2]` on the concrete generator. -/
theorem selected_g1_generator_avx_basis :
    (cppG1AvxEndomorphism selectedG1Generator).equivalent
      (cppBinaryMul 4 selectedG1Generator curveRadixSquared) := by
  right
  constructor <;> native_decide

/-- The selected AVX G2 `minus_psi` map is `[z]` on the concrete generator. -/
theorem selected_g2_generator_avx_basis :
    (cppG2AvxMinusPsi selectedG2Generator).equivalent
      (cppBinaryMul selectedG2CurveB selectedG2Generator curveRadix) := by
  right
  constructor <;> native_decide

noncomputable def selectedG1GeneratorPoint :
    selectedG1Weierstrass.toProjective.Point :=
  toSelectedG1Point selectedG1Generator selected_g1_generator_valid

noncomputable def selectedG2GeneratorPoint :
    selectedG2Weierstrass.toProjective.Point :=
  toSelectedG2Point selectedG2Generator selected_g2_generator_valid

theorem selected_g1_generator_point_order_divides_r :
    blsScalarOrder • selectedG1GeneratorPoint = 0 := by
  unfold selectedG1GeneratorPoint
  rw [← to_selected_g1_point_cpp_binary_mul
    selectedG1Generator selected_g1_generator_valid blsScalarOrder]
  apply (to_selected_g1_point_eq_iff _ _
    (cpp_binary_mul_projectively_valid 4 selected_g1_no_affine_two_torsion
      selected_g1_two_ne_zero selectedG1Generator selected_g1_generator_valid
      blsScalarOrder)
    selected_g1_zero_valid).2
  exact selected_g1_generator_order_divides_r

theorem selected_g2_generator_point_order_divides_r :
    blsScalarOrder • selectedG2GeneratorPoint = 0 := by
  unfold selectedG2GeneratorPoint
  rw [← to_selected_g2_point_cpp_binary_mul
    selectedG2Generator selected_g2_generator_valid blsScalarOrder]
  apply (to_selected_g2_point_eq_iff _ _
    (cpp_binary_mul_projectively_valid selectedG2CurveB
      selected_g2_no_affine_two_torsion selected_g2_two_ne_zero
      selectedG2Generator selected_g2_generator_valid blsScalarOrder)
    selected_g2_zero_valid).2
  exact selected_g2_generator_order_divides_r

theorem selected_g1_generator_avx_basis_group :
    toSelectedG1Point (cppG1AvxEndomorphism selectedG1Generator)
        (by
          constructor
          · exact cpp_g1_avx_endomorphism_preserves_curve 4 selectedG1Generator
              selected_g1_generator_valid.1
          · native_decide) =
      curveRadixSquared • selectedG1GeneratorPoint := by
  unfold selectedG1GeneratorPoint
  rw [← to_selected_g1_point_cpp_binary_mul
    selectedG1Generator selected_g1_generator_valid curveRadixSquared]
  apply (to_selected_g1_point_eq_iff _ _ _ _).2
  exact selected_g1_generator_avx_basis

theorem selected_g2_generator_avx_basis_group :
    toSelectedG2Point (cppG2AvxMinusPsi selectedG2Generator)
        (by
          constructor
          · exact cpp_g2_avx_minus_psi_preserves_curve selectedG2Generator
              selected_g2_generator_valid.1
          · native_decide) =
      curveRadix • selectedG2GeneratorPoint := by
  unfold selectedG2GeneratorPoint
  rw [← to_selected_g2_point_cpp_binary_mul
    selectedG2Generator selected_g2_generator_valid curveRadix]
  apply (to_selected_g2_point_eq_iff _ _ _ _).2
  exact selected_g2_generator_avx_basis

end CatBLST.AVX.MatrixNoK
