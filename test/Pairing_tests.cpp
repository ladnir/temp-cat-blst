#include "gtest/gtest.h"

extern "C" {
#include "src/point.h"
#include "src/fields.h"
#include "src/e1.h"
#include "src/e2.h"
#include "src/pairing.h"
#include "src/fp12_tower.h"
}


TEST(PairingTest, BasicPairing) {
    // Test basic pairing functionality
    POINTonE1 g1_point = *blst_p1_generator();
    POINTonE2 g2_point = *blst_p2_generator();

    // Create affine versions of the points
    POINTonE1_affine g1_affine;
    POINTonE2_affine g2_affine;
    blst_p1_to_affine(&g1_affine, &g1_point);
    blst_p2_to_affine(&g2_affine, &g2_point);

    // Compute pairing e(G1, G2)
    vec384fp12 pairing_result;
    blst_miller_loop(pairing_result, &g2_affine, &g1_affine);
    blst_final_exp(pairing_result, pairing_result);

    // The pairing should not be equal to 1 (non-degeneracy)
    EXPECT_FALSE(blst_fp12_is_one(pairing_result))
        << "Pairing should be non-degenerate";

    // Test bilinearity: e(aP, Q) = e(P, aQ) = e(P,Q)^a
    // Let's use a = 2 for simplicity
    POINTonE1 g1_double;
    POINTonE2 g2_double;
    blst_p1_double(&g1_double, &g1_point);
    blst_p2_double(&g2_double, &g2_point);

    // Convert to affine
    POINTonE1_affine g1_double_affine;
    POINTonE2_affine g2_double_affine;
    blst_p1_to_affine(&g1_double_affine, &g1_double);
    blst_p2_to_affine(&g2_double_affine, &g2_double);

    // Compute e(2*G1, G2)
    vec384fp12 pairing_2g1_g2;
    blst_miller_loop(pairing_2g1_g2, &g2_affine, &g1_double_affine);
    blst_final_exp(pairing_2g1_g2, pairing_2g1_g2);

    // Compute e(G1, 2*G2)
    vec384fp12 pairing_g1_2g2;
    blst_miller_loop(pairing_g1_2g2, &g2_double_affine, &g1_affine);
    blst_final_exp(pairing_g1_2g2, pairing_g1_2g2);

    // Compute e(G1, G2)^2
    vec384fp12 pairing_squared;
    blst_fp12_sqr(pairing_squared, pairing_result);

    // All three results should be equal
    EXPECT_TRUE(blst_fp12_is_equal(pairing_2g1_g2, pairing_g1_2g2))
        << "e(2*G1, G2) != e(G1, 2*G2)";

    EXPECT_TRUE(blst_fp12_is_equal(pairing_2g1_g2, pairing_squared))
        << "e(2*G1, G2) != e(G1, G2)^2";

    // Test pairing with identity points
    POINTonE1 identity_e1;
    POINTonE2 identity_e2;
    memset(&identity_e1, 0, sizeof(identity_e1));
    memset(&identity_e2, 0, sizeof(identity_e2));

    POINTonE1_affine identity_e1_affine;
    POINTonE2_affine identity_e2_affine;
    blst_p1_to_affine(&identity_e1_affine, &identity_e1);
    blst_p2_to_affine(&identity_e2_affine, &identity_e2);

    // e(G1, O) should equal 1
    vec384fp12 pairing_with_identity;
    blst_miller_loop(pairing_with_identity, &identity_e2_affine, &g1_affine);
    blst_final_exp(pairing_with_identity, pairing_with_identity);

    EXPECT_TRUE(blst_fp12_is_one(pairing_with_identity))
        << "e(G1, O) should equal 1";

    // e(O, G2) should equal 1
    blst_miller_loop(pairing_with_identity, &g2_affine, &identity_e1_affine);
    blst_final_exp(pairing_with_identity, pairing_with_identity);

    EXPECT_TRUE(blst_fp12_is_one(pairing_with_identity))
        << "e(O, G2) should equal 1";
}

TEST(PairingTest, MultiPairing) {
    // Test multi-pairing functionality

    // Create some sample points
    POINTonE1 g1 = *blst_p1_generator();
    POINTonE2 g2 = *blst_p2_generator();

    POINTonE1 g1_double;
    POINTonE2 g2_double;
    blst_p1_double(&g1_double, &g1);
    blst_p2_double(&g2_double, &g2);

    // Convert to affine
    POINTonE1_affine g1_affine, g1_double_affine;
    POINTonE2_affine g2_affine, g2_double_affine;
    blst_p1_to_affine(&g1_affine, &g1);
    blst_p1_to_affine(&g1_double_affine, &g1_double);
    blst_p2_to_affine(&g2_affine, &g2);
    blst_p2_to_affine(&g2_double_affine, &g2_double);

    // Test that e(P1,Q1)·e(P2,Q2) = e(P1+P2,Q1+Q2)·e(P1,Q2)⁻¹·e(P2,Q1)⁻¹

    // Compute e(G1, G2)
    vec384fp12 e_g1_g2;
    blst_miller_loop(e_g1_g2, &g2_affine, &g1_affine);
    blst_final_exp(e_g1_g2, e_g1_g2);

    // Compute e(2*G1, 2*G2)
    vec384fp12 e_2g1_2g2;
    blst_miller_loop(e_2g1_2g2, &g2_double_affine, &g1_double_affine);
    blst_final_exp(e_2g1_2g2, e_2g1_2g2);

    // Compute e(G1, G2) * e(2*G1, 2*G2) using vec384fp12_mul
    vec384fp12 product;
    blst_fp12_mul(product, e_g1_g2, e_2g1_2g2);

    // Test multi-pairing API if available
    POINTonE1_affine* g1_points[2] = { &g1_affine, &g1_double_affine };
    POINTonE2_affine* g2_points[2] = { &g2_affine, &g2_double_affine };

    vec384fp12 multi_result;
    blst_miller_loop_n(multi_result, g2_points, g1_points, 2);
    blst_final_exp(multi_result, multi_result);

    // The results should be equal
    EXPECT_TRUE(blst_fp12_is_equal(product, multi_result))
        << "Multi-pairing result does not match individual pairings";
}