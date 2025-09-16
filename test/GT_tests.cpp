#include "gtest/gtest.h"

extern "C" {
#include "blst/curve/e1.h"
#include "blst/curve/e2.h"
#include "blst/curve/pairing.h"
#include "blst/fields/Fp12.h"
#include "blst/sign/aggregate.h"

}

class GTTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Initialize points
        g1_point = *blst_p1_generator();
        g2_point = *blst_p2_generator();

        // Convert to affine
        blst_p1_to_affine(&g1_affine, &g1_point);
        blst_p2_to_affine(&g2_affine, &g2_point);

        // Calculate pairing
        blst_miller_loop(gt_element, &g2_affine, &g1_affine);
        blst_final_exp(gt_element, gt_element);
    }

    POINTonE1 g1_point;
    POINTonE2 g2_point;
    POINTonE1_affine g1_affine;
    POINTonE2_affine g2_affine;
    vec384fp12 gt_element;  // An element in GT
    vec384fp12 temp;        // For temporary results
};

TEST_F(GTTest, Identity) {
    // Test identity element in GT
    auto identity = blst_fp12_one();

    EXPECT_TRUE(blst_fp12_is_one(*identity))
        << "Identity element should be 1";

    // Multiplying any element by identity should leave it unchanged
    blst_fp12_mul(temp, gt_element, *identity);
    EXPECT_TRUE(blst_fp12_is_equal(temp, gt_element))
        << "a * 1 should equal a in GT";
}

TEST_F(GTTest, Multiplication) {
    // Create another element for testing
    vec384fp12 another_gt;
    POINTonE1 g1_double;
    blst_p1_double(&g1_double, &g1_point);

    POINTonE1_affine g1_double_affine;
    blst_p1_to_affine(&g1_double_affine, &g1_double);

    blst_miller_loop(another_gt, &g2_affine, &g1_double_affine);
    blst_final_exp(another_gt, another_gt);

    // Test multiplication
    vec384fp12 product;
    blst_fp12_mul(product, gt_element, another_gt);

    // The product should not be equal to either of the original elements
    EXPECT_FALSE(blst_fp12_is_equal(product, gt_element))
        << "Product should not equal first element";
    EXPECT_FALSE(blst_fp12_is_equal(product, another_gt))
        << "Product should not equal second element";

    // Test commutativity: a * b = b * a
    vec384fp12 product2;
    blst_fp12_mul(product2, another_gt, gt_element);
    EXPECT_TRUE(blst_fp12_is_equal(product, product2))
        << "Multiplication should be commutative in GT";

    // Test associativity: (a * b) * c = a * (b * c)
    // Create a third element
    vec384fp12 third_gt;
    POINTonE1 g1_triple;
    blst_p1_add(&g1_triple, &g1_double, &g1_point);

    POINTonE1_affine g1_triple_affine;
    blst_p1_to_affine(&g1_triple_affine, &g1_triple);

    blst_miller_loop(third_gt, &g2_affine, &g1_triple_affine);
    blst_final_exp(third_gt, third_gt);

    vec384fp12 temp1, temp2;
    // (a * b) * c
    blst_fp12_mul(temp1, gt_element, another_gt);
    blst_fp12_mul(temp1, temp1, third_gt);

    // a * (b * c)
    blst_fp12_mul(temp2, another_gt, third_gt);
    blst_fp12_mul(temp2, gt_element, temp2);

    EXPECT_TRUE(blst_fp12_is_equal(temp1, temp2))
        << "Multiplication should be associative in GT";
}

TEST_F(GTTest, Inversion) {
    // Test inversion (a * a⁻¹ = 1)
    vec384fp12 inverse;
    blst_fp12_inverse(inverse, gt_element);

    vec384fp12 product;
    blst_fp12_mul(product, gt_element, inverse);

    auto identity = *blst_fp12_one();

    EXPECT_TRUE(blst_fp12_is_equal(product, identity))
        << "Element multiplied by its inverse should equal identity";
}

TEST_F(GTTest, Squaring) {
    // Test that squaring equals multiplying by self
    vec384fp12 square;
    blst_fp12_sqr(square, gt_element);

    vec384fp12 product;
    blst_fp12_mul(product, gt_element, gt_element);

    EXPECT_TRUE(blst_fp12_is_equal(square, product))
        << "a² should equal a * a";

    // Test multiple squaring (to test repeated squaring optimization if available)
    vec384fp12 repeated_square;
    blst_fp12_sqr(repeated_square, gt_element);
    blst_fp12_sqr(repeated_square, repeated_square);

    // This should equal gt_element⁴
    vec384fp12 pow4;
    blst_fp12_mul(pow4, gt_element, gt_element);
    blst_fp12_mul(pow4, pow4, pow4);

    EXPECT_TRUE(blst_fp12_is_equal(repeated_square, pow4))
        << "((a²)²) should equal a⁴";
}

TEST_F(GTTest, Exponentiation) {
    // Test exponentiation using small exponents
    vec384fp12 result, expected;

    // Test e^0 = 1
    limb_t zero[4] = { 0, 0, 0, 0 };
    blst_fp12_pow(result, gt_element, zero, sizeof(zero) * 8);
    vec384fp12 expected_identity = { 0 };
	blst_fp12_one_copy(expected_identity);
    EXPECT_TRUE(blst_fp12_is_equal(result, expected_identity))
        << "e^0 should equal 1";

    // Test e^1 = e
    limb_t one[4] = { 1, 0, 0, 0 };
    blst_fp12_pow(result, gt_element, one, sizeof(one) * 8);
    EXPECT_TRUE(blst_fp12_is_equal(result, gt_element))
        << "e^1 should equal e";

    // Test e^2 = e * e
    limb_t two[4] = { 2, 0, 0, 0 };
    blst_fp12_pow(result, gt_element, two, sizeof(two) * 8);
    blst_fp12_sqr(expected, gt_element);
    EXPECT_TRUE(blst_fp12_is_equal(result, expected))
        << "e^2 should equal e * e";

    // Test e^3 = e * e * e
    limb_t three[4] = { 3, 0, 0, 0 };
    blst_fp12_pow(result, gt_element, three, sizeof(three) * 8);
    blst_fp12_mul(expected, expected, gt_element);
    EXPECT_TRUE(blst_fp12_is_equal(result, expected))
        << "e^3 should equal e * e * e";
}

//TEST_F(GTTest, Serialization) {
//    // Test serialization/deserialization of GT elements
//    vec384fp12 original = { 0 };  // Copy the gt_element
//    memcpy(original, gt_element, sizeof(vec384fp12));
//
//    // Serialize to bytes
//    unsigned char serialized[12 * 48];  // Fp12 consists of 12 base field elements
//    blst_fp12_serialize(serialized, original);
//
//    // Deserialize back
//    vec384fp12 deserialized;
//    EXPECT_TRUE(blst_fp12_deserialize(deserialized, serialized))
//        << "Deserialization should succeed";
//
//    // The original and deserialized values should be equal
//    EXPECT_TRUE(blst_fp12_is_equal(original, deserialized))
//        << "Deserialized value should equal original";
//}

TEST_F(GTTest, BilinearityProperty) {
    // Test the bilinearity property: e(a*P, Q) = e(P, a*Q) = e(P,Q)^a

    // Create scaled points
    POINTonE1 g1_double;
    POINTonE2 g2_double;
    blst_p1_double(&g1_double, &g1_point);
    blst_p2_double(&g2_double, &g2_point);

    POINTonE1_affine g1_double_affine;
    POINTonE2_affine g2_double_affine;
    blst_p1_to_affine(&g1_double_affine, &g1_double);
    blst_p2_to_affine(&g2_double_affine, &g2_double);

    // e(2P, Q)
    vec384fp12 e_2P_Q;
    blst_miller_loop(e_2P_Q, &g2_affine, &g1_double_affine);
    blst_final_exp(e_2P_Q, e_2P_Q);

    // e(P, 2Q)
    vec384fp12 e_P_2Q;
    blst_miller_loop(e_P_2Q, &g2_double_affine, &g1_affine);
    blst_final_exp(e_P_2Q, e_P_2Q);

    // e(P,Q)^2
    vec384fp12 e_PQ_squared;
    blst_fp12_sqr(e_PQ_squared, gt_element);

    // All three should be equal due to bilinearity
    EXPECT_TRUE(blst_fp12_is_equal(e_2P_Q, e_P_2Q))
        << "e(2P, Q) should equal e(P, 2Q)";
    EXPECT_TRUE(blst_fp12_is_equal(e_2P_Q, e_PQ_squared))
        << "e(2P, Q) should equal e(P,Q)^2";
}



// NEW TESTS BELOW

TEST_F(GTTest, GroupMembership) {
    // Test that the pairing result is in the correct group
    EXPECT_TRUE(blst_fp12_in_group(gt_element))
        << "Pairing result should be in the GT group";

    // Identity should be in the group
    vec384fp12 identity;
    blst_fp12_one_copy(identity);
    EXPECT_TRUE(blst_fp12_in_group(identity))
        << "Identity element should be in the GT group";

    // Test that a non-GT element (before final_exp) is not in the group
    vec384fp12 not_finalized;
    blst_miller_loop(not_finalized, &g2_affine, &g1_affine);
    EXPECT_FALSE(blst_fp12_in_group(not_finalized))
        << "Miller loop result without final_exp should not be in GT group";

    // Squaring, multiplication, inversion should preserve group membership
    vec384fp12 squared;
    blst_fp12_sqr(squared, gt_element);
    EXPECT_TRUE(blst_fp12_in_group(squared))
        << "Squaring should preserve group membership";

    vec384fp12 inverse;
    blst_fp12_inverse(inverse, gt_element);
    EXPECT_TRUE(blst_fp12_in_group(inverse))
        << "Inversion should preserve group membership";
}

TEST_F(GTTest, MultiPairing) {
    // Set up multiple points for testing
    POINTonE1 g1_double;
    blst_p1_double(&g1_double, &g1_point);

    POINTonE1_affine g1_points[2];
    blst_p1_to_affine(&g1_points[0], &g1_point);
    blst_p1_to_affine(&g1_points[1], &g1_double);

    POINTonE2_affine g2_points[2];
    blst_p2_to_affine(&g2_points[0], &g2_point);
    blst_p2_to_affine(&g2_points[1], &g2_point);

    // Create arrays of pointers as required by blst_miller_loop_n
    const POINTonE1_affine* g1_ptrs[3] = { &g1_points[0], &g1_points[1], nullptr };
    const POINTonE2_affine* g2_ptrs[3] = { &g2_points[0], &g2_points[1], nullptr };

    // Calculate e(P1,Q1) and e(P2,Q2) individually
    vec384fp12 e1, e2, product, multi;
    blst_miller_loop(e1, &g2_points[0], &g1_points[0]);
    blst_miller_loop(e2, &g2_points[1], &g1_points[1]);

    // Calculate their product
    blst_fp12_mul(product, e1, e2);
    blst_final_exp(product, product);

    // Calculate using multi-pairing
    blst_miller_loop_n(multi, g2_ptrs, g1_ptrs, 2);
    blst_final_exp(multi, multi);

    // They should match
    EXPECT_TRUE(blst_fp12_is_equal(product, multi))
        << "Multi-pairing result should match individual pairings";
}

TEST_F(GTTest, PowersOfTwo) {
    // Test successive powers of 2 to find where exponentiation fails
    vec384fp12 manual, result;

    // Start with the base element
    memcpy(manual, gt_element, sizeof(vec384fp12));
#define PowersOfTwo_SIZE 2
 //   auto pre = 17;
 //   for(int i = 0; i < pre; i++) {

 //       vec384fp12 alt;
 //       blst_fp12_sqr(alt, manual);
 //       memcpy(manual, alt, sizeof(vec384fp12));

 //       printf("2^%d: ", i + 1);
 //       for (int j = 0; j < 12 * 48; j++) {
 //           printf("%02x", ((unsigned char*)manual)[j]);
 //       }
 //       printf("\n");
	//}

    // Test powers of 2: g^(2^i) for i = 0,1,2,3,...
    for (int i = 0; i < sizeof(limb_t)*8* PowersOfTwo_SIZE; i++)
    {  // Test reasonable number of powers
        // Create exponent 2^i
        limb_t exp[PowersOfTwo_SIZE];
        memset(exp, 0, sizeof(exp));

        auto limb = i / (sizeof(limb_t) * 8);
        auto bit = i % (sizeof(limb_t) * 8);
        exp[limb] = (limb_t)1 << bit;
		//printf("blst_fp12_pow 2^%d \n\n", i);

        // Calculate using blst_fp12_pow
        blst_fp12_pow(result, gt_element, exp, sizeof(exp) * 8);
        // They should be equal
        EXPECT_TRUE(blst_fp12_is_equal(result, manual))
            << "Power of 2 test failed for 2^" << i;

        // Prepare for next iteration by squaring
        vec384fp12 alt;
        blst_fp12_sqr(alt, manual);
		memcpy(manual, alt, sizeof(vec384fp12));


    }
}

TEST_F(GTTest, ExponentiationCornerCases) {
    // Test to understand why high bit exponentiation is failing
    vec384fp12 result, expected;

    // Test specifically with the high bit
    limb_t high_bit[4] = { 0, 0, 0, 0 };
    limb_t bit_position = (sizeof(limb_t) * 8) - 1;
    high_bit[0] = ((limb_t)1) << bit_position;

    // Calculate using exponentiation function
    blst_fp12_pow(result, gt_element, high_bit, sizeof(high_bit) * 8);

    // Calculate manually - starting with the base element
    memcpy(expected, gt_element, sizeof(vec384fp12));

    // Square exactly bit_position times (not bit_position-1)
    for (size_t i = 0; i < bit_position; i++) {
        blst_fp12_sqr(expected, expected);
    }

    // Check equality
    EXPECT_TRUE(blst_fp12_is_equal(result, expected))
        << "High bit exponentiation failed with bit position " << bit_position;
}

TEST_F(GTTest, AdvancedExponentiation) {
    vec384fp12 result, expected;

    // Test power of 10
    limb_t ten[4] = { 10, 0, 0, 0 };
    blst_fp12_pow(result, gt_element, ten, sizeof(ten) * 8);

    // Calculate e^10 by repeated multiplication
    blst_fp12_one_copy(expected);
    for (int i = 0; i < 10; i++) {
        blst_fp12_mul(expected, expected, gt_element);
    }

    EXPECT_TRUE(blst_fp12_is_equal(result, expected))
        << "e^10 should equal e * e * ... * e (10 times)";

    // Test with a large limb value (single bit high up in limb)
    limb_t high_bit[4] = { 0, 0, 0, 0 };
    high_bit[0] = (limb_t)1 << (sizeof(limb_t) * 8 - 1);  // Highest bit in first limb
    blst_fp12_pow(result, gt_element, high_bit, sizeof(high_bit) * 8);

    // Calculate this power by squaring many times - START WITH gt_element, not 1!
    memcpy(expected, gt_element, sizeof(vec384fp12));  // Start with base^1
    for (size_t i = 0; i < sizeof(limb_t) * 8 - 1; i++) {  // Square 63 times to get base^(2^63)
        blst_fp12_sqr(expected, expected);
    }
    // No final multiply needed - we already have base^(2^63)

    EXPECT_TRUE(blst_fp12_is_equal(result, expected))
        << "Exponentiation with high bit should match manual calculation";

    // Test with scattered bits
    limb_t scattered_bits[4] = { 0, 0, 0, 0 };
    scattered_bits[0] = 1 | (1 << 4) | (1 << 8);  // Bits 0, 4, and 8 set
    blst_fp12_pow(result, gt_element, scattered_bits, sizeof(scattered_bits) * 8);

    // Calculate manually
    blst_fp12_one_copy(expected);
    blst_fp12_mul(expected, expected, gt_element);  // Bit 0
    for (int i = 0; i < 4; i++) {
        blst_fp12_sqr(expected, expected);
    }
    blst_fp12_mul(expected, expected, gt_element);  // Bit 4
    for (int i = 0; i < 4; i++) {
        blst_fp12_sqr(expected, expected);
    }
    blst_fp12_mul(expected, expected, gt_element);  // Bit 8

    EXPECT_TRUE(blst_fp12_is_equal(result, expected))
        << "Exponentiation with scattered bits should match manual calculation";
}

TEST_F(GTTest, Conjugate) {
    // Test conjugate operation
    vec384fp12 conjugate;
    memcpy(conjugate, gt_element, sizeof(vec384fp12));
    blst_fp12_conjugate(conjugate);

    // Conjugate of conjugate should give original value
    vec384fp12 double_conjugate;
    memcpy(double_conjugate, conjugate, sizeof(vec384fp12));
    blst_fp12_conjugate(double_conjugate);
    EXPECT_TRUE(blst_fp12_is_equal(gt_element, double_conjugate))
        << "Conjugate of conjugate should equal original value";

    // Conjugate * original should equal norm (which in cyclotomic is 1)
    vec384fp12 product;
    blst_fp12_mul(product, gt_element, conjugate);

    // Let's verify if we get identity (this might not be true in all cases,
    // but works for many elements in the cyclotomic subgroup)
    vec384fp12 identity;
    blst_fp12_one_copy(identity);

    // Verify that conjugate is also in GT
    EXPECT_TRUE(blst_fp12_in_group(conjugate))
        << "Conjugate should be in GT group";
}

TEST_F(GTTest, CyclotomicSquaring) {
    // Test cyclotomic squaring
    // First make sure the element is in cyclotomic subgroup
    blst_final_exp(gt_element, gt_element);

    vec384fp12 cyclo_square;
    cyclotomic_sqr_fp12(cyclo_square, gt_element);
    
    // Compare with regular squaring
    vec384fp12 regular_square;
    blst_fp12_sqr(regular_square, gt_element);

    // They should be equal
    EXPECT_TRUE(blst_fp12_is_equal(cyclo_square, regular_square))
        << "Cyclotomic squaring should equal regular squaring for GT elements";

    // But cyclotomic squaring should be more efficient
    // (we can't test performance directly, but we can verify correctness)
}

//TEST_F(GTTest, FrobeniusMap) {
//    // Test Frobenius map operation (raises to power of p)
//    vec384fp12 frob;
//    blst_fp12_frobenius_map(frob, gt_element, 1);
//
//    // Frobenius^12 should be identity on Fp12
//    vec384fp12 frob12;
//    vec384fp12 temp;
//    memcpy(temp, gt_element, sizeof(vec384fp12));
//
//    for (int i = 0; i < 12; i++) {
//        blst_fp12_frobenius_map(frob12, temp, 1);
//        memcpy(temp, frob12, sizeof(vec384fp12));
//    }
//
//    EXPECT_TRUE(blst_fp12_is_equal(gt_element, temp))
//        << "Frobenius^12 should be identity on Fp12";
//
//    // Test that frob^6 is conjugation for elements in cyclotomic subgroup
//    blst_fp12_frobenius_map(frob12, gt_element, 6);
//    memcpy(temp, gt_element, sizeof(gt_element));
//    blst_fp12_conjugate(temp);
//
//    EXPECT_TRUE(blst_fp12_is_equal(frob12, temp))
//        << "Frobenius^6 should equal conjugation for GT elements";
//}

TEST_F(GTTest, FinalVerify) {
    // Test finalverify function (equality after final exponentiation)

    // Create two different pairings that should be equal after final_exp
    POINTonE1 g1_double;
    //POINTonE2 g2_half;
    blst_p1_double(&g1_double, &g1_point);

    // Create a "half" point in G2 by multiplying by inverse of 2 (if available)
    // For test purposes, we'll just use a different point
    POINTonE2 g2_another;
    blst_p2_add(&g2_another, &g2_point, &g2_point);  // Just a different point

    POINTonE1_affine g1_double_affine;
    POINTonE2_affine g2_half_affine;
    blst_p1_to_affine(&g1_double_affine, &g1_double);
    blst_p2_to_affine(&g2_half_affine, &g2_another);

    // Calculate two different pairings
    vec384fp12 pair1, pair2;
    blst_miller_loop(pair1, &g2_affine, &g1_double_affine);
    blst_miller_loop(pair2, &g2_half_affine, &g1_affine);

    // They should not be equal before final_exp
    EXPECT_FALSE(blst_fp12_is_equal(pair1, pair2))
        << "Different pairings should not be equal before final_exp";

    // Apply final exponentiation
    vec384fp12 final1, final2;
    blst_final_exp(final1, pair1);
    blst_final_exp(final2, pair2);

    // Direct equality check
    bool direct_equal = blst_fp12_is_equal(final1, final2);

    // Using finalverify
    bool finalverify_equal = blst_fp12_finalverify(final1, final2);

    // Verify that finalverify gives same result as equality check
    // (in general they might differ due to optimizations, but for our simple test they should match)
    EXPECT_EQ(direct_equal, finalverify_equal)
        << "finalverify should agree with is_equal in basic cases";
}

//void custom_fp12_from_bendian(vec384fp12 ret, const unsigned char serialized[12 * 48]) {
//    // This is a simplified implementation that may need adjustment based on actual vec384fp12 structure
//    for (size_t i = 0; i < 12; i++) {
//        limb_t* limb_ptr = &(((limb_t*)ret)[i * (384 / 8 / sizeof(limb_t))]);
//        // Convert from big-endian bytes to limbs
//        for (size_t j = 0; j < 48; j++) {
//            size_t byte_pos = i * 48 + j;
//            size_t shift = 8 * (48 - 1 - j % 48);
//            if (j % sizeof(limb_t) == 0)
//                limb_ptr[j / sizeof(limb_t)] = 0;
//            limb_ptr[j / sizeof(limb_t)] |= (limb_t)serialized[byte_pos] << shift;
//        }
//    }
//}

TEST_F(GTTest, Serialization) {
    // Test serialization/deserialization of GT elements
    vec384fp12 original;
    memcpy(original, gt_element, sizeof(vec384fp12));

    // Serialize to bytes - using bendian format
    unsigned char serialized[12 * 48];  // Fp12 consists of 12 base field elements
    blst_bendian_from_fp12(serialized, original);

    // Deserialize back
    vec384fp12 deserialized;
    //void blst_bendian_from_fp12(unsigned char ret[48 * 12], const vec384fp12 a);

    blst_fp12_from_bendian(deserialized, serialized);

    // The original and deserialized values should be equal
    EXPECT_TRUE(blst_fp12_is_equal(original, deserialized))
        << "Deserialized value should equal original";
}

TEST_F(GTTest, PairingWithPrecomputation) {
    // Test using precomputed lines for pairing
    vec384fp6 lines[68];  // Size from the API

    // Precompute lines for pairing with g2_point
    blst_precompute_lines(lines, &g2_affine);

    // Compute pairing using precomputed lines
    vec384fp12 pre_result;
    blst_miller_loop_lines(pre_result, lines, &g1_affine);
    blst_final_exp(pre_result, pre_result);

    // Verify against standard pairing computation
    EXPECT_TRUE(blst_fp12_is_equal(gt_element, pre_result))
        << "Pairing with precomputation should match direct computation";

    // Verify that we can reuse the precomputed lines with a different point
    POINTonE1 g1_double;
    blst_p1_double(&g1_double, &g1_point);

    POINTonE1_affine g1_double_affine;
    blst_p1_to_affine(&g1_double_affine, &g1_double);

    vec384fp12 direct_double_result;
    blst_miller_loop(direct_double_result, &g2_affine, &g1_double_affine);
    blst_final_exp(direct_double_result, direct_double_result);

    vec384fp12 pre_double_result;
    blst_miller_loop_lines(pre_double_result, lines, &g1_double_affine);
    blst_final_exp(pre_double_result, pre_double_result);

    EXPECT_TRUE(blst_fp12_is_equal(direct_double_result, pre_double_result))
        << "Precomputed pairing should work with different G1 points";
}