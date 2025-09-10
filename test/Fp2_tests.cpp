#include <gtest/gtest.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>

extern "C" {
#include "src/fields.h"
}

// Helper: Print Fp2 element for debugging
void print_fp2(const vec384x a, const char* name) {
    printf("%s:\n", name);
    printf("  Re: ");
    for (int i = NLIMBS(384) - 1; i >= 0; i--) {
        printf("%016llx ", (unsigned long long)a[0][i]);
    }
    printf("\n  Im: ");
    for (int i = NLIMBS(384) - 1; i >= 0; i--) {
        printf("%016llx ", (unsigned long long)a[1][i]);
    }
    printf("\n");
}




TEST(Fp2Test, ZeroAndOne) {
    vec384x zero, one;

    // Initialize zero and one
    vec_zero(zero, sizeof(zero));
    vec_zero(one, sizeof(one));
    vec_copy(one[0], BLS12_381_Rx.p, sizeof(vec384));  // Re(1) = 1

    // Check zero properties
    vec384x result;
    add_fp2(result, zero, zero);
    EXPECT_TRUE(vec_is_equal(result, zero, sizeof(vec384x))) << "0 + 0 = 0";

    mul_fp2(result, zero, one);
    EXPECT_TRUE(vec_is_equal(result, zero, sizeof(vec384x))) << "0 * 1 = 0";

    // Check one properties
    mul_fp2(result, one, one);
    EXPECT_TRUE(vec_is_equal(result, one, sizeof(vec384x))) << "1 * 1 = 1";

    vec384x temp;
    vec_copy(temp, one, sizeof(vec384x));
    add_fp2(result, temp, temp);

    // 1 + 1 should be 2 (in Montgomery representation)
    vec384x expected_two;
    vec_zero(expected_two, sizeof(expected_two));
    lshift_fp(expected_two[0], one[0], 1);  // Re(2) = 2

    EXPECT_TRUE(vec_is_equal(result, expected_two, sizeof(vec384x))) << "1 + 1 = 2";
}

TEST(Fp2Test, Addition) {
    vec384x a, b, c, expected;

    // Initialize with some values
    // a = 2 + 3i (in Montgomery representation)
    vec_zero(a, sizeof(a));
    lshift_fp(a[0], BLS12_381_Rx.p, 1);     // Re(a) = 2
    mul_by_3_fp(a[1], BLS12_381_Rx.p);      // Im(a) = 3

    // b = 4 + 5i (in Montgomery representation)
    vec_zero(b, sizeof(b));
    lshift_fp(b[0], BLS12_381_Rx.p, 2);     // Re(b) = 4
    add_fp(b[1], BLS12_381_Rx.p, b[0]);     // Im(b) = 1 + 4 = 5

    // c = a + b should be 6 + 8i
    add_fp2(c, a, b);

    vec_zero(expected, sizeof(expected));
    // Calculate Re(expected) = 6
    mul_by_3_fp(expected[0], BLS12_381_Rx.p);       // Re(expected) = 3
    lshift_fp(expected[0], expected[0], 1);         // Re(expected) = 6
    // Calculate Im(expected) = 8
    lshift_fp(expected[1], BLS12_381_Rx.p, 3);      // Im(expected) = 8

    EXPECT_TRUE(vec_is_equal(c, expected, sizeof(vec384x)))
        << "(2+3i) + (4+5i) = (6+8i)";

    // Test commutativity
    add_fp2(c, b, a);
    EXPECT_TRUE(vec_is_equal(c, expected, sizeof(vec384x)))
        << "(4+5i) + (2+3i) = (6+8i)";
}

TEST(Fp2Test, Subtraction) {
    vec384x a, b, c, expected;

    // Initialize with some values
    // a = 7 + 4i (in Montgomery representation)
    vec_zero(a, sizeof(a));
    mul_by_3_fp(a[0], BLS12_381_Rx.p);      // Re(a) = 3
    add_fp(a[0], a[0], BLS12_381_Rx.p);     // Re(a) = 3+1 = 4
    add_fp(a[0], a[0], BLS12_381_Rx.p);     // Re(a) = 4+1 = 5
    add_fp(a[0], a[0], BLS12_381_Rx.p);     // Re(a) = 5+1 = 6
    add_fp(a[0], a[0], BLS12_381_Rx.p);     // Re(a) = 6+1 = 7
    lshift_fp(a[1], BLS12_381_Rx.p, 2);     // Im(a) = 4

    // b = 2 + 3i (in Montgomery representation)
    vec_zero(b, sizeof(b));
    lshift_fp(b[0], BLS12_381_Rx.p, 1);     // Re(b) = 2
    mul_by_3_fp(b[1], BLS12_381_Rx.p);      // Im(b) = 3

    // c = a - b should be 5 + 1i
    sub_fp2(c, a, b);

    vec_zero(expected, sizeof(expected));
    add_fp(expected[0], BLS12_381_Rx.p, BLS12_381_Rx.p);  // Re(expected) = 2
    add_fp(expected[0], expected[0], BLS12_381_Rx.p);     // Re(expected) = 2+1 = 3
    add_fp(expected[0], expected[0], BLS12_381_Rx.p);     // Re(expected) = 3+1 = 4
    add_fp(expected[0], expected[0], BLS12_381_Rx.p);     // Re(expected) = 4+1 = 5
    vec_copy(expected[1], BLS12_381_Rx.p, sizeof(vec384));  // Im(expected) = 1

    EXPECT_TRUE(vec_is_equal(c, expected, sizeof(vec384x)))
        << "(7+4i) - (2+3i) = (5+1i)";

    // c = b - a should be -5 - 1i
    sub_fp2(c, b, a);

    // Negate expected to get -5 - 1i
    cneg_fp2(expected, expected, 1);

    EXPECT_TRUE(vec_is_equal(c, expected, sizeof(vec384x)))
        << "(2+3i) - (7+4i) = (-5-1i)";
}

TEST(Fp2Test, Multiplication) {
    vec384x a, b, c, expected;

    // a = 2 + 3i
    vec_zero(a, sizeof(a));
    lshift_fp(a[0], BLS12_381_Rx.p, 1);     // Re(a) = 2
    mul_by_3_fp(a[1], BLS12_381_Rx.p);      // Im(a) = 3

    // b = 4 + 5i
    vec_zero(b, sizeof(b));
    lshift_fp(b[0], BLS12_381_Rx.p, 2);     // Re(b) = 4
    add_fp(b[1], BLS12_381_Rx.p, b[0]);     // Im(b) = 1 + 4 = 5

    // c = a * b should be (2*4 - 3*5) + (2*5 + 3*4)i = (8 - 15) + (10 + 12)i = -7 + 22i
    mul_fp2(c, a, b);

    vec_zero(expected, sizeof(expected));
    // Re(expected) = -7
    mul_by_3_fp(expected[0], BLS12_381_Rx.p);      // 3
    lshift_fp(expected[0], expected[0], 1);        // 6
    add_fp(expected[0], expected[0], BLS12_381_Rx.p); // 7
    cneg_fp(expected[0], expected[0], 1);           // -7

    // Im(expected) = 22
    vec384 temp;
    lshift_fp(temp, BLS12_381_Rx.p, 1);             // 2
    mul_by_3_fp(expected[1], BLS12_381_Rx.p);       // 3
    lshift_fp(expected[1], expected[1], 3);         // 24
    sub_fp(expected[1], expected[1], temp);         // 22

    EXPECT_TRUE(vec_is_equal(c, expected, sizeof(vec384x)))
        << "(2+3i) * (4+5i) = (-7+22i)";

    // Test commutativity
    mul_fp2(c, b, a);
    EXPECT_TRUE(vec_is_equal(c, expected, sizeof(vec384x)))
        << "(4+5i) * (2+3i) = (-7+22i)";
}

TEST(Fp2Test, SquareAndMultiply) {
    vec384x a, a_squared1, a_squared2;

    // a = 2 + 3i
    vec_zero(a, sizeof(a));
    lshift_fp(a[0], BLS12_381_Rx.p, 1);     // Re(a) = 2
    mul_by_3_fp(a[1], BLS12_381_Rx.p);      // Im(a) = 3

    // Compute a² using multiplication
    mul_fp2(a_squared1, a, a);

    // Compute a² using squaring
    sqr_fp2(a_squared2, a);

    // Results should be equal: (2+3i)² = 4-9 + 12i = -5 + 12i
    EXPECT_TRUE(vec_is_equal(a_squared1, a_squared2, sizeof(vec384x)))
        << "(2+3i)² computed both ways should be equal";
}

TEST(Fp2Test, Inversion) {
    vec384x a, a_inv, result;

    // a = 2 + 3i
    vec_zero(a, sizeof(a));
    lshift_fp(a[0], BLS12_381_Rx.p, 1);     // Re(a) = 2
    mul_by_3_fp(a[1], BLS12_381_Rx.p);      // Im(a) = 3

    // Compute 1/a
    reciprocal_fp2(a_inv, a);

    // Check that a * (1/a) = 1
    mul_fp2(result, a, a_inv);

    // Should equal 1 (real part 1, imaginary part 0)
    vec384x one;
    vec_zero(one, sizeof(one));
    vec_copy(one[0], BLS12_381_Rx.p, sizeof(vec384));

    EXPECT_TRUE(vec_is_equal(result, one, sizeof(vec384x)))
        << "a * (1/a) should equal 1";
}

TEST(Fp2Test, Frobenius) {
    vec384x a, frob_a, expected;

    // a = 2 + 3i
    vec_zero(a, sizeof(a));
    lshift_fp(a[0], BLS12_381_Rx.p, 1);     // Re(a) = 2
    mul_by_3_fp(a[1], BLS12_381_Rx.p);      // Im(a) = 3

    // Frobenius map in Fp2 is complex conjugation: (a+bi) -> (a-bi)
    // expected = 2 - 3i
    vec_copy(expected, a, sizeof(vec384x));
    cneg_fp(expected[1], expected[1], 1);  // Negate imaginary part

    // Apply Frobenius
    conjugate_fp2(frob_a, a);

    EXPECT_TRUE(vec_is_equal(frob_a, expected, sizeof(vec384x)))
        << "Frobenius(2+3i) should be 2-3i";
}

TEST(Fp2Test, SquareRoot) {
    vec384x a, root, root_squared;

    // Let's try with a perfect square: a = (2+3i)²
    vec384x temp;
    vec_zero(temp, sizeof(temp));
    lshift_fp(temp[0], BLS12_381_Rx.p, 1);    // Re(temp) = 2
    mul_by_3_fp(temp[1], BLS12_381_Rx.p);     // Im(temp) = 3

    sqr_fp2(a, temp);  // a = (2+3i)²

    // Find square root
    EXPECT_TRUE(sqrt_fp2(root, a)) << "Square root should succeed for a perfect square";

    // Verify by squaring the root
    sqr_fp2(root_squared, root);

    // Should get back our original value
    EXPECT_TRUE(vec_is_equal(root_squared, a, sizeof(vec384x)))
        << "sqrt(a)² should equal a";

    // Verify that either root or -root equals our original value
    bool matches_orig = vec_is_equal(root, temp, sizeof(vec384x));
    cneg_fp2(root, root, 1);
    bool matches_neg = vec_is_equal(root, temp, sizeof(vec384x));

    EXPECT_TRUE(matches_orig || matches_neg)
        << "Either sqrt(a) or -sqrt(a) should equal the original number";
}

TEST(Fp2Test, Serialization) {
    vec384x a, b;
    unsigned char serialized[96]; // 2 * 48 bytes

    // a = 2 + 3i
    vec_zero(a, sizeof(a));
    lshift_fp(a[0], BLS12_381_Rx.p, 1);     // Re(a) = 2
    mul_by_3_fp(a[1], BLS12_381_Rx.p);      // Im(a) = 3

    // Serialize to bytes: convert from Montgomery to standard form
    vec384 temp_re, temp_im;
    from_fp(temp_re, a[0]);  // Convert real part from Montgomery
    from_fp(temp_im, a[1]);  // Convert imaginary part from Montgomery

    // Copy to serialized buffer
    memcpy(serialized, temp_re, 48);
    memcpy(serialized + 48, temp_im, 48);

    // Deserialize back: copy from buffer then convert to Montgomery form
    vec384 raw_re, raw_im;
    memcpy(raw_re, serialized, 48);
    memcpy(raw_im, serialized + 48, 48);

    // Convert to Montgomery form by multiplying by R^2
    mul_fp(b[0], raw_re, BLS12_381_RR);
    mul_fp(b[1], raw_im, BLS12_381_RR);

    // Should get back our original value
    EXPECT_TRUE(vec_is_equal(a, b, sizeof(vec384x)))
        << "Deserialized value should equal original";
}