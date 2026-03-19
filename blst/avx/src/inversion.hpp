

// A 425 step addition chain for inversion via exponentiation in the BLS381 field
// Computes z^(-1) = z^(q - 2) = z^(3 * q - 4), where q is the BLS381 modulus.
// 3 * q - 4 = x**6 + 2*x**5 - 2*x**3 - x - 3 for the BLS12-381 X value
template<class Ring>
class BLS381AddChainInversion {
    using StandardElement = typename Ring::StandardElement;

    const Ring &ring;

    public:
    BLS381AddChainInversion(const Ring &ring) : ring(ring) {}

    StandardElement mul(const StandardElement &a, const StandardElement &b) const {
        return ring.modmul(a, b);
    }

    StandardElement cube(const StandardElement &a) const {
        return mul(mul(a, a), a);
    }

    StandardElement shift(const StandardElement &a, int bits) const {
        StandardElement result = a;
        for (int i = 0; i < bits; i++) {
            result = mul(result, result);
        }
        return result;
    }

    StandardElement exp_x(const StandardElement &base) const {
        // From addchain search
        auto base10 = mul(base, base);
        auto base11 = mul(base10, base);
        auto base1100 = shift(base11, 2);
        auto base1101 = mul(base1100, base);
        auto base1101000 = shift(base1101, 3);
        auto base1101001 = mul(base1101000, base);
        return shift(mul(shift(mul(shift(base1101001, 9), base), 32), base), 16);
    }

    StandardElement exp_x_m3(const StandardElement &base) const {
        // From addchain search
        auto base2 = mul(base, base);
        auto base3 = mul(base2, base);
        auto base12 = shift(base3, 2);
        auto base13 = mul(base12, base);
        auto base15 = mul(base13, base2);
        auto base104 = shift(base13, 3);
        auto base105 = mul(base104, base);
        auto i610 = shift(base105, 9);
        auto i611 = mul(i610, base);
        auto i612 = shift(i611, 36);
        auto i614 = mul(i612, base15);
        auto i61 = shift(i614, 4);
        auto r1 = mul(base15, i61);
        auto r2 = shift(r1, 4);
        auto r3 = mul(r2, base15);
        auto r4 = shift(r3, 4);
        return mul(r4, base13);

    }

    StandardElement invert(const StandardElement &base, const Ring &ring) const {
        (void)ring;  // Suppress unused parameter warning - using member ring instead
        // Evaluate the polynomial x**6 + 2*x**5 - 2*x**3 - x - 3 using the exp_x and exp_x_m3 functions
        auto a = base;
        auto ax_m3 = exp_x_m3(a);
        auto ax_m2 = mul(ax_m3, a);
        auto ax2_m2x = exp_x(ax_m2);
        auto ax = mul(mul(ax_m2, a), a);
        auto ax2_mx = mul(ax2_m2x, ax);
        auto ax3_mx2 = exp_x(ax2_mx);
        auto ax4_mx3 = exp_x(ax3_mx2);
        auto ax5_mx4 = exp_x(ax4_mx3);
        auto ax6_mx5 = exp_x(ax5_mx4);
        auto ax2_mx_m3 = mul(ax2_m2x, ax_m3);
        auto ax3_mx_m3 = mul(ax2_mx_m3, ax3_mx2);
        auto a3x4_m3x2 = cube(ax4_mx3);
        auto a3x4_m2x3_mx_m3 = mul(a3x4_m3x2, ax3_mx_m3);
        auto a3x5_m3x4 = cube(ax5_mx4);
        auto a3x5_m2x3_mx_m2 = mul(a3x5_m3x4, a3x4_m2x3_mx_m3);
        auto ax6_2x5_m2x3_mx_m2 = mul(ax6_mx5, a3x5_m2x3_mx_m2);
        return ax6_2x5_m2x3_mx_m2;
    }
};