#include "BTE/transparent_bte/core.h"

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cstddef>
#include <memory>
#include <stdexcept>

#include "BTE/building_blocks/blst_cpp.h"
#include "BTE/building_blocks/lagrange.h"
#include "blst/fft/context.h"
#include "blst/fft/ctFft.h"
#include "blst/fft/domain.h"
#include "blst/scalar/curve/GT.h"
#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
#include "bindings/blst_avx.h"
#include "blst/avx/curve/GT.h"
#endif

namespace bte::transparent_bte
{
namespace
{
    blst::Fr fr_zero() noexcept
    {
        return blst::Fr{};
    }

    bool fr_is_zero(const blst::Fr& value) noexcept
    {
        return blst::Fr::equal(value, fr_zero()) != 0;
    }

    blst::Fr fr_add(const blst::Fr& lhs, const blst::Fr& rhs) noexcept
    {
        blst::Fr out{};
        blst::Fr::add(out, lhs, rhs);
        return out;
    }

    blst::Fr fr_sub(const blst::Fr& lhs, const blst::Fr& rhs) noexcept
    {
        blst::Fr out{};
        blst::Fr::sub(out, lhs, rhs);
        return out;
    }

    blst::Fr fr_mul(const blst::Fr& lhs, const blst::Fr& rhs) noexcept
    {
        blst::Fr out{};
        blst::Fr::mul(out, lhs, rhs);
        return out;
    }

    blst::Fr fr_neg(const blst::Fr& value) noexcept
    {
        blst::Fr out{};
        blst::Fr::neg(out, value);
        return out;
    }

    blst::Fr fr_inverse(const blst::Fr& value) noexcept
    {
        blst::Fr out{};
        blst::Fr::reciprocal(out, value);
        return out;
    }

    std::array<std::uint8_t, 32> fr_bytes(const blst::Fr& value) noexcept
    {
        std::array<std::uint8_t, 32> out{};
        blst::Fr::to_scaler(out.data(), value);
        return out;
    }

    blst::G1Point add_g1(const blst::G1Point& lhs, const blst::G1Point& rhs)
    {
        blst::G1Point out{};
        blst::G1Point::dadd(out, lhs, rhs);
        return out;
    }

    blst::G2Point add_g2(const blst::G2Point& lhs, const blst::G2Point& rhs)
    {
        blst::G2Point out{};
        blst::G2Point::dadd(out, lhs, rhs);
        return out;
    }

    blst::GT add_gt(const blst::GT& lhs, const blst::GT& rhs)
    {
        blst::GT out{};
        blst::GT::add(out, lhs, rhs);
        return out;
    }

    blst::GT sub_gt(const blst::GT& lhs, const blst::GT& rhs)
    {
        blst::GT negated{};
        blst::GT::neg(negated, rhs);
        return add_gt(lhs, negated);
    }

    template<typename Group>
    Group group_identity();

    template<>
    blst::G1Point group_identity<blst::G1Point>()
    {
        return blst::G1Point::infinity();
    }

    template<>
    blst::G2Point group_identity<blst::G2Point>()
    {
        return blst::G2Point::infinity();
    }

    template<typename Group>
    Group group_add(const Group& lhs, const Group& rhs)
    {
        Group out{};
        Group::dadd(out, lhs, rhs);
        return out;
    }

    template<typename Group>
    Group group_sub(const Group& lhs, const Group& rhs)
    {
        Group negated{};
        Group::neg(negated, rhs);
        return group_add(lhs, negated);
    }

    template<typename Group>
    Group scale_public(const Group& point, const blst::Fr& scalar)
    {
        if (fr_is_zero(scalar)) {
            return group_identity<Group>();
        }
        const auto bytes = fr_bytes(scalar);
        Group out{};
        Group::mult_vartime(out, point, bytes.data(), 255);
        return out;
    }

    template<typename Group>
    Group scale_secret(const Group& point, const blst::Fr& scalar)
    {
        if (fr_is_zero(scalar)) {
            return group_identity<Group>();
        }
        auto bytes = fr_bytes(scalar);
        Group out{};
        Group::mult(out, point, bytes.data(), 255);
        std::fill(bytes.begin(), bytes.end(), 0);
        return out;
    }

    blst::GT scale_secret_gt(const blst::GT& value, const blst::Fr& scalar)
    {
        if (fr_is_zero(scalar)) {
            return blst::GT::identity();
        }
        auto bytes = fr_bytes(scalar);
        blst::GT out{};
        blst::GT::mult(out, value, bytes.data(), 255);
        std::fill(bytes.begin(), bytes.end(), 0);
        return out;
    }

    blst::GT pairing(const blst::G1Point& p, const blst::G2Point& q)
    {
        blst::G1AffinePoint p_affine{};
        blst::G2AffinePoint q_affine{};
        blst::G1Point::to_affine(p_affine, p);
        blst::G2Point::to_affine(q_affine, q);
        blst::MillerResult miller{};
        blst::GT::miller_loop(miller, &q_affine, &p_affine);
        blst::GT out{};
        blst::GT::final_exp(out, miller);
        return out;
    }

    blst::GT pairing_sum_two(
        const blst::G1AffinePoint& p0,
        const blst::G2AffinePoint& q0,
        const blst::G1AffinePoint& p1,
        const blst::G2AffinePoint& q1)
    {
        const blst::G1AffinePoint* ps[2] = { &p0, &p1 };
        const blst::G2AffinePoint* qs[2] = { &q0, &q1 };
        blst::MillerResult miller{};
        blst::GT::miller_loop_n(miller, qs, ps, 2);
        blst::GT out{};
        blst::GT::final_exp(out, miller);
        return out;
    }

    blst::Fr evaluate_scalar_polynomial(
        const std::span<const blst::Fr> coefficients,
        const blst::Fr& point)
    {
        blst::Fr out = fr_zero();
        for (auto it = coefficients.rbegin(); it != coefficients.rend(); ++it) {
            out = fr_add(fr_mul(out, point), *it);
        }
        return out;
    }

    void require_distinct(const std::span<const blst::Fr> points)
    {
        for (std::size_t i = 0; i < points.size(); ++i) {
            for (std::size_t j = 0; j < i; ++j) {
                if (blst::Fr::equal(points[i], points[j]) != 0) {
                    throw std::invalid_argument("transparent BTE ciphertext hash collision");
                }
            }
        }
    }

    std::vector<blst::Fr> ciphertext_points(
        const std::span<const Ciphertext> ciphertexts)
    {
        if (ciphertexts.empty()) {
            throw std::invalid_argument("transparent BTE batch must be nonempty");
        }
        std::vector<blst::Fr> points;
        points.reserve(ciphertexts.size());
        for (const Ciphertext& ciphertext : ciphertexts) {
            points.push_back(CiphertextPoint(ciphertext));
        }
        require_distinct(points);
        return points;
    }

    std::vector<blst::Fr> batch_invert_nonzero(
        const std::span<const blst::Fr> values)
    {
        if (values.empty()) {
            return {};
        }
        std::vector<blst::Fr> prefixes(values.size() + 1, blst::Fr::one());
        for (std::size_t i = 0; i < values.size(); ++i) {
            if (fr_is_zero(values[i])) {
                throw std::invalid_argument("transparent BTE attempted to invert zero");
            }
            prefixes[i + 1] = fr_mul(prefixes[i], values[i]);
        }

        blst::Fr suffix = fr_inverse(prefixes.back());
        std::vector<blst::Fr> inverses(values.size(), fr_zero());
        for (std::size_t i = values.size(); i-- > 0;) {
            inverses[i] = fr_mul(suffix, prefixes[i]);
            suffix = fr_mul(suffix, values[i]);
        }
        return inverses;
    }

    std::size_t next_power_of_two(const std::size_t value)
    {
        if (value <= 1) {
            return 1;
        }
        return std::size_t{1} << std::bit_width(value - 1);
    }

    class PolynomialArithmetic
    {
      public:
        explicit PolynomialArithmetic(const std::size_t maximum_product_size)
        {
            std::size_t size = 2;
            const std::size_t maximum_fft = next_power_of_two(maximum_product_size);
            while (size <= maximum_fft) {
                domains_.push_back(Domain{
                    size,
                    blst::fft_detail::make_roots<blst::Fr, blst::FrCoeffCtx>(size, false),
                    fr_inverse(blst::Fr::from_uint64(size)) });
                size <<= 1;
            }
        }

        std::vector<blst::Fr> scalar_multiply(
            const std::span<const blst::Fr> lhs,
            const std::span<const blst::Fr> rhs) const
        {
            if (lhs.empty() || rhs.empty()) {
                return {};
            }
            constexpr std::size_t naive_threshold = 32;
            const std::size_t result_size = lhs.size() + rhs.size() - 1;
            if (std::min(lhs.size(), rhs.size()) <= naive_threshold) {
                std::vector<blst::Fr> out(result_size, fr_zero());
                for (std::size_t i = 0; i < lhs.size(); ++i) {
                    for (std::size_t j = 0; j < rhs.size(); ++j) {
                        out[i + j] = fr_add(out[i + j], fr_mul(lhs[i], rhs[j]));
                    }
                }
                return out;
            }

            const std::size_t fft_size = next_power_of_two(result_size);
            const Domain& domain = get_domain(fft_size);
            std::vector<blst::Fr> a(fft_size, fr_zero());
            std::vector<blst::Fr> b(fft_size, fr_zero());
            std::copy(lhs.begin(), lhs.end(), a.begin());
            std::copy(rhs.begin(), rhs.end(), b.begin());
            fft<blst::Fr, blst::Fr, blst::FrCoeffCtx>(a, domain.roots, false);
            fft<blst::Fr, blst::Fr, blst::FrCoeffCtx>(b, domain.roots, false);
            for (std::size_t i = 0; i < fft_size; ++i) {
                a[i] = fr_mul(a[i], b[i]);
            }
            fft_with_forward_roots<blst::Fr, blst::Fr, blst::FrCoeffCtx>(
                a, domain.roots, true, domain.inverse_size);
            a.resize(result_size);
            return a;
        }

        template<typename Group>
        std::vector<Group> group_scalar_multiply(
            const std::span<const Group> lhs,
            const std::span<const blst::Fr> rhs) const
        {
            if (lhs.empty() || rhs.empty()) {
                return {};
            }
            constexpr std::size_t naive_threshold = 16;
            const std::size_t result_size = lhs.size() + rhs.size() - 1;
            if (std::min(lhs.size(), rhs.size()) <= naive_threshold) {
                std::vector<Group> out(result_size, group_identity<Group>());
                for (std::size_t i = 0; i < lhs.size(); ++i) {
                    for (std::size_t j = 0; j < rhs.size(); ++j) {
                        out[i + j] = group_add(out[i + j], scale_public(lhs[i], rhs[j]));
                    }
                }
                return out;
            }

            const std::size_t fft_size = next_power_of_two(result_size);
            const Domain& domain = get_domain(fft_size);
            std::vector<Group> a(fft_size, group_identity<Group>());
            std::vector<blst::Fr> b(fft_size, fr_zero());
            std::copy(lhs.begin(), lhs.end(), a.begin());
            std::copy(rhs.begin(), rhs.end(), b.begin());
            fft<Group, blst::Fr, blst::CoeffCtx>(a, domain.roots, false);
            fft<blst::Fr, blst::Fr, blst::FrCoeffCtx>(b, domain.roots, false);
            for (std::size_t i = 0; i < fft_size; ++i) {
                a[i] = scale_public(a[i], b[i]);
            }
            fft_with_forward_roots<Group, blst::Fr, blst::CoeffCtx>(
                a, domain.roots, true, domain.inverse_size);
            a.resize(result_size);
            return a;
        }

        std::vector<blst::Fr> scalar_remainder(
            const std::span<const blst::Fr> dividend,
            const std::span<const blst::Fr> divisor) const
        {
            if (divisor.empty()) {
                throw std::invalid_argument("polynomial divisor must be nonempty");
            }
            if (dividend.size() < divisor.size()) {
                return { dividend.begin(), dividend.end() };
            }
            constexpr std::size_t naive_threshold = 32;
            if (divisor.size() <= naive_threshold) {
                std::vector<blst::Fr> remainder(dividend.begin(), dividend.end());
                const blst::Fr inverse_lead = fr_inverse(divisor.back());
                const std::size_t quotient_size = dividend.size() - divisor.size() + 1;
                for (std::size_t i = quotient_size; i-- > 0;) {
                    const blst::Fr factor = fr_mul(
                        remainder[i + divisor.size() - 1], inverse_lead);
                    for (std::size_t j = 0; j < divisor.size(); ++j) {
                        remainder[i + j] = fr_sub(
                            remainder[i + j], fr_mul(factor, divisor[j]));
                    }
                }
                remainder.resize(divisor.size() - 1);
                return remainder;
            }

            const std::size_t quotient_size = dividend.size() - divisor.size() + 1;
            const std::vector<blst::Fr> reversed_dividend = reverse(dividend);
            const std::vector<blst::Fr> reversed_divisor = reverse(divisor);
            std::vector<blst::Fr> dividend_prefix(
                reversed_dividend.begin(), reversed_dividend.begin() + quotient_size);
            const std::vector<blst::Fr> divisor_inverse =
                inverse_series(reversed_divisor, quotient_size);
            std::vector<blst::Fr> quotient_reversed =
                scalar_multiply(dividend_prefix, divisor_inverse);
            quotient_reversed.resize(quotient_size, fr_zero());
            const std::vector<blst::Fr> quotient = reverse(
                std::span<const blst::Fr>(quotient_reversed));
            const std::vector<blst::Fr> product = scalar_multiply(quotient, divisor);
            std::vector<blst::Fr> remainder(dividend.size(), fr_zero());
            for (std::size_t i = 0; i < dividend.size(); ++i) {
                remainder[i] = fr_sub(dividend[i], product[i]);
            }
            remainder.resize(divisor.size() - 1);
            return remainder;
        }

        template<typename Group>
        std::vector<Group> group_remainder(
            const std::span<const Group> dividend,
            const std::span<const blst::Fr> divisor) const
        {
            if (divisor.empty()) {
                throw std::invalid_argument("group polynomial divisor must be nonempty");
            }
            if (dividend.size() < divisor.size()) {
                return { dividend.begin(), dividend.end() };
            }
            constexpr std::size_t naive_threshold = 16;
            if (divisor.size() <= naive_threshold) {
                std::vector<Group> remainder(dividend.begin(), dividend.end());
                const blst::Fr inverse_lead = fr_inverse(divisor.back());
                const std::size_t quotient_size = dividend.size() - divisor.size() + 1;
                for (std::size_t i = quotient_size; i-- > 0;) {
                    const Group factor = scale_public(
                        remainder[i + divisor.size() - 1], inverse_lead);
                    for (std::size_t j = 0; j < divisor.size(); ++j) {
                        remainder[i + j] = group_sub(
                            remainder[i + j], scale_public(factor, divisor[j]));
                    }
                }
                remainder.resize(divisor.size() - 1);
                return remainder;
            }

            const std::size_t quotient_size = dividend.size() - divisor.size() + 1;
            const std::vector<Group> reversed_dividend = reverse(dividend);
            const std::vector<blst::Fr> reversed_divisor = reverse(divisor);
            std::vector<Group> dividend_prefix(
                reversed_dividend.begin(), reversed_dividend.begin() + quotient_size);
            const std::vector<blst::Fr> divisor_inverse =
                inverse_series(reversed_divisor, quotient_size);
            std::vector<Group> quotient_reversed =
                group_scalar_multiply<Group>(dividend_prefix, divisor_inverse);
            quotient_reversed.resize(quotient_size, group_identity<Group>());
            const std::vector<Group> quotient = reverse(
                std::span<const Group>(quotient_reversed));
            const std::vector<Group> product = group_scalar_multiply<Group>(quotient, divisor);
            std::vector<Group> remainder(dividend.size(), group_identity<Group>());
            for (std::size_t i = 0; i < dividend.size(); ++i) {
                remainder[i] = group_sub(dividend[i], product[i]);
            }
            remainder.resize(divisor.size() - 1);
            return remainder;
        }

      private:
        struct Domain
        {
            std::size_t size = 0;
            std::vector<blst::Fr> roots;
            blst::Fr inverse_size{};
        };

        std::vector<Domain> domains_;

        const Domain& get_domain(const std::size_t size) const
        {
            const auto it = std::find_if(
                domains_.begin(), domains_.end(),
                [size](const Domain& domain) { return domain.size == size; });
            if (it == domains_.end()) {
                throw std::length_error("unsupported polynomial FFT size");
            }
            return *it;
        }

        template<typename T>
        static std::vector<T> reverse(const std::span<const T> values)
        {
            return std::vector<T>(values.rbegin(), values.rend());
        }

        std::vector<blst::Fr> inverse_series(
            const std::span<const blst::Fr> polynomial,
            const std::size_t size) const
        {
            if (size == 0) {
                return {};
            }
            if (polynomial.empty() || fr_is_zero(polynomial.front())) {
                throw std::invalid_argument("polynomial series has zero constant");
            }
            std::vector<blst::Fr> inverse(1, fr_inverse(polynomial.front()));
            while (inverse.size() < size) {
                const std::size_t next = std::min(inverse.size() * 2, size);
                std::vector<blst::Fr> prefix(
                    polynomial.begin(),
                    polynomial.begin() + std::min(polynomial.size(), next));
                std::vector<blst::Fr> product = scalar_multiply(prefix, inverse);
                product.resize(next, fr_zero());
                product[0] = fr_sub(blst::Fr::from_uint64(2), product[0]);
                for (std::size_t i = 1; i < next; ++i) {
                    product[i] = fr_neg(product[i]);
                }
                inverse = scalar_multiply(inverse, product);
                inverse.resize(next, fr_zero());
            }
            return inverse;
        }
    };

    struct ProductTreeNode
    {
        std::size_t begin = 0;
        std::size_t end = 0;
        std::vector<blst::Fr> polynomial;
        std::unique_ptr<ProductTreeNode> left;
        std::unique_ptr<ProductTreeNode> right;
    };

    std::unique_ptr<ProductTreeNode> build_product_tree(
        const PolynomialArithmetic& arithmetic,
        const std::span<const blst::Fr> points,
        const std::size_t begin,
        const std::size_t end)
    {
        auto node = std::make_unique<ProductTreeNode>();
        node->begin = begin;
        node->end = end;
        if (end - begin == 1) {
            node->polynomial = { fr_neg(points[begin]), blst::Fr::one() };
            return node;
        }
        const std::size_t middle = begin + (end - begin) / 2;
        node->left = build_product_tree(arithmetic, points, begin, middle);
        node->right = build_product_tree(arithmetic, points, middle, end);
        node->polynomial = arithmetic.scalar_multiply(
            node->left->polynomial, node->right->polynomial);
        return node;
    }

    template<typename Group>
    std::vector<Group> build_numerator_polynomial(
        const PolynomialArithmetic& arithmetic,
        const ProductTreeNode& node,
        const std::span<const Group> values)
    {
        if (node.end - node.begin == 1) {
            return { values[node.begin] };
        }
        std::vector<Group> left =
            build_numerator_polynomial(arithmetic, *node.left, values);
        std::vector<Group> right =
            build_numerator_polynomial(arithmetic, *node.right, values);
        std::vector<Group> left_product = arithmetic.group_scalar_multiply<Group>(
            left, node.right->polynomial);
        std::vector<Group> right_product = arithmetic.group_scalar_multiply<Group>(
            right, node.left->polynomial);
        for (std::size_t i = 0; i < left_product.size(); ++i) {
            left_product[i] = group_add(left_product[i], right_product[i]);
        }
        return left_product;
    }

    std::vector<blst::Fr> derivative(
        const std::span<const blst::Fr> polynomial)
    {
        if (polynomial.size() <= 1) {
            return {};
        }
        std::vector<blst::Fr> out(polynomial.size() - 1, fr_zero());
        for (std::size_t i = 1; i < polynomial.size(); ++i) {
            out[i - 1] = fr_mul(polynomial[i], blst::Fr::from_uint64(i));
        }
        return out;
    }

    template<typename Group>
    std::vector<Group> derivative(const std::span<const Group> polynomial)
    {
        if (polynomial.size() <= 1) {
            return {};
        }
        std::vector<Group> out(polynomial.size() - 1, group_identity<Group>());
        for (std::size_t i = 1; i < polynomial.size(); ++i) {
            out[i - 1] = scale_public(polynomial[i], blst::Fr::from_uint64(i));
        }
        return out;
    }

    void evaluate_scalar_recursive(
        const PolynomialArithmetic& arithmetic,
        const ProductTreeNode& node,
        const std::span<const blst::Fr> polynomial,
        std::span<blst::Fr> out)
    {
        if (node.end - node.begin == 1) {
            out[node.begin] = polynomial.empty() ? fr_zero() : polynomial.front();
            return;
        }
        const std::vector<blst::Fr> left =
            arithmetic.scalar_remainder(polynomial, node.left->polynomial);
        const std::vector<blst::Fr> right =
            arithmetic.scalar_remainder(polynomial, node.right->polynomial);
        evaluate_scalar_recursive(arithmetic, *node.left, left, out);
        evaluate_scalar_recursive(arithmetic, *node.right, right, out);
    }

    template<typename Group>
    void evaluate_group_recursive(
        const PolynomialArithmetic& arithmetic,
        const ProductTreeNode& node,
        const std::span<const Group> polynomial,
        std::span<Group> out)
    {
        if (node.end - node.begin == 1) {
            out[node.begin] = polynomial.empty()
                ? group_identity<Group>()
                : polynomial.front();
            return;
        }
        const std::vector<Group> left =
            arithmetic.group_remainder<Group>(polynomial, node.left->polynomial);
        const std::vector<Group> right =
            arithmetic.group_remainder<Group>(polynomial, node.right->polynomial);
        evaluate_group_recursive<Group>(
            arithmetic, *node.left, std::span<const Group>(left), out);
        evaluate_group_recursive<Group>(
            arithmetic, *node.right, std::span<const Group>(right), out);
    }

    template<typename Group>
    std::vector<typename Group::Field_type> affine_prefix_storage(const std::size_t count)
    {
        return std::vector<typename Group::Field_type>(count);
    }

    template<typename Point, typename AffinePoint>
    void to_affine_batch(
        std::vector<AffinePoint>& out,
        const std::vector<Point>& in,
        std::vector<typename Point::Field_type>& prefix_products)
    {
        if constexpr (requires {
            Point::to_affine_batch(
                std::span<AffinePoint>(out),
                std::span<const Point>(in),
                std::span<typename Point::Field_type>(prefix_products));
        }) {
            Point::to_affine_batch(
                std::span<AffinePoint>(out),
                std::span<const Point>(in),
                std::span<typename Point::Field_type>(prefix_products));
        }
        else {
            // The scalar CatBLST point API has no batch conversion primitive.
            // Keep that backend functional without imposing conversions on the
            // native AVX path, which takes the single-inversion batch branch.
            for (std::size_t i = 0; i < in.size(); ++i) {
                Point::to_affine(out[i], in[i]);
            }
        }
    }

    BatchPrecomputation make_precomputation_shell(
        const std::span<const Ciphertext> ciphertexts)
    {
        BatchPrecomputation out{};
        out.cross_terms.resize(ciphertexts.size(), blst::GT::identity());
        out.c2_affine.resize(ciphertexts.size());

        std::vector<blst::G2Point> c2_points;
        c2_points.reserve(ciphertexts.size());
        for (const Ciphertext& ciphertext : ciphertexts) {
            c2_points.push_back(ciphertext.c2);
        }
        auto prefixes = affine_prefix_storage<blst::G2Point>(ciphertexts.size());
        to_affine_batch(out.c2_affine, c2_points, prefixes);
        return out;
    }

    void compute_cross_terms_from_evaluations(
        BatchPrecomputation& out,
        const std::span<const Ciphertext> ciphertexts,
        const std::span<const blst::G2Point> n2_prime_values,
        const std::span<const blst::G1Point> n3_prime_values,
        const std::span<const blst::Fr> inverse_f_prime_values)
    {
        const std::size_t count = ciphertexts.size();
        std::vector<blst::G1Point> p(2 * count, blst::G1Point::infinity());
        std::vector<blst::G2Point> q(2 * count, blst::G2Point::infinity());
        for (std::size_t i = 0; i < count; ++i) {
            const blst::Fr negative_inverse = fr_neg(inverse_f_prime_values[i]);
            p[2 * i] = scale_public(n3_prime_values[i], negative_inverse);
            q[2 * i] = ciphertexts[i].c2;
            p[2 * i + 1] = ciphertexts[i].c3;
            q[2 * i + 1] = scale_public(n2_prime_values[i], inverse_f_prime_values[i]);
        }

        std::vector<blst::G1AffinePoint> p_affine(p.size());
        std::vector<blst::G2AffinePoint> q_affine(q.size());
        auto p_prefixes = affine_prefix_storage<blst::G1Point>(p.size());
        auto q_prefixes = affine_prefix_storage<blst::G2Point>(q.size());
        to_affine_batch(p_affine, p, p_prefixes);
        to_affine_batch(q_affine, q, q_prefixes);

        for (std::size_t i = 0; i < count; ++i) {
            out.cross_terms[i] = pairing_sum_two(
                p_affine[2 * i], q_affine[2 * i],
                p_affine[2 * i + 1], q_affine[2 * i + 1]);
        }
    }
}

SetupResult SetupFromSharingPolynomial(
    const std::span<const blst::Fr> sharing_polynomial,
    const std::size_t num_servers,
    const blst::G1Point& independent_generator)
{
    if (sharing_polynomial.empty() || num_servers < sharing_polynomial.size()) {
        throw std::invalid_argument("transparent BTE sharing shape is invalid");
    }
    if (fr_is_zero(sharing_polynomial.front())) {
        throw std::invalid_argument("transparent BTE secret must be nonzero");
    }

    SetupResult out{};
    out.public_parameters.threshold = sharing_polynomial.size();
    out.public_parameters.independent_generator = independent_generator;
    out.public_parameters.public_key = scale_secret(
        blst::G1Point::generator(), sharing_polynomial.front());
    out.public_parameters.public_key_in_gt = pairing(
        out.public_parameters.public_key, blst::G2Point::generator());
    out.public_parameters.verification_keys.reserve(num_servers);
    out.secret_shares.reserve(num_servers);
    for (std::size_t i = 0; i < num_servers; ++i) {
        const std::uint32_t index = static_cast<std::uint32_t>(i + 1);
        const blst::Fr share = evaluate_scalar_polynomial(
            sharing_polynomial, blst::Fr::from_uint64(index));
        out.secret_shares.push_back(SecretShare{ index, share });
        out.public_parameters.verification_keys.push_back(
            scale_secret(blst::G1Point::generator(), share));
    }
    return out;
}

Ciphertext Encrypt(
    const PublicParameters& public_parameters,
    const blst::Fr& randomness,
    const blst::GT& message)
{
    if (fr_is_zero(randomness)) {
        throw std::invalid_argument("transparent BTE encryption randomness must be nonzero");
    }
    Ciphertext out{};
    out.c1 = scale_secret(blst::G1Point::generator(), randomness);
    out.c2 = scale_secret(blst::G2Point::generator(), randomness);
    const blst::Fr x = CiphertextPoint(out);
    const blst::G1Point x_pk = scale_public(public_parameters.public_key, x);
    const blst::G1Point base = add_g1(public_parameters.independent_generator, x_pk);
    out.c3 = scale_secret(base, randomness);
    const blst::Fr randomness_squared = fr_mul(randomness, randomness);
    const blst::GT mask = scale_secret_gt(
        public_parameters.public_key_in_gt, randomness_squared);
    out.c4 = add_gt(mask, message);
    return out;
}

blst::Fr CiphertextPoint(const Ciphertext& ciphertext)
{
    std::array<std::uint8_t, 48> compressed{};
    blst::G1Point::compress(compressed.data(), ciphertext.c1);
    const auto scalar_bytes = bte::building_blocks::hash_to_scalar_le(
        compressed, "BTE_TRANSPARENT_2026_2205_X_V1");
    blst::Fr out{};
    blst::Fr::from_scalar(out, scalar_bytes.data());
    return out;
}

PartialDecryptionShare PartialDecrypt(
    const SecretShare& secret_share,
    const std::span<const Ciphertext> ciphertexts)
{
    if (ciphertexts.empty()) {
        throw std::invalid_argument("transparent BTE batch must be nonempty");
    }
    blst::G1Point sum = blst::G1Point::infinity();
    for (const Ciphertext& ciphertext : ciphertexts) {
        sum = add_g1(sum, ciphertext.c1);
    }
    return PartialDecryptionShare{
        secret_share.index,
        scale_secret(sum, secret_share.value)
    };
}

bool VerifyPartialDecryptionShare(
    const PublicParameters& public_parameters,
    const PartialDecryptionShare& share,
    const std::span<const Ciphertext> ciphertexts)
{
    if (share.index == 0 || share.index > public_parameters.verification_keys.size() ||
        ciphertexts.empty()) {
        return false;
    }
    blst::G2Point sum_c2 = blst::G2Point::infinity();
    for (const Ciphertext& ciphertext : ciphertexts) {
        sum_c2 = add_g2(sum_c2, ciphertext.c2);
    }
    blst::G1Point negative_vk{};
    blst::G1Point::neg(
        negative_vk,
        public_parameters.verification_keys[share.index - 1]);

    blst::G1AffinePoint p[2]{};
    blst::G2AffinePoint q[2]{};
    blst::G1Point::to_affine(p[0], share.value);
    blst::G2Point::to_affine(q[0], blst::G2Point::generator());
    blst::G1Point::to_affine(p[1], negative_vk);
    blst::G2Point::to_affine(q[1], sum_c2);
    return blst::GT::is_identity(pairing_sum_two(p[0], q[0], p[1], q[1])) != 0;
}

blst::G1Point CombinePartialDecryptionShares(
    const PublicParameters& public_parameters,
    const std::span<const PartialDecryptionShare> shares)
{
    if (shares.size() < public_parameters.threshold) {
        throw std::invalid_argument("transparent BTE has too few shares");
    }
    std::vector<std::uint32_t> indices;
    indices.reserve(shares.size());
    for (const PartialDecryptionShare& share : shares) {
        indices.push_back(share.index);
    }
    const std::vector<blst::Fr> lambdas =
        bte::building_blocks::lagrange_at_zero(indices);
    blst::G1Point out = blst::G1Point::infinity();
    for (std::size_t i = 0; i < shares.size(); ++i) {
        out = add_g1(out, scale_public(shares[i].value, lambdas[i]));
    }
    return out;
}

BatchPrecomputation PrecomputeBatchDecryptionNaive(
    const std::span<const Ciphertext> ciphertexts)
{
    std::vector<blst::Fr> points = ciphertext_points(ciphertexts);
    BatchPrecomputation out = make_precomputation_shell(ciphertexts);
    const std::size_t count = ciphertexts.size();

    std::vector<blst::Fr> differences;
    differences.reserve(count * (count - 1));
    for (std::size_t i = 0; i < count; ++i) {
        for (std::size_t j = 0; j < count; ++j) {
            if (i != j) {
                differences.push_back(fr_sub(points[j], points[i]));
            }
        }
    }
    const std::vector<blst::Fr> inverses = batch_invert_nonzero(differences);
    std::vector<blst::G1Point> w(count, blst::G1Point::infinity());
    std::vector<blst::G2Point> u(count, blst::G2Point::infinity());
    std::size_t offset = 0;
    for (std::size_t i = 0; i < count; ++i) {
        for (std::size_t j = 0; j < count; ++j) {
            if (i == j) {
                continue;
            }
            u[i] = add_g2(u[i], scale_public(ciphertexts[j].c2, inverses[offset]));
            w[i] = add_g1(w[i], scale_public(ciphertexts[j].c3, inverses[offset]));
            ++offset;
        }
    }

    std::vector<blst::G1Point> p(2 * count, blst::G1Point::infinity());
    std::vector<blst::G2Point> q(2 * count, blst::G2Point::infinity());
    for (std::size_t i = 0; i < count; ++i) {
        p[2 * i] = w[i];
        q[2 * i] = ciphertexts[i].c2;
        blst::G1Point::neg(p[2 * i + 1], ciphertexts[i].c3);
        q[2 * i + 1] = u[i];
    }
    std::vector<blst::G1AffinePoint> p_affine(p.size());
    std::vector<blst::G2AffinePoint> q_affine(q.size());
    auto p_prefixes = affine_prefix_storage<blst::G1Point>(p.size());
    auto q_prefixes = affine_prefix_storage<blst::G2Point>(q.size());
    to_affine_batch(p_affine, p, p_prefixes);
    to_affine_batch(q_affine, q, q_prefixes);
    for (std::size_t i = 0; i < count; ++i) {
        out.cross_terms[i] = pairing_sum_two(
            p_affine[2 * i], q_affine[2 * i],
            p_affine[2 * i + 1], q_affine[2 * i + 1]);
    }
    return out;
}

BatchPrecomputation PrecomputeBatchDecryptionMsm(
    const std::span<const Ciphertext> ciphertexts)
{
    return BenchmarkPrecomputeBatchDecryptionMsm(
        ciphertexts, MsmWeightMode::Cauchy, nullptr);
}

BatchPrecomputation BenchmarkPrecomputeBatchDecryptionMsm(
    const std::span<const Ciphertext> ciphertexts,
    const MsmWeightMode weight_mode,
    MsmPrecomputeTimings* const timings)
{
#if BLST_DEFAULT_AVX
    using Clock = std::chrono::steady_clock;
    const auto elapsed_ns = [](const Clock::time_point begin) {
        return static_cast<std::uint64_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                Clock::now() - begin).count());
    };
    if (timings != nullptr) {
        *timings = {};
    }
    const auto total_begin = Clock::now();

    static_assert(sizeof(blst::G1Point) == sizeof(blst_avx_p1));
    static_assert(alignof(blst::G1Point) == alignof(blst_avx_p1));
    static_assert(sizeof(blst::G2Point) == sizeof(blst_avx_p2));
    static_assert(alignof(blst::G2Point) == alignof(blst_avx_p2));
    static_assert(sizeof(blst::Fr) == sizeof(blst_fr));
    static_assert(alignof(blst::Fr) == alignof(blst_fr));

    auto phase_begin = Clock::now();
    std::vector<blst::Fr> points = ciphertext_points(ciphertexts);
    if (timings != nullptr) {
        timings->point_hash_ns = elapsed_ns(phase_begin);
    }

    phase_begin = Clock::now();
    BatchPrecomputation out = make_precomputation_shell(ciphertexts);
    if (timings != nullptr) {
        timings->output_shell_ns = elapsed_ns(phase_begin);
    }
    const std::size_t count = ciphertexts.size();

    phase_begin = Clock::now();
    std::vector<blst::Fr> differences;
    differences.reserve(count * (count - 1));
    for (std::size_t i = 0; i < count; ++i) {
        for (std::size_t j = 0; j < count; ++j) {
            if (i != j) {
                differences.push_back(fr_sub(points[j], points[i]));
            }
        }
    }
    std::vector<blst::Fr> inverses = batch_invert_nonzero(differences);
    if (weight_mode == MsmWeightMode::Unit) {
        std::fill(inverses.begin(), inverses.end(), blst::Fr::one());
    }
    if (timings != nullptr) {
        timings->cauchy_weights_ns = elapsed_ns(phase_begin);
    }

    // Feed the native projective arrays directly to CatBLST.  Keeping a zero
    // scalar in the diagonal slot avoids rebuilding B nearly-identical point
    // arrays; the public-scalar AVX kernels choose shared-doubling schedules
    // for small rows and Pippenger for larger rows.
    phase_begin = Clock::now();
    std::vector<blst::G1Point> c3(count);
    std::vector<blst::G2Point> c2(count);
    for (std::size_t i = 0; i < count; ++i) {
        c3[i] = ciphertexts[i].c3;
        c2[i] = ciphertexts[i].c2;
    }

    struct alignas(BLST_AVX_ALIGNMENT) ScratchBlock
    {
        std::byte bytes[BLST_AVX_ALIGNMENT];
    };
    const std::size_t g1_scratch_bytes =
        blst_avx_p1s_mult_fr_vartime_scratch_sizeof(count);
    const std::size_t g2_scratch_bytes =
        blst_avx_p2s_mult_fr_vartime_scratch_sizeof(count);
    const std::size_t scratch_bytes = std::max(g1_scratch_bytes, g2_scratch_bytes);
    std::vector<ScratchBlock> scratch(
        (scratch_bytes + sizeof(ScratchBlock) - 1) / sizeof(ScratchBlock));
    void* const scratch_data = scratch.empty() ? nullptr : scratch.data();
    const std::size_t scratch_capacity = scratch.size() * sizeof(ScratchBlock);

    std::vector<blst::Fr> row_scalars(count, fr_zero());
    std::vector<blst::G1Point> w(count);
    std::vector<blst::G2Point> u(count);
    if (timings != nullptr) {
        timings->buffer_setup_ns = elapsed_ns(phase_begin);
    }
    std::size_t inverse_offset = 0;
    for (std::size_t i = 0; i < count; ++i) {
        std::fill(row_scalars.begin(), row_scalars.end(), fr_zero());
        for (std::size_t j = 0; j < count; ++j) {
            if (i != j) {
                row_scalars[j] = inverses[inverse_offset++];
            }
        }
        const auto g1_begin = Clock::now();
        const bool g1_ok = blst_avx_p1s_mult_fr_vartime(
            reinterpret_cast<blst_avx_p1*>(&w[i]),
            reinterpret_cast<const blst_avx_p1*>(c3.data()),
            reinterpret_cast<const blst_fr*>(row_scalars.data()),
            count,
            scratch_data,
            scratch_capacity);
        if (timings != nullptr) {
            timings->g1_msm_ns += elapsed_ns(g1_begin);
        }
        const auto g2_begin = Clock::now();
        const bool g2_ok = blst_avx_p2s_mult_fr_vartime(
            reinterpret_cast<blst_avx_p2*>(&u[i]),
            reinterpret_cast<const blst_avx_p2*>(c2.data()),
            reinterpret_cast<const blst_fr*>(row_scalars.data()),
            count,
            scratch_data,
            scratch_capacity);
        if (timings != nullptr) {
            timings->g2_msm_ns += elapsed_ns(g2_begin);
        }
        if (!g1_ok || !g2_ok) {
            throw std::runtime_error("transparent BTE native MSM failed");
        }
    }

    phase_begin = Clock::now();
    std::vector<blst::G1Point> p(2 * count, blst::G1Point::infinity());
    std::vector<blst::G2Point> q(2 * count, blst::G2Point::infinity());
    for (std::size_t i = 0; i < count; ++i) {
        p[2 * i] = w[i];
        q[2 * i] = ciphertexts[i].c2;
        blst::G1Point::neg(p[2 * i + 1], ciphertexts[i].c3);
        q[2 * i + 1] = u[i];
    }
    std::vector<blst::G1AffinePoint> p_affine(p.size());
    std::vector<blst::G2AffinePoint> q_affine(q.size());
    auto p_prefixes = affine_prefix_storage<blst::G1Point>(p.size());
    auto q_prefixes = affine_prefix_storage<blst::G2Point>(q.size());
    to_affine_batch(p_affine, p, p_prefixes);
    to_affine_batch(q_affine, q, q_prefixes);
    if (timings != nullptr) {
        timings->affine_ns = elapsed_ns(phase_begin);
    }

    phase_begin = Clock::now();
    for (std::size_t i = 0; i < count; ++i) {
        out.cross_terms[i] = pairing_sum_two(
            p_affine[2 * i], q_affine[2 * i],
            p_affine[2 * i + 1], q_affine[2 * i + 1]);
    }
    if (timings != nullptr) {
        timings->pairing_ns = elapsed_ns(phase_begin);
        timings->total_ns = elapsed_ns(total_begin);
    }
    return out;
#else
    // The experiment's scalar comparison retains the direct oracle.  The
    // native projective row-MSM API under test is intentionally AVX-specific.
    if (weight_mode != MsmWeightMode::Cauchy) {
        throw std::invalid_argument(
            "transparent BTE unit-weight MSM benchmark requires AVX");
    }
    const auto begin = std::chrono::steady_clock::now();
    BatchPrecomputation out = PrecomputeBatchDecryptionNaive(ciphertexts);
    if (timings != nullptr) {
        *timings = {};
        timings->total_ns = static_cast<std::uint64_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - begin).count());
    }
    return out;
#endif
}

BatchPrecomputation PrecomputeBatchDecryptionQuadraticPairing(
    const std::span<const Ciphertext> ciphertexts)
{
    std::vector<blst::Fr> points = ciphertext_points(ciphertexts);
    BatchPrecomputation out = make_precomputation_shell(ciphertexts);
    const std::size_t count = ciphertexts.size();
    if (count == 1) {
        return out;
    }

    // For i < j, both output slots contain the same contribution
    //
    //   (e(c3_j, c2_i) / e(c3_i, c2_j))^(1 / (x_j - x_i)).
    //
    // Compute its two scaled G1 inputs once, then reuse them in the large
    // prepared-line Miller loop for both incident output rows.  This removes
    // the direct oracle's quadratic G2 scalar multiplications while retaining
    // only B final exponentiations.
    const std::size_t pair_count = count * (count - 1) / 2;
    std::vector<blst::Fr> differences;
    differences.reserve(pair_count);
    for (std::size_t i = 0; i < count; ++i) {
        for (std::size_t j = i + 1; j < count; ++j) {
            differences.push_back(fr_sub(points[j], points[i]));
        }
    }
    const std::vector<blst::Fr> inverses = batch_invert_nonzero(differences);

    std::vector<blst::G1Point> scaled_points(2 * pair_count);
    std::size_t pair_index = 0;
    for (std::size_t i = 0; i < count; ++i) {
        for (std::size_t j = i + 1; j < count; ++j, ++pair_index) {
            scaled_points[2 * pair_index] =
                scale_public(ciphertexts[j].c3, inverses[pair_index]);
            scaled_points[2 * pair_index + 1] =
                scale_public(ciphertexts[i].c3, fr_neg(inverses[pair_index]));
        }
    }

    std::vector<blst::G1AffinePoint> scaled_affine(scaled_points.size());
    auto scaled_prefixes =
        affine_prefix_storage<blst::G1Point>(scaled_points.size());
    to_affine_batch(scaled_affine, scaled_points, scaled_prefixes);

    const std::size_t row_width = 2 * (count - 1);
    std::vector<blst::G1AffinePoint> row_points(count * row_width);
    std::vector<std::size_t> row_offsets(count, 0);

#if BLST_DEFAULT_AVX
    std::vector<blst::avx::PreparedLineTable> prepared_lines(count);
    for (std::size_t i = 0; i < count; ++i) {
        blst::avx::GT::precompute_lines(
            prepared_lines[i].lines, out.c2_affine[i]);
    }
    std::vector<const blst::avx::PreparedLine*> row_lines(
        count * row_width);
#else
    struct ScalarPreparedLineTable
    {
        blst::scalar::Fp6 lines[68];
    };
    std::vector<ScalarPreparedLineTable> prepared_lines(count);
    for (std::size_t i = 0; i < count; ++i) {
        blst::scalar::GT::precompute_lines(
            prepared_lines[i].lines, out.c2_affine[i]);
    }
    std::vector<const blst::scalar::Fp6*> row_lines(count * row_width);
#endif

    pair_index = 0;
    for (std::size_t i = 0; i < count; ++i) {
        for (std::size_t j = i + 1; j < count; ++j, ++pair_index) {
            for (const std::size_t row : { i, j }) {
                const std::size_t offset = row * row_width + row_offsets[row];
                row_points[offset] = scaled_affine[2 * pair_index];
                row_points[offset + 1] = scaled_affine[2 * pair_index + 1];
                row_lines[offset] = prepared_lines[i].lines;
                row_lines[offset + 1] = prepared_lines[j].lines;
                row_offsets[row] += 2;
            }
        }
    }

#if BLST_DEFAULT_AVX
    std::vector<blst::avx::MillerResult> millers(
        count, blst::avx::MillerResult::identity());
    for (std::size_t i = 0; i < count; ++i) {
        blst::avx::GT::miller_loop_lines_n_preprocessed(
            millers[i],
            row_lines.data() + i * row_width,
            row_points.data() + i * row_width,
            row_width);
    }
    blst::avx::GT::FinalExpScratch final_exp_scratch;
    blst::avx::GT::final_exp_n(
        out.cross_terms.data(),
        millers.data(),
        count,
        final_exp_scratch);
#else
    std::vector<blst::scalar::PreparedG1LineScale> row_scales(
        count * row_width);
    for (std::size_t i = 0; i < row_points.size(); ++i) {
        blst::scalar::Fp::add(
            row_scales[i].neg_x2, row_points[i].x, row_points[i].x);
        blst::scalar::Fp::neg(
            row_scales[i].neg_x2, row_scales[i].neg_x2);
        blst::scalar::Fp::add(
            row_scales[i].y2, row_points[i].y, row_points[i].y);
    }
    for (std::size_t i = 0; i < count; ++i) {
        blst::scalar::MillerResult miller =
            blst::scalar::MillerResult::identity();
        blst::scalar::GT::miller_loop_lines_n_preprocessed(
            miller,
            row_lines.data() + i * row_width,
            row_scales.data() + i * row_width,
            row_width);
        blst::scalar::GT::final_exp(out.cross_terms[i], miller);
    }
#endif
    return out;
}

BatchPrecomputation PrecomputeBatchDecryptionFast(
    const std::span<const Ciphertext> ciphertexts)
{
    std::vector<blst::Fr> points = ciphertext_points(ciphertexts);
    BatchPrecomputation out = make_precomputation_shell(ciphertexts);
    const std::size_t count = ciphertexts.size();
    if (count == 1) {
        out.cross_terms[0] = blst::GT::identity();
        return out;
    }

    PolynomialArithmetic arithmetic(next_power_of_two(2 * count));
    const std::unique_ptr<ProductTreeNode> tree =
        build_product_tree(arithmetic, points, 0, count);

    std::vector<blst::G2Point> c2;
    std::vector<blst::G1Point> c3;
    c2.reserve(count);
    c3.reserve(count);
    for (const Ciphertext& ciphertext : ciphertexts) {
        c2.push_back(ciphertext.c2);
        c3.push_back(ciphertext.c3);
    }
    const std::vector<blst::G2Point> n2 =
        build_numerator_polynomial<blst::G2Point>(
            arithmetic, *tree, std::span<const blst::G2Point>(c2));
    const std::vector<blst::G1Point> n3 =
        build_numerator_polynomial<blst::G1Point>(
            arithmetic, *tree, std::span<const blst::G1Point>(c3));
    const std::vector<blst::Fr> f_prime = derivative(tree->polynomial);
    const std::vector<blst::G2Point> n2_prime = derivative<blst::G2Point>(n2);
    const std::vector<blst::G1Point> n3_prime = derivative<blst::G1Point>(n3);

    std::vector<blst::Fr> f_prime_values(count, fr_zero());
    std::vector<blst::G2Point> n2_prime_values(count, blst::G2Point::infinity());
    std::vector<blst::G1Point> n3_prime_values(count, blst::G1Point::infinity());
    evaluate_scalar_recursive(arithmetic, *tree, f_prime, f_prime_values);
    evaluate_group_recursive<blst::G2Point>(
        arithmetic,
        *tree,
        std::span<const blst::G2Point>(n2_prime),
        std::span<blst::G2Point>(n2_prime_values));
    evaluate_group_recursive<blst::G1Point>(
        arithmetic,
        *tree,
        std::span<const blst::G1Point>(n3_prime),
        std::span<blst::G1Point>(n3_prime_values));
    const std::vector<blst::Fr> inverse_f_prime_values =
        batch_invert_nonzero(f_prime_values);

    compute_cross_terms_from_evaluations(
        out,
        ciphertexts,
        n2_prime_values,
        n3_prime_values,
        inverse_f_prime_values);
    return out;
}

std::vector<blst::GT> FinalizeBatchDecryption(
    const BatchPrecomputation& precomputation,
    const blst::G1Point& combined_share,
    const std::span<const Ciphertext> ciphertexts)
{
    if (precomputation.cross_terms.size() != ciphertexts.size() ||
        precomputation.c2_affine.size() != ciphertexts.size() ||
        ciphertexts.empty()) {
        throw std::invalid_argument("transparent BTE precomputation shape mismatch");
    }

    blst::G1AffinePoint combined_affine{};
    blst::G1Point::to_affine(combined_affine, combined_share);
    std::vector<blst::GT> messages(ciphertexts.size(), blst::GT::identity());
    for (std::size_t i = 0; i < ciphertexts.size(); ++i) {
        blst::MillerResult miller{};
        blst::GT::miller_loop(miller, &precomputation.c2_affine[i], &combined_affine);
        blst::GT mask{};
        blst::GT::final_exp(mask, miller);
        messages[i] = add_gt(
            sub_gt(ciphertexts[i].c4, mask),
            precomputation.cross_terms[i]);
    }
    return messages;
}
}
