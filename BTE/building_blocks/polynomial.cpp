#include "BTE/building_blocks/polynomial.h"

#include <algorithm>
#include <cstring>
#include <memory>
#include <stdexcept>

#include "blst/fft/context.h"
#include "blst/fft/ctFft.h"
#include "blst/fft/primitive_root.h"

namespace bte::building_blocks
{
    namespace
    {
        blst::Fr fr_zero()
        {
            blst::Fr zero{};
            std::memset(&zero, 0, sizeof(zero));
            return zero;
        }

        bool fr_is_zero(const blst::Fr& value)
        {
            return blst::Fr::equal(value, fr_zero());
        }

        blst::Fr fr_add(const blst::Fr& lhs, const blst::Fr& rhs)
        {
            blst::Fr out{};
            blst::Fr::add(out, lhs, rhs);
            return out;
        }

        blst::Fr fr_sub(const blst::Fr& lhs, const blst::Fr& rhs)
        {
            blst::Fr out{};
            blst::Fr::sub(out, lhs, rhs);
            return out;
        }

        blst::Fr fr_mul(const blst::Fr& lhs, const blst::Fr& rhs)
        {
            blst::Fr out{};
            blst::Fr::mul(out, lhs, rhs);
            return out;
        }

        blst::Fr fr_neg(const blst::Fr& value)
        {
            blst::Fr out{};
            blst::Fr::neg(out, value);
            return out;
        }

        blst::Fr fr_inverse(const blst::Fr& value)
        {
            blst::Fr out{};
            blst::Fr::reciprocal(out, value);
            return out;
        }

        std::size_t next_power_of_two(std::size_t value)
        {
            std::size_t out = 1;
            while (out < value) {
                out <<= 1;
            }
            return out;
        }

        std::vector<blst::Fr> trim_polynomial(std::vector<blst::Fr> polynomial)
        {
            while (!polynomial.empty() && fr_is_zero(polynomial.back())) {
                polynomial.pop_back();
            }
            return polynomial;
        }

        std::vector<blst::Fr> reverse_polynomial(
            const std::span<const blst::Fr> polynomial,
            const std::size_t size = 0)
        {
            const std::size_t out_size = size == 0 ? polynomial.size() : size;
            std::vector<blst::Fr> reversed(out_size, fr_zero());
            const std::size_t copy = std::min(polynomial.size(), out_size);
            for (std::size_t i = 0; i < copy; ++i) {
                reversed[i] = polynomial[polynomial.size() - 1 - i];
            }
            return reversed;
        }

        std::vector<blst::Fr> take_prefix(
            const std::span<const blst::Fr> polynomial,
            const std::size_t count)
        {
            const std::size_t actual = std::min(polynomial.size(), count);
            return std::vector<blst::Fr>(polynomial.begin(), polynomial.begin() + actual);
        }

        std::vector<blst::Fr> add_polynomials(
            const std::span<const blst::Fr> lhs,
            const std::span<const blst::Fr> rhs)
        {
            std::vector<blst::Fr> out(std::max(lhs.size(), rhs.size()), fr_zero());
            for (std::size_t i = 0; i < lhs.size(); ++i) {
                out[i] = fr_add(out[i], lhs[i]);
            }
            for (std::size_t i = 0; i < rhs.size(); ++i) {
                out[i] = fr_add(out[i], rhs[i]);
            }
            return trim_polynomial(std::move(out));
        }

        std::vector<blst::Fr> subtract_polynomials(
            const std::span<const blst::Fr> lhs,
            const std::span<const blst::Fr> rhs)
        {
            std::vector<blst::Fr> out(std::max(lhs.size(), rhs.size()), fr_zero());
            for (std::size_t i = 0; i < lhs.size(); ++i) {
                out[i] = fr_add(out[i], lhs[i]);
            }
            for (std::size_t i = 0; i < rhs.size(); ++i) {
                out[i] = fr_sub(out[i], rhs[i]);
            }
            return trim_polynomial(std::move(out));
        }

        std::vector<blst::Fr> naive_multiply(
            const std::span<const blst::Fr> lhs,
            const std::span<const blst::Fr> rhs)
        {
            if (lhs.empty() || rhs.empty()) {
                return {};
            }

            std::vector<blst::Fr> out(lhs.size() + rhs.size() - 1, fr_zero());
            for (std::size_t i = 0; i < lhs.size(); ++i) {
                for (std::size_t j = 0; j < rhs.size(); ++j) {
                    out[i + j] = fr_add(out[i + j], fr_mul(lhs[i], rhs[j]));
                }
            }
            return trim_polynomial(std::move(out));
        }

        std::vector<blst::Fr> make_fft_roots(const std::size_t n)
        {
            const blst::Fr prim = blst::primRootOfUnity<blst::Fr>(n);
            std::vector<blst::Fr> roots(n);
            blst::Fr current = blst::Fr::one();
            for (std::size_t i = 0; i < n; ++i) {
                roots[i] = current;
                blst::Fr::mul(current, current, prim);
            }
            return roots;
        }

        std::vector<blst::Fr> make_ifft_roots(const std::span<const blst::Fr> roots)
        {
            std::vector<blst::Fr> inv_roots(roots.size());
            for (std::size_t i = 0; i < roots.size(); ++i) {
                blst::Fr::reciprocal(inv_roots[i], roots[i]);
            }
            return inv_roots;
        }

        std::vector<blst::Fr> inverse_series(
            const std::span<const blst::Fr> polynomial,
            const std::size_t size)
        {
            if (size == 0) {
                return {};
            }
            if (polynomial.empty() || fr_is_zero(polynomial.front())) {
                throw std::invalid_argument("Power series inversion requires a nonzero constant term");
            }

            std::vector<blst::Fr> inverse(1, fr_inverse(polynomial.front()));
            std::size_t width = 1;
            while (width < size) {
                const std::size_t next_width = std::min(width << 1, size);
                const std::vector<blst::Fr> truncated = take_prefix(polynomial, next_width);
                std::vector<blst::Fr> product = multiply_polynomials(truncated, inverse);
                product.resize(next_width, fr_zero());
                product[0] = fr_sub(blst::Fr::from_uint64(2), product[0]);
                for (std::size_t i = 1; i < next_width; ++i) {
                    product[i] = fr_neg(product[i]);
                }
                inverse = multiply_polynomials(inverse, product);
                inverse.resize(next_width, fr_zero());
                width = next_width;
            }
            inverse.resize(size, fr_zero());
            return inverse;
        }

        blst::Fr evaluate_polynomial(
            const std::span<const blst::Fr> polynomial,
            const blst::Fr& point)
        {
            blst::Fr acc = fr_zero();
            for (auto it = polynomial.rbegin(); it != polynomial.rend(); ++it) {
                acc = fr_add(fr_mul(acc, point), *it);
            }
            return acc;
        }

        struct ProductTreeNode
        {
            std::size_t begin = 0;
            std::size_t end = 0;
            std::vector<blst::Fr> polynomial;
            std::unique_ptr<ProductTreeNode> left;
            std::unique_ptr<ProductTreeNode> right;
        };

        std::unique_ptr<ProductTreeNode> build_product_tree(
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

            const std::size_t mid = begin + ((end - begin) >> 1);
            node->left = build_product_tree(points, begin, mid);
            node->right = build_product_tree(points, mid, end);
            node->polynomial = multiply_polynomials(node->left->polynomial, node->right->polynomial);
            return node;
        }

        void evaluate_many_recursive(
            const ProductTreeNode& node,
            const std::span<const blst::Fr> polynomial,
            std::vector<blst::Fr>& out)
        {
            if (node.end - node.begin == 1) {
                out[node.begin] = polynomial.empty() ? fr_zero() : polynomial.front();
                return;
            }

            const auto [ignored_left_q, left_remainder] = divide_remainder(polynomial, node.left->polynomial);
            const auto [ignored_right_q, right_remainder] = divide_remainder(polynomial, node.right->polynomial);
            (void)ignored_left_q;
            (void)ignored_right_q;
            evaluate_many_recursive(*node.left, left_remainder, out);
            evaluate_many_recursive(*node.right, right_remainder, out);
        }

        void quotient_accumulate_recursive(
            const ProductTreeNode& node,
            const std::vector<blst::Fr>& accumulator,
            std::vector<std::vector<blst::Fr>>& out)
        {
            if (node.end - node.begin == 1) {
                out[node.begin] = trim_polynomial(accumulator);
                return;
            }

            std::vector<blst::Fr> left_accumulator = multiply_polynomials(accumulator, node.right->polynomial);
            std::vector<blst::Fr> right_accumulator = multiply_polynomials(accumulator, node.left->polynomial);
            quotient_accumulate_recursive(*node.left, left_accumulator, out);
            quotient_accumulate_recursive(*node.right, right_accumulator, out);
        }
    }

    std::vector<blst::Fr> multiply_polynomials(
        const std::span<const blst::Fr> lhs,
        const std::span<const blst::Fr> rhs)
    {
        if (lhs.empty() || rhs.empty()) {
            return {};
        }

        constexpr std::size_t naive_threshold = 32;
        if (std::min(lhs.size(), rhs.size()) <= naive_threshold) {
            return naive_multiply(lhs, rhs);
        }

        const std::size_t result_size = lhs.size() + rhs.size() - 1;
        const std::size_t fft_size = next_power_of_two(result_size);
        std::vector<blst::Fr> a(fft_size, fr_zero());
        std::vector<blst::Fr> b(fft_size, fr_zero());
        std::vector<blst::Fr> c(fft_size, fr_zero());

        std::copy(lhs.begin(), lhs.end(), a.begin());
        std::copy(rhs.begin(), rhs.end(), b.begin());

        const std::vector<blst::Fr> roots = make_fft_roots(fft_size);
        const std::vector<blst::Fr> inv_roots = make_ifft_roots(roots);
        poly_mult_fft<blst::Fr, blst::Fr, blst::Fr, blst::Fr, blst::FrCoeffCtx>(
            a,
            b,
            false,
            c,
            roots,
            inv_roots,
            nullptr,
            lhs.size(),
            rhs.size());

        c.resize(result_size);
        return trim_polynomial(std::move(c));
    }

    std::vector<blst::Fr> derivative_polynomial(
        const std::span<const blst::Fr> polynomial)
    {
        if (polynomial.size() <= 1) {
            return {};
        }

        std::vector<blst::Fr> derivative(polynomial.size() - 1, fr_zero());
        for (std::size_t i = 1; i < polynomial.size(); ++i) {
            derivative[i - 1] = fr_mul(polynomial[i], blst::Fr::from_uint64(i));
        }
        return trim_polynomial(std::move(derivative));
    }

    std::pair<std::vector<blst::Fr>, std::vector<blst::Fr>> divide_remainder(
        const std::span<const blst::Fr> dividend,
        const std::span<const blst::Fr> divisor)
    {
        if (divisor.empty()) {
            throw std::invalid_argument("Polynomial division requires a non-empty divisor");
        }

        if (dividend.size() < divisor.size()) {
            return { {}, trim_polynomial(std::vector<blst::Fr>(dividend.begin(), dividend.end())) };
        }

        constexpr std::size_t naive_threshold = 32;
        if (divisor.size() <= naive_threshold) {
            std::vector<blst::Fr> remainder(dividend.begin(), dividend.end());
            std::vector<blst::Fr> quotient(dividend.size() - divisor.size() + 1, fr_zero());
            const blst::Fr divisor_lead_inverse = fr_inverse(divisor.back());

            for (std::size_t idx = quotient.size(); idx-- > 0;) {
                const std::size_t degree = idx + divisor.size() - 1;
                if (fr_is_zero(remainder[degree])) {
                    continue;
                }

                quotient[idx] = fr_mul(remainder[degree], divisor_lead_inverse);
                for (std::size_t j = 0; j < divisor.size(); ++j) {
                    remainder[idx + j] = fr_sub(
                        remainder[idx + j],
                        fr_mul(quotient[idx], divisor[j]));
                }
            }

            remainder.resize(divisor.size() - 1);
            return {
                trim_polynomial(std::move(quotient)),
                trim_polynomial(std::move(remainder))
            };
        }

        const std::size_t quotient_size = dividend.size() - divisor.size() + 1;
        const std::vector<blst::Fr> reversed_dividend = reverse_polynomial(dividend);
        const std::vector<blst::Fr> reversed_divisor = reverse_polynomial(divisor);
        const std::vector<blst::Fr> reversed_dividend_prefix =
            take_prefix(reversed_dividend, quotient_size);
        const std::vector<blst::Fr> inverse_reversed_divisor =
            inverse_series(reversed_divisor, quotient_size);

        std::vector<blst::Fr> quotient_reversed =
            multiply_polynomials(reversed_dividend_prefix, inverse_reversed_divisor);
        quotient_reversed.resize(quotient_size, fr_zero());

        std::vector<blst::Fr> quotient = reverse_polynomial(quotient_reversed);
        quotient = trim_polynomial(std::move(quotient));

        const std::vector<blst::Fr> product = multiply_polynomials(quotient, divisor);
        std::vector<blst::Fr> remainder = subtract_polynomials(dividend, product);
        if (remainder.size() >= divisor.size()) {
            remainder.resize(divisor.size() - 1);
            remainder = trim_polynomial(std::move(remainder));
        }

        return { quotient, remainder };
    }

    std::vector<blst::Fr> evaluate_many(
        const std::span<const blst::Fr> polynomial,
        const std::span<const blst::Fr> points)
    {
        if (points.empty()) {
            return {};
        }

        if (points.size() == 1) {
            return { evaluate_polynomial(polynomial, points.front()) };
        }

        const std::unique_ptr<ProductTreeNode> tree = build_product_tree(points, 0, points.size());
        const auto [ignored_quotient, remainder] = divide_remainder(polynomial, tree->polynomial);
        (void)ignored_quotient;

        std::vector<blst::Fr> values(points.size(), fr_zero());
        evaluate_many_recursive(*tree, remainder, values);
        return values;
    }

    std::vector<std::vector<blst::Fr>> exclusion_quotients(
        const std::span<const blst::Fr> roots)
    {
        if (roots.empty()) {
            return {};
        }

        if (roots.size() == 1) {
            return { { blst::Fr::one() } };
        }

        const std::unique_ptr<ProductTreeNode> tree = build_product_tree(roots, 0, roots.size());
        std::vector<std::vector<blst::Fr>> quotients(roots.size());
        quotient_accumulate_recursive(*tree, { blst::Fr::one() }, quotients);
        return quotients;
    }

    std::vector<blst::Fr> build_vanishing_polynomial(const std::vector<blst::Fr>& roots)
    {
        if (roots.empty()) {
            return { blst::Fr::one() };
        }
        return build_product_tree(std::span<const blst::Fr>(roots.data(), roots.size()), 0, roots.size())->polynomial;
    }

    std::vector<blst::Fr> synthetic_divide(
        const std::vector<blst::Fr>& polynomial,
        const blst::Fr& root,
        blst::Fr* const remainder)
    {
        if (polynomial.size() <= 1) {
            if (remainder != nullptr) {
                *remainder = polynomial.empty() ? fr_zero() : polynomial.front();
            }
            return {};
        }

        const std::size_t degree = polynomial.size() - 1;
        std::vector<blst::Fr> quotient(degree, fr_zero());
        quotient[degree - 1] = polynomial[degree];
        for (std::size_t k = degree - 1; k > 0; --k) {
            blst::Fr term{};
            blst::Fr::mul(term, quotient[k], root);
            blst::Fr::add(quotient[k - 1], polynomial[k], term);
        }

        if (remainder != nullptr) {
            blst::Fr term{};
            blst::Fr::mul(term, quotient.front(), root);
            blst::Fr::add(*remainder, polynomial.front(), term);
        }

        return quotient;
    }
}
