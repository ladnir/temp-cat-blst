#include "blst/fft/fr_interpolation.h"

#include "blst/fft/context.h"
#include "blst/fft/ctFft.h"
#include "blst/fft/domain.h"

#include <algorithm>
#include <bit>
#include <limits>
#include <memory>
#include <new>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace blst::fr_interpolation
{
namespace
{
    struct DuplicatePointError
    {};

    bool checked_add(
        const std::size_t lhs,
        const std::size_t rhs,
        std::size_t& out) noexcept
    {
        if (rhs > std::numeric_limits<std::size_t>::max() - lhs) {
            return false;
        }
        out = lhs + rhs;
        return true;
    }

    std::size_t polynomial_product_size(
        const std::size_t lhs_size,
        const std::size_t rhs_size)
    {
        std::size_t sum = 0;
        if (!checked_add(lhs_size, rhs_size, sum)) {
            throw std::length_error("polynomial product size overflow");
        }
        return sum - 1;
    }

    std::size_t next_power_of_two(const std::size_t value) noexcept
    {
        if (value <= 1) {
            return 1;
        }

        const unsigned int shift = std::bit_width(value - 1);
        if (shift >= std::numeric_limits<std::size_t>::digits) {
            return 0;
        }

        const std::size_t out = std::size_t{1} << shift;
        return fft_detail::valid_size(out) ? out : 0;
    }

    Fr fr_zero() noexcept
    {
        return Fr{};
    }

    bool fr_is_zero(const Fr& value) noexcept
    {
        return Fr::equal(value, fr_zero());
    }

    Fr fr_add(const Fr& lhs, const Fr& rhs) noexcept
    {
        Fr out{};
        Fr::add(out, lhs, rhs);
        return out;
    }

    Fr fr_sub(const Fr& lhs, const Fr& rhs) noexcept
    {
        Fr out{};
        Fr::sub(out, lhs, rhs);
        return out;
    }

    Fr fr_mul(const Fr& lhs, const Fr& rhs) noexcept
    {
        Fr out{};
        Fr::mul(out, lhs, rhs);
        return out;
    }

    Fr fr_neg(const Fr& value) noexcept
    {
        Fr out{};
        Fr::neg(out, value);
        return out;
    }

    Fr fr_inverse(const Fr& value) noexcept
    {
        Fr out{};
        Fr::reciprocal(out, value);
        return out;
    }

    std::vector<Fr> make_roots(const std::size_t n, const bool inverse)
    {
        return fft_detail::make_roots<Fr, FrCoeffCtx>(n, inverse);
    }

    std::vector<Fr> trim_polynomial(std::vector<Fr> polynomial)
    {
        while (!polynomial.empty() && fr_is_zero(polynomial.back())) {
            polynomial.pop_back();
        }
        return polynomial;
    }

    std::vector<Fr> take_prefix(
        const std::span<const Fr> polynomial,
        const std::size_t count)
    {
        const std::size_t actual = std::min(polynomial.size(), count);
        return std::vector<Fr>(polynomial.begin(), polynomial.begin() + actual);
    }

    std::vector<Fr> reverse_polynomial(
        const std::span<const Fr> polynomial,
        const std::size_t size = 0)
    {
        const std::size_t out_size = size == 0 ? polynomial.size() : size;
        std::vector<Fr> reversed(out_size, fr_zero());
        const std::size_t copy = std::min(polynomial.size(), out_size);
        for (std::size_t i = 0; i < copy; ++i) {
            reversed[i] = polynomial[polynomial.size() - 1 - i];
        }
        return reversed;
    }

    std::vector<Fr> add_polynomials(
        const std::span<const Fr> lhs,
        const std::span<const Fr> rhs)
    {
        std::vector<Fr> out(std::max(lhs.size(), rhs.size()), fr_zero());
        for (std::size_t i = 0; i < lhs.size(); ++i) {
            out[i] = fr_add(out[i], lhs[i]);
        }
        for (std::size_t i = 0; i < rhs.size(); ++i) {
            out[i] = fr_add(out[i], rhs[i]);
        }
        return trim_polynomial(std::move(out));
    }

    std::vector<Fr> subtract_polynomials(
        const std::span<const Fr> lhs,
        const std::span<const Fr> rhs)
    {
        std::vector<Fr> out(std::max(lhs.size(), rhs.size()), fr_zero());
        for (std::size_t i = 0; i < lhs.size(); ++i) {
            out[i] = fr_add(out[i], lhs[i]);
        }
        for (std::size_t i = 0; i < rhs.size(); ++i) {
            out[i] = fr_sub(out[i], rhs[i]);
        }
        return trim_polynomial(std::move(out));
    }

    std::vector<Fr> naive_multiply(
        const std::span<const Fr> lhs,
        const std::span<const Fr> rhs)
    {
        if (lhs.empty() || rhs.empty()) {
            return {};
        }

        const std::size_t result_size = polynomial_product_size(lhs.size(), rhs.size());
        std::vector<Fr> out(result_size, fr_zero());
        for (std::size_t i = 0; i < lhs.size(); ++i) {
            for (std::size_t j = 0; j < rhs.size(); ++j) {
                out[i + j] = fr_add(out[i + j], fr_mul(lhs[i], rhs[j]));
            }
        }
        return trim_polynomial(std::move(out));
    }

    std::vector<Fr> multiply_polynomials(
        const std::span<const Fr> lhs,
        const std::span<const Fr> rhs)
    {
        if (lhs.empty() || rhs.empty()) {
            return {};
        }

        constexpr std::size_t naive_threshold = 32;
        if (std::min(lhs.size(), rhs.size()) <= naive_threshold) {
            return naive_multiply(lhs, rhs);
        }

        const std::size_t result_size = polynomial_product_size(lhs.size(), rhs.size());
        const std::size_t fft_size = next_power_of_two(result_size);
        if (fft_size == 0) {
            throw std::length_error("polynomial product exceeds supported FFT size");
        }

        std::vector<Fr> a(fft_size, fr_zero());
        std::vector<Fr> b(fft_size, fr_zero());
        std::vector<Fr> c(fft_size, fr_zero());
        std::copy(lhs.begin(), lhs.end(), a.begin());
        std::copy(rhs.begin(), rhs.end(), b.begin());

        const std::vector<Fr> roots = make_roots(fft_size, false);
        const std::vector<Fr> inv_roots = make_roots(fft_size, true);
        poly_mult_fft<Fr, Fr, Fr, Fr, FrCoeffCtx>(
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

    std::vector<Fr> derivative_polynomial(const std::span<const Fr> polynomial)
    {
        if (polynomial.size() <= 1) {
            return {};
        }

        std::vector<Fr> derivative(polynomial.size() - 1, fr_zero());
        for (std::size_t i = 1; i < polynomial.size(); ++i) {
            derivative[i - 1] = fr_mul(polynomial[i], Fr::from_uint64(i));
        }
        return trim_polynomial(std::move(derivative));
    }

    std::vector<Fr> inverse_series(
        const std::span<const Fr> polynomial,
        const std::size_t size)
    {
        if (size == 0) {
            return {};
        }
        if (polynomial.empty() || fr_is_zero(polynomial.front())) {
            throw std::invalid_argument("power series inverse requires nonzero constant");
        }

        std::vector<Fr> inverse(1, fr_inverse(polynomial.front()));
        std::size_t width = 1;
        while (width < size) {
            const std::size_t remaining = size - width;
            const std::size_t next_width = width > remaining ? size : width + width;
            const std::vector<Fr> truncated = take_prefix(polynomial, next_width);
            std::vector<Fr> product = multiply_polynomials(truncated, inverse);
            product.resize(next_width, fr_zero());
            product[0] = fr_sub(Fr::from_uint64(2), product[0]);
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

    std::pair<std::vector<Fr>, std::vector<Fr>> divide_remainder(
        const std::span<const Fr> dividend,
        const std::span<const Fr> divisor)
    {
        if (divisor.empty()) {
            throw std::invalid_argument("polynomial division requires divisor");
        }

        if (dividend.size() < divisor.size()) {
            return { {}, trim_polynomial(std::vector<Fr>(dividend.begin(), dividend.end())) };
        }

        constexpr std::size_t naive_threshold = 32;
        if (divisor.size() <= naive_threshold) {
            std::vector<Fr> remainder(dividend.begin(), dividend.end());
            std::vector<Fr> quotient(dividend.size() - divisor.size() + 1, fr_zero());
            const Fr divisor_lead_inverse = fr_inverse(divisor.back());

            for (std::size_t idx = quotient.size(); idx-- > 0;) {
                const std::size_t degree = idx + divisor.size() - 1;
                if (fr_is_zero(remainder[degree])) {
                    continue;
                }

                quotient[idx] = fr_mul(remainder[degree], divisor_lead_inverse);
                for (std::size_t j = 0; j < divisor.size(); ++j) {
                    remainder[idx + j] =
                        fr_sub(remainder[idx + j], fr_mul(quotient[idx], divisor[j]));
                }
            }

            remainder.resize(divisor.size() - 1);
            return {
                trim_polynomial(std::move(quotient)),
                trim_polynomial(std::move(remainder))
            };
        }

        const std::size_t quotient_size = dividend.size() - divisor.size() + 1;
        const std::vector<Fr> reversed_dividend = reverse_polynomial(dividend);
        const std::vector<Fr> reversed_divisor = reverse_polynomial(divisor);
        const std::vector<Fr> reversed_prefix =
            take_prefix(reversed_dividend, quotient_size);
        const std::vector<Fr> inverse_reversed_divisor =
            inverse_series(reversed_divisor, quotient_size);

        std::vector<Fr> quotient_reversed =
            multiply_polynomials(reversed_prefix, inverse_reversed_divisor);
        quotient_reversed.resize(quotient_size, fr_zero());

        std::vector<Fr> quotient =
            trim_polynomial(reverse_polynomial(quotient_reversed));
        const std::vector<Fr> product = multiply_polynomials(quotient, divisor);
        std::vector<Fr> remainder = subtract_polynomials(dividend, product);
        if (remainder.size() >= divisor.size()) {
            remainder.resize(divisor.size() - 1);
            remainder = trim_polynomial(std::move(remainder));
        }

        return { quotient, remainder };
    }

    struct ProductTreeNode
    {
        std::size_t begin = 0;
        std::size_t end = 0;
        std::vector<Fr> polynomial;
        std::unique_ptr<ProductTreeNode> left;
        std::unique_ptr<ProductTreeNode> right;
    };

    std::unique_ptr<ProductTreeNode> build_product_tree(
        const std::span<const Fr> points,
        const std::size_t begin,
        const std::size_t end)
    {
        auto node = std::make_unique<ProductTreeNode>();
        node->begin = begin;
        node->end = end;

        if (end - begin == 1) {
            node->polynomial = { fr_neg(points[begin]), Fr::one() };
            return node;
        }

        const std::size_t mid = begin + ((end - begin) >> 1);
        node->left = build_product_tree(points, begin, mid);
        node->right = build_product_tree(points, mid, end);
        node->polynomial =
            multiply_polynomials(node->left->polynomial, node->right->polynomial);
        return node;
    }

    void evaluate_many_recursive(
        const ProductTreeNode& node,
        const std::span<const Fr> polynomial,
        std::vector<Fr>& out)
    {
        if (node.end - node.begin == 1) {
            out[node.begin] = polynomial.empty() ? fr_zero() : polynomial.front();
            return;
        }

        const auto [ignored_left_q, left_remainder] =
            divide_remainder(polynomial, node.left->polynomial);
        const auto [ignored_right_q, right_remainder] =
            divide_remainder(polynomial, node.right->polynomial);
        (void)ignored_left_q;
        (void)ignored_right_q;
        evaluate_many_recursive(*node.left, left_remainder, out);
        evaluate_many_recursive(*node.right, right_remainder, out);
    }

    std::vector<Fr> evaluate_many_with_tree(
        const ProductTreeNode& tree,
        const std::span<const Fr> polynomial)
    {
        const auto [ignored_quotient, remainder] =
            divide_remainder(polynomial, tree.polynomial);
        (void)ignored_quotient;

        std::vector<Fr> values(tree.end - tree.begin, fr_zero());
        evaluate_many_recursive(tree, remainder, values);
        return values;
    }

    std::vector<Fr> batch_invert_nonzero(const std::span<const Fr> values)
    {
        std::size_t prefix_count = 0;
        if (!checked_add(values.size(), 1, prefix_count)) {
            throw std::length_error("batch inversion size overflow");
        }

        std::vector<Fr> prefixes(prefix_count, Fr::one());
        for (std::size_t i = 0; i < values.size(); ++i) {
            if (fr_is_zero(values[i])) {
                throw DuplicatePointError{};
            }
            prefixes[i + 1] = fr_mul(prefixes[i], values[i]);
        }

        Fr all_inverse = fr_inverse(prefixes.back());
        std::vector<Fr> inverses(values.size(), Fr::one());
        for (std::size_t i = values.size(); i-- > 0;) {
            inverses[i] = fr_mul(all_inverse, prefixes[i]);
            all_inverse = fr_mul(all_inverse, values[i]);
        }
        return inverses;
    }

    void batch_invert_nonzero_into(
        const std::span<const Fr> values,
        const std::span<Fr> prefixes,
        const std::span<Fr> inverses)
    {
        if (prefixes.size() != values.size() + 1 ||
            inverses.size() != values.size()) {
            throw std::invalid_argument("batch inversion scratch shape mismatch");
        }

        prefixes[0] = Fr::one();
        for (std::size_t i = 0; i < values.size(); ++i) {
            if (fr_is_zero(values[i])) {
                throw DuplicatePointError{};
            }
            prefixes[i + 1] = fr_mul(prefixes[i], values[i]);
        }

        Fr all_inverse = fr_inverse(prefixes.back());
        for (std::size_t i = values.size(); i-- > 0;) {
            inverses[i] = fr_mul(all_inverse, prefixes[i]);
            all_inverse = fr_mul(all_inverse, values[i]);
        }
    }

    std::vector<Fr> interpolate_recursive(
        const ProductTreeNode& node,
        const std::span<const Fr> weights)
    {
        if (node.end - node.begin == 1) {
            return { weights[node.begin] };
        }

        std::vector<Fr> left = interpolate_recursive(*node.left, weights);
        std::vector<Fr> right = interpolate_recursive(*node.right, weights);
        left = multiply_polynomials(left, node.right->polynomial);
        right = multiply_polynomials(right, node.left->polynomial);
        return add_polynomials(left, right);
    }

    void interpolate_polynomial(
        const std::span<const Fr> points,
        const std::span<const Fr> values,
        const std::span<Fr> out)
    {
        const std::unique_ptr<ProductTreeNode> tree =
            build_product_tree(points, 0, points.size());
        const std::vector<Fr> derivative = derivative_polynomial(tree->polynomial);
        const std::vector<Fr> denominators =
            evaluate_many_with_tree(*tree, derivative);
        const std::vector<Fr> denominator_inverses =
            batch_invert_nonzero(denominators);

        std::vector<Fr> weights(points.size(), fr_zero());
        for (std::size_t i = 0; i < points.size(); ++i) {
            weights[i] = fr_mul(values[i], denominator_inverses[i]);
        }

        std::vector<Fr> coefficients = interpolate_recursive(*tree, weights);
        coefficients.resize(points.size(), fr_zero());
        std::copy(coefficients.begin(), coefficients.end(), out.begin());
    }

    Status exception_to_status() noexcept
    {
        try {
            throw;
        } catch (const DuplicatePointError&) {
            return Status::duplicate_point;
        } catch (const std::bad_alloc&) {
            return Status::allocation_failure;
        } catch (const std::length_error&) {
            return Status::invalid_size;
        } catch (...) {
            return Status::internal_failure;
        }
    }
}

Status interpolate(
    const std::span<Fr> out_coefficients,
    const std::span<const Fr> points,
    const std::span<const Fr> values) noexcept
{
    if (points.empty() || points.size() != values.size() ||
        points.size() != out_coefficients.size()) {
        return Status::invalid_shape;
    }

    try {
        interpolate_polynomial(points, values, out_coefficients);
        return Status::success;
    } catch (...) {
        return exception_to_status();
    }
}

Status interpolate_domain_suffix_with_scratch(
    const std::span<Fr> out_suffix,
    const std::span<const std::size_t> indices,
    const std::span<const Fr> values,
    const std::size_t ell,
    const std::size_t fft_len,
    const std::span<Fr> scratch) noexcept
{
    if (!fft_detail::valid_size(fft_len)) {
        return Status::invalid_size;
    }

    const std::size_t n = indices.size();
    if (n == 0 || n != values.size() || n > fft_len || ell > n ||
        out_suffix.size() != n - ell) {
        return Status::invalid_shape;
    }

    const std::size_t required = domain_suffix_scratch_len(n, fft_len);
    if (required == 0) {
        return Status::invalid_size;
    }
    if (scratch.size() < required) {
        return Status::insufficient_scratch;
    }
    for (const std::size_t index : indices) {
        if (index >= fft_len) {
            return Status::invalid_index;
        }
    }

    try {
        Fr* cursor = scratch.data();
        const std::span<Fr> roots(cursor, fft_len);
        cursor += fft_len;
        const std::span<Fr> inv_roots(cursor, fft_len);
        cursor += fft_len;
        const std::span<Fr> product(cursor, n + 1);
        cursor += n + 1;
        const std::span<Fr> domain_work(cursor, fft_len);
        cursor += fft_len;
        const std::span<Fr> denominators(cursor, n);
        cursor += n;
        const std::span<Fr> prefixes(cursor, n + 1);

        fft_detail::fill_root_pair<Fr, FrCoeffCtx>(roots, inv_roots);

        std::fill(product.begin(), product.end(), fr_zero());
        product[0] = Fr::one();
        std::size_t degree = 0;
        for (std::size_t sample = 0; sample < n; ++sample) {
            const Fr& point = roots[indices[sample]];
            product[degree + 1] = product[degree];
            for (std::size_t coeff = degree; coeff > 0; --coeff) {
                product[coeff] =
                    fr_sub(product[coeff - 1], fr_mul(point, product[coeff]));
            }
            product[0] = fr_neg(fr_mul(point, product[0]));
            ++degree;
        }

        std::fill(domain_work.begin(), domain_work.end(), fr_zero());
        for (std::size_t coeff = 1; coeff <= n; ++coeff) {
            domain_work[coeff - 1] =
                fr_mul(product[coeff], Fr::from_uint64(coeff));
        }
        fft<Fr, Fr, FrCoeffCtx>(domain_work, roots, false);

        for (std::size_t sample = 0; sample < n; ++sample) {
            const Fr denominator = domain_work[indices[sample]];
            if (fr_is_zero(denominator)) {
                return Status::duplicate_point;
            }
            denominators[sample] = denominator;
        }
        batch_invert_nonzero_into(denominators, prefixes, prefixes.subspan(0, n));

        std::fill(domain_work.begin(), domain_work.end(), fr_zero());
        const Fr fft_len_scalar = Fr::from_uint64(fft_len);
        for (std::size_t sample = 0; sample < n; ++sample) {
            const std::size_t index = indices[sample];
            const Fr& inverse_root = roots[(fft_len - index) % fft_len];
            const Fr domain_derivative = fr_mul(fft_len_scalar, inverse_root);
            const Fr weighted_value = fr_mul(values[sample], prefixes[sample]);
            domain_work[index] = fr_mul(weighted_value, domain_derivative);
        }
        fft<Fr, Fr, FrCoeffCtx>(domain_work, inv_roots, true);

        const std::size_t suffix_len = n - ell;
        for (std::size_t suffix_index = 0; suffix_index < suffix_len; ++suffix_index) {
            const std::size_t target_degree = fft_len + ell + suffix_index;
            const std::size_t min_product_index = ell + suffix_index + 1;
            Fr acc = fr_zero();
            for (std::size_t product_index = min_product_index;
                 product_index <= n;
                 ++product_index) {
                const std::size_t quotient_index = target_degree - product_index;
                acc = fr_add(
                    acc,
                    fr_mul(product[product_index], domain_work[quotient_index]));
            }
            out_suffix[suffix_index] = acc;
        }

        return Status::success;
    } catch (...) {
        return exception_to_status();
    }
}
}
