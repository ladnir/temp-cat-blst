// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs
//
// CatBLST scalar-field interpolation implementation.

#include "blst/fft/fr_interpolation.h"

#include "blst/common/defines.h"
#include "blst/fft/context.h"
#include "blst/fft/ctFft.h"
#include "blst/fft/domain.h"

#include <algorithm>
#include <bit>
#include <limits>
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

    void naive_multiply_into(
        const std::span<const Fr> lhs,
        const std::span<const Fr> rhs,
        const std::span<Fr> out)
    {
        const std::size_t result_size = polynomial_product_size(lhs.size(), rhs.size());
        if (out.size() != result_size) {
            throw std::invalid_argument("polynomial product output shape mismatch");
        }
        std::fill(out.begin(), out.end(), fr_zero());
        for (std::size_t i = 0; i < lhs.size(); ++i) {
            for (std::size_t j = 0; j < rhs.size(); ++j) {
                out[i + j] = fr_add(out[i + j], fr_mul(lhs[i], rhs[j]));
            }
        }
    }

    void multiply_polynomials_into(
        const std::span<const Fr> lhs,
        const std::span<const Fr> rhs,
        const std::span<Fr> out)
    {
        if (lhs.empty() || rhs.empty()) {
            if (!out.empty()) {
                throw std::invalid_argument("empty polynomial product output mismatch");
            }
            return;
        }

        const std::size_t result_size = polynomial_product_size(lhs.size(), rhs.size());
        if (out.size() != result_size) {
            throw std::invalid_argument("polynomial product output shape mismatch");
        }
        constexpr std::size_t naive_threshold = 32;
        if (std::min(lhs.size(), rhs.size()) <= naive_threshold) {
            naive_multiply_into(lhs, rhs, out);
            return;
        }

        const std::size_t fft_size = next_power_of_two(result_size);
        if (fft_size == 0) {
            throw std::length_error("polynomial product exceeds supported FFT size");
        }

        std::vector<Fr> a(fft_size, fr_zero());
        std::vector<Fr> b(fft_size, fr_zero());
        std::copy(lhs.begin(), lhs.end(), a.begin());
        std::copy(rhs.begin(), rhs.end(), b.begin());

        const std::vector<Fr> roots = make_roots(fft_size, false);
        fft<Fr, Fr, FrCoeffCtx>(a, roots, false);
        fft<Fr, Fr, FrCoeffCtx>(b, roots, false);
        for (std::size_t i = 0; i < fft_size; ++i) {
            a[i] = fr_mul(a[i], b[i]);
        }
        const Fr inverse_size = fr_inverse(Fr::from_uint64(fft_size));
        fft_with_forward_roots<Fr, Fr, FrCoeffCtx>(a, roots, true, inverse_size);
        std::copy_n(a.begin(), result_size, out.begin());
    }

    std::vector<Fr> multiply_polynomials_exact(
        const std::span<const Fr> lhs,
        const std::span<const Fr> rhs)
    {
        if (lhs.empty() || rhs.empty()) {
            return {};
        }
        std::vector<Fr> out(polynomial_product_size(lhs.size(), rhs.size()));
        multiply_polynomials_into(lhs, rhs, out);
        return out;
    }

    std::vector<Fr> multiply_polynomials(
        const std::span<const Fr> lhs,
        const std::span<const Fr> rhs)
    {
        return trim_polynomial(multiply_polynomials_exact(lhs, rhs));
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

    std::vector<Fr> polynomial_remainder(
        const std::span<const Fr> dividend,
        const std::span<const Fr> divisor)
    {
        if (divisor.empty()) {
            throw std::invalid_argument("polynomial division requires divisor");
        }

        if (dividend.size() < divisor.size()) {
            return trim_polynomial(std::vector<Fr>(dividend.begin(), dividend.end()));
        }

        constexpr std::size_t naive_threshold = 32;
        if (divisor.size() <= naive_threshold) {
            std::vector<Fr> remainder(dividend.begin(), dividend.end());
            const Fr divisor_lead_inverse = fr_inverse(divisor.back());
            const std::size_t quotient_size = dividend.size() - divisor.size() + 1;

            for (std::size_t idx = quotient_size; idx-- > 0;) {
                const std::size_t degree = idx + divisor.size() - 1;
                if (fr_is_zero(remainder[degree])) {
                    continue;
                }

                const Fr factor = fr_mul(remainder[degree], divisor_lead_inverse);
                for (std::size_t j = 0; j < divisor.size(); ++j) {
                    remainder[idx + j] =
                        fr_sub(remainder[idx + j], fr_mul(factor, divisor[j]));
                }
            }

            remainder.resize(divisor.size() - 1);
            return trim_polynomial(std::move(remainder));
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
            BLST_ASSERT(remainder.size() < divisor.size());
            throw std::logic_error("polynomial remainder degree invariant failed");
        }

        return remainder;
    }

    struct ProductTree
    {
        struct Node
        {
            static constexpr std::size_t no_child =
                std::numeric_limits<std::size_t>::max();

            std::size_t begin = 0;
            std::size_t end = 0;
            std::size_t polynomial_offset = 0;
            std::size_t polynomial_size = 0;
            std::size_t left = no_child;
            std::size_t right = no_child;
        };

        std::vector<Node> nodes;
        std::vector<Fr> coefficients;
        std::size_t root = Node::no_child;

        [[nodiscard]] std::span<const Fr> polynomial(
            const std::size_t node_index) const noexcept
        {
            const Node& node = nodes[node_index];
            return std::span<const Fr>(coefficients).subspan(
                node.polynomial_offset,
                node.polynomial_size);
        }
    };

    std::size_t product_tree_coefficient_count(const std::size_t point_count)
    {
        if (point_count == 1) {
            return 2;
        }
        const std::size_t left_count = point_count >> 1;
        const std::size_t right_count = point_count - left_count;
        std::size_t child_coefficients = 0;
        if (!checked_add(
                product_tree_coefficient_count(left_count),
                product_tree_coefficient_count(right_count),
                child_coefficients)) {
            throw std::length_error("product-tree coefficient size overflow");
        }
        std::size_t node_coefficients = 0;
        std::size_t total_coefficients = 0;
        if (!checked_add(point_count, 1, node_coefficients) ||
            !checked_add(
                child_coefficients,
                node_coefficients,
                total_coefficients)) {
            throw std::length_error("product-tree coefficient size overflow");
        }
        return total_coefficients;
    }

    std::size_t build_product_tree_node(
        ProductTree& tree,
        const std::span<const Fr> points,
        const std::size_t begin,
        const std::size_t end)
    {
        const std::size_t node_index = tree.nodes.size();
        tree.nodes.emplace_back();
        ProductTree::Node node;
        node.begin = begin;
        node.end = end;

        if (end - begin == 1) {
            node.polynomial_offset = tree.coefficients.size();
            node.polynomial_size = 2;
            tree.coefficients.push_back(fr_neg(points[begin]));
            tree.coefficients.push_back(Fr::one());
            tree.nodes[node_index] = node;
            return node_index;
        }

        const std::size_t mid = begin + ((end - begin) >> 1);
        node.left = build_product_tree_node(tree, points, begin, mid);
        node.right = build_product_tree_node(tree, points, mid, end);
        node.polynomial_offset = tree.coefficients.size();
        node.polynomial_size = end - begin + 1;
        tree.coefficients.resize(
            node.polynomial_offset + node.polynomial_size,
            fr_zero());
        multiply_polynomials_into(
            tree.polynomial(node.left),
            tree.polynomial(node.right),
            std::span<Fr>(tree.coefficients).subspan(
                node.polynomial_offset,
                node.polynomial_size));
        tree.nodes[node_index] = node;
        return node_index;
    }

    ProductTree build_product_tree(const std::span<const Fr> points)
    {
        std::size_t node_count = 0;
        if (!checked_add(points.size(), points.size(), node_count)) {
            throw std::length_error("product-tree node count overflow");
        }
        --node_count;

        // One topology allocation and one coefficient allocation replace the
        // former owning allocation at every node. Polynomial spans are recovered
        // from offsets only after any arena growth, so no pointer is retained.
        ProductTree tree;
        tree.nodes.reserve(node_count);
        const std::size_t coefficient_count =
            product_tree_coefficient_count(points.size());
        tree.coefficients.reserve(coefficient_count);
        tree.root = build_product_tree_node(tree, points, 0, points.size());
        BLST_ASSERT(tree.nodes.size() == node_count);
        BLST_ASSERT(tree.coefficients.size() == coefficient_count);
        return tree;
    }

    void evaluate_many_recursive(
        const ProductTree& tree,
        const std::size_t node_index,
        const std::span<const Fr> polynomial,
        std::vector<Fr>& out)
    {
        const ProductTree::Node& node = tree.nodes[node_index];
        if (node.end - node.begin == 1) {
            out[node.begin] = polynomial.empty() ? fr_zero() : polynomial.front();
            return;
        }

        const std::vector<Fr> left_remainder =
            polynomial_remainder(polynomial, tree.polynomial(node.left));
        const std::vector<Fr> right_remainder =
            polynomial_remainder(polynomial, tree.polynomial(node.right));
        evaluate_many_recursive(tree, node.left, left_remainder, out);
        evaluate_many_recursive(tree, node.right, right_remainder, out);
    }

    std::vector<Fr> evaluate_many_with_tree(
        const ProductTree& tree,
        const std::span<const Fr> polynomial)
    {
        const ProductTree::Node& root = tree.nodes[tree.root];
        const std::vector<Fr> remainder =
            polynomial_remainder(polynomial, tree.polynomial(tree.root));

        std::vector<Fr> values(root.end - root.begin, fr_zero());
        evaluate_many_recursive(tree, tree.root, remainder, values);
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

    std::span<Fr> reconstruction_span(
        const ProductTree& tree,
        const std::size_t node_index,
        const std::span<Fr> reconstruction) noexcept
    {
        const ProductTree::Node& node = tree.nodes[node_index];
        return reconstruction.subspan(
            node.polynomial_offset,
            node.end - node.begin);
    }

    void interpolate_recursive_into(
        const ProductTree& tree,
        const std::size_t node_index,
        const std::span<const Fr> weights,
        const std::span<Fr> reconstruction,
        const std::span<Fr> temporary)
    {
        const ProductTree::Node& node = tree.nodes[node_index];
        const std::span<Fr> out =
            reconstruction_span(tree, node_index, reconstruction);
        if (node.end - node.begin == 1) {
            out.front() = weights[node.begin];
            return;
        }

        interpolate_recursive_into(
            tree,
            node.left,
            weights,
            reconstruction,
            temporary);
        interpolate_recursive_into(
            tree,
            node.right,
            weights,
            reconstruction,
            temporary);

        const std::span<const Fr> left =
            reconstruction_span(tree, node.left, reconstruction);
        const std::span<const Fr> right =
            reconstruction_span(tree, node.right, reconstruction);
        multiply_polynomials_into(
            left,
            tree.polynomial(node.right),
            out);
        const std::span<Fr> other = temporary.first(out.size());
        multiply_polynomials_into(
            right,
            tree.polynomial(node.left),
            other);
        for (std::size_t i = 0; i < out.size(); ++i) {
            out[i] = fr_add(out[i], other[i]);
        }
    }

    void interpolate_polynomial(
        const std::span<const Fr> points,
        const std::span<const Fr> values,
        const std::span<Fr> out)
    {
        const ProductTree tree = build_product_tree(points);
        const std::vector<Fr> derivative =
            derivative_polynomial(tree.polynomial(tree.root));
        const std::vector<Fr> denominators =
            evaluate_many_with_tree(tree, derivative);
        const std::vector<Fr> denominator_inverses =
            batch_invert_nonzero(denominators);

        std::vector<Fr> weights(points.size(), fr_zero());
        for (std::size_t i = 0; i < points.size(); ++i) {
            weights[i] = fr_mul(values[i], denominator_inverses[i]);
        }

        // Reuse the public tree layout as offsets into one reconstruction arena.
        // Every node retains its full mathematical width: no value-dependent
        // trimming, allocation size, or naive/FFT dispatch is possible here.
        std::vector<Fr> reconstruction(tree.coefficients.size(), fr_zero());
        std::vector<Fr> temporary(points.size(), fr_zero());
        interpolate_recursive_into(
            tree,
            tree.root,
            weights,
            reconstruction,
            temporary);
        const std::span<Fr> coefficients =
            reconstruction_span(tree, tree.root, reconstruction);
        BLST_ASSERT(coefficients.size() == out.size());
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
    constexpr std::size_t max_points =
        std::bit_floor(static_cast<std::size_t>(std::numeric_limits<int>::max())) - 1;
    if (points.size() > max_points) {
        return Status::invalid_size;
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
        const std::span<Fr> product(cursor, n + 1);
        cursor += n + 1;
        const std::span<Fr> domain_work(cursor, fft_len);
        cursor += fft_len;
        const std::span<Fr> denominators(cursor, n);
        cursor += n;
        const std::span<Fr> prefixes(cursor, n + 1);

        fft_detail::fill_roots<Fr, FrCoeffCtx>(roots, false);

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
        const Fr inverse_fft_len = fr_inverse(fft_len_scalar);
        fft_with_forward_roots<Fr, Fr, FrCoeffCtx>(
            domain_work,
            roots,
            true,
            inverse_fft_len);

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
