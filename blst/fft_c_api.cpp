#include "bindings/blst_fft.h"

#include "blst/fft/context.h"
#include "blst/fft/ctFft.h"
#include "blst/fft/primitveRoot.h"

#include <algorithm>
#include <memory>
#include <span>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace {

using NativeFr = blst::Fr;

struct BlstFrLayout
{
    blst::limb_t l[4];
};

static_assert(sizeof(BlstFrLayout) == sizeof(NativeFr));
static_assert(alignof(BlstFrLayout) == alignof(NativeFr));
static_assert(sizeof(blst_fr) == sizeof(NativeFr));
static_assert(alignof(blst_fr) == alignof(NativeFr));
static_assert(std::is_standard_layout_v<NativeFr>);
static_assert(std::is_trivially_copyable_v<NativeFr>);

bool valid_fft_size(const std::size_t n)
{
    if (n == 0 || (n & (n - 1)) != 0) {
        return false;
    }
    return n <= static_cast<std::size_t>(std::numeric_limits<int>::max());
}

std::vector<NativeFr> make_roots(const std::size_t n, const bool inverse)
{
    NativeFr primitive = blst::primRootOfUnity<NativeFr>(n);
    if (inverse) {
        NativeFr::reciprocal(primitive, primitive);
    }
    std::vector<NativeFr> roots(n);
    NativeFr current = NativeFr::one();

    for (std::size_t i = 0; i < n; ++i) {
        roots[i] = current;
        NativeFr::mul(current, current, primitive);
    }

    return roots;
}

NativeFr fr_zero()
{
    return NativeFr{};
}

bool fr_is_zero(const NativeFr& value)
{
    return NativeFr::equal(value, fr_zero());
}

NativeFr fr_add(const NativeFr& lhs, const NativeFr& rhs)
{
    NativeFr out{};
    NativeFr::add(out, lhs, rhs);
    return out;
}

NativeFr fr_sub(const NativeFr& lhs, const NativeFr& rhs)
{
    NativeFr out{};
    NativeFr::sub(out, lhs, rhs);
    return out;
}

NativeFr fr_mul(const NativeFr& lhs, const NativeFr& rhs)
{
    NativeFr out{};
    NativeFr::mul(out, lhs, rhs);
    return out;
}

NativeFr fr_neg(const NativeFr& value)
{
    NativeFr out{};
    NativeFr::neg(out, value);
    return out;
}

NativeFr fr_inverse(const NativeFr& value)
{
    NativeFr out{};
    NativeFr::reciprocal(out, value);
    return out;
}

std::size_t next_power_of_two(const std::size_t value)
{
    std::size_t out = 1;
    while (out < value) {
        out <<= 1;
    }
    return out;
}

std::vector<NativeFr> trim_polynomial(std::vector<NativeFr> polynomial)
{
    while (!polynomial.empty() && fr_is_zero(polynomial.back())) {
        polynomial.pop_back();
    }
    return polynomial;
}

std::vector<NativeFr> take_prefix(
    const std::span<const NativeFr> polynomial,
    const std::size_t count)
{
    const std::size_t actual = std::min(polynomial.size(), count);
    return std::vector<NativeFr>(polynomial.begin(), polynomial.begin() + actual);
}

std::vector<NativeFr> reverse_polynomial(
    const std::span<const NativeFr> polynomial,
    const std::size_t size = 0)
{
    const std::size_t out_size = size == 0 ? polynomial.size() : size;
    std::vector<NativeFr> reversed(out_size, fr_zero());
    const std::size_t copy = std::min(polynomial.size(), out_size);
    for (std::size_t i = 0; i < copy; ++i) {
        reversed[i] = polynomial[polynomial.size() - 1 - i];
    }
    return reversed;
}

std::vector<NativeFr> add_polynomials(
    const std::span<const NativeFr> lhs,
    const std::span<const NativeFr> rhs)
{
    std::vector<NativeFr> out(std::max(lhs.size(), rhs.size()), fr_zero());
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        out[i] = fr_add(out[i], lhs[i]);
    }
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        out[i] = fr_add(out[i], rhs[i]);
    }
    return trim_polynomial(std::move(out));
}

std::vector<NativeFr> subtract_polynomials(
    const std::span<const NativeFr> lhs,
    const std::span<const NativeFr> rhs)
{
    std::vector<NativeFr> out(std::max(lhs.size(), rhs.size()), fr_zero());
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        out[i] = fr_add(out[i], lhs[i]);
    }
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        out[i] = fr_sub(out[i], rhs[i]);
    }
    return trim_polynomial(std::move(out));
}

std::vector<NativeFr> naive_multiply(
    const std::span<const NativeFr> lhs,
    const std::span<const NativeFr> rhs)
{
    if (lhs.empty() || rhs.empty()) {
        return {};
    }

    std::vector<NativeFr> out(lhs.size() + rhs.size() - 1, fr_zero());
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        for (std::size_t j = 0; j < rhs.size(); ++j) {
            out[i + j] = fr_add(out[i + j], fr_mul(lhs[i], rhs[j]));
        }
    }
    return trim_polynomial(std::move(out));
}

std::vector<NativeFr> multiply_polynomials(
    const std::span<const NativeFr> lhs,
    const std::span<const NativeFr> rhs)
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
    std::vector<NativeFr> a(fft_size, fr_zero());
    std::vector<NativeFr> b(fft_size, fr_zero());
    std::vector<NativeFr> c(fft_size, fr_zero());

    std::copy(lhs.begin(), lhs.end(), a.begin());
    std::copy(rhs.begin(), rhs.end(), b.begin());

    const std::vector<NativeFr> roots = make_roots(fft_size, false);
    const std::vector<NativeFr> inv_roots = make_roots(fft_size, true);
    poly_mult_fft<NativeFr, NativeFr, NativeFr, NativeFr, blst::FrCoeffCtx>(
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

std::vector<NativeFr> derivative_polynomial(const std::span<const NativeFr> polynomial)
{
    if (polynomial.size() <= 1) {
        return {};
    }

    std::vector<NativeFr> derivative(polynomial.size() - 1, fr_zero());
    for (std::size_t i = 1; i < polynomial.size(); ++i) {
        derivative[i - 1] = fr_mul(polynomial[i], NativeFr::from_uint64(i));
    }
    return trim_polynomial(std::move(derivative));
}

std::vector<NativeFr> inverse_series(
    const std::span<const NativeFr> polynomial,
    const std::size_t size)
{
    if (size == 0) {
        return {};
    }
    if (polynomial.empty() || fr_is_zero(polynomial.front())) {
        throw std::invalid_argument("power series inverse requires nonzero constant");
    }

    std::vector<NativeFr> inverse(1, fr_inverse(polynomial.front()));
    std::size_t width = 1;
    while (width < size) {
        const std::size_t next_width = std::min(width << 1, size);
        const std::vector<NativeFr> truncated = take_prefix(polynomial, next_width);
        std::vector<NativeFr> product = multiply_polynomials(truncated, inverse);
        product.resize(next_width, fr_zero());
        product[0] = fr_sub(NativeFr::from_uint64(2), product[0]);
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

std::pair<std::vector<NativeFr>, std::vector<NativeFr>> divide_remainder(
    const std::span<const NativeFr> dividend,
    const std::span<const NativeFr> divisor)
{
    if (divisor.empty()) {
        throw std::invalid_argument("polynomial division requires divisor");
    }

    if (dividend.size() < divisor.size()) {
        return { {}, trim_polynomial(std::vector<NativeFr>(dividend.begin(), dividend.end())) };
    }

    constexpr std::size_t naive_threshold = 32;
    if (divisor.size() <= naive_threshold) {
        std::vector<NativeFr> remainder(dividend.begin(), dividend.end());
        std::vector<NativeFr> quotient(dividend.size() - divisor.size() + 1, fr_zero());
        const NativeFr divisor_lead_inverse = fr_inverse(divisor.back());

        for (std::size_t idx = quotient.size(); idx-- > 0;) {
            const std::size_t degree = idx + divisor.size() - 1;
            if (fr_is_zero(remainder[degree])) {
                continue;
            }

            quotient[idx] = fr_mul(remainder[degree], divisor_lead_inverse);
            for (std::size_t j = 0; j < divisor.size(); ++j) {
                remainder[idx + j] = fr_sub(remainder[idx + j], fr_mul(quotient[idx], divisor[j]));
            }
        }

        remainder.resize(divisor.size() - 1);
        return { trim_polynomial(std::move(quotient)), trim_polynomial(std::move(remainder)) };
    }

    const std::size_t quotient_size = dividend.size() - divisor.size() + 1;
    const std::vector<NativeFr> reversed_dividend = reverse_polynomial(dividend);
    const std::vector<NativeFr> reversed_divisor = reverse_polynomial(divisor);
    const std::vector<NativeFr> reversed_prefix = take_prefix(reversed_dividend, quotient_size);
    const std::vector<NativeFr> inverse_reversed_divisor =
        inverse_series(reversed_divisor, quotient_size);

    std::vector<NativeFr> quotient_reversed =
        multiply_polynomials(reversed_prefix, inverse_reversed_divisor);
    quotient_reversed.resize(quotient_size, fr_zero());

    std::vector<NativeFr> quotient = trim_polynomial(reverse_polynomial(quotient_reversed));
    const std::vector<NativeFr> product = multiply_polynomials(quotient, divisor);
    std::vector<NativeFr> remainder = subtract_polynomials(dividend, product);
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
    std::vector<NativeFr> polynomial;
    std::unique_ptr<ProductTreeNode> left;
    std::unique_ptr<ProductTreeNode> right;
};

std::unique_ptr<ProductTreeNode> build_product_tree(
    const std::span<const NativeFr> points,
    const std::size_t begin,
    const std::size_t end)
{
    auto node = std::make_unique<ProductTreeNode>();
    node->begin = begin;
    node->end = end;

    if (end - begin == 1) {
        node->polynomial = { fr_neg(points[begin]), NativeFr::one() };
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
    const std::span<const NativeFr> polynomial,
    std::vector<NativeFr>& out)
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

std::vector<NativeFr> evaluate_many_with_tree(
    const ProductTreeNode& tree,
    const std::span<const NativeFr> polynomial)
{
    const auto [ignored_quotient, remainder] = divide_remainder(polynomial, tree.polynomial);
    (void)ignored_quotient;

    std::vector<NativeFr> values(tree.end - tree.begin, fr_zero());
    evaluate_many_recursive(tree, remainder, values);
    return values;
}

std::vector<NativeFr> batch_invert_nonzero(const std::span<const NativeFr> values)
{
    std::vector<NativeFr> prefixes(values.size() + 1, NativeFr::one());
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (fr_is_zero(values[i])) {
            throw std::invalid_argument("interpolation points must be distinct");
        }
        prefixes[i + 1] = fr_mul(prefixes[i], values[i]);
    }

    NativeFr all_inverse = fr_inverse(prefixes.back());
    std::vector<NativeFr> inverses(values.size(), NativeFr::one());
    for (std::size_t i = values.size(); i-- > 0;) {
        inverses[i] = fr_mul(all_inverse, prefixes[i]);
        all_inverse = fr_mul(all_inverse, values[i]);
    }
    return inverses;
}

void batch_invert_nonzero_into(
    const std::span<const NativeFr> values,
    const std::span<NativeFr> prefixes,
    const std::span<NativeFr> inverses)
{
    if (prefixes.size() != values.size() + 1 || inverses.size() != values.size()) {
        throw std::invalid_argument("batch inversion scratch shape mismatch");
    }

    prefixes[0] = NativeFr::one();
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (fr_is_zero(values[i])) {
            throw std::invalid_argument("interpolation points must be distinct");
        }
        prefixes[i + 1] = fr_mul(prefixes[i], values[i]);
    }

    NativeFr all_inverse = fr_inverse(prefixes.back());
    for (std::size_t i = values.size(); i-- > 0;) {
        inverses[i] = fr_mul(all_inverse, prefixes[i]);
        all_inverse = fr_mul(all_inverse, values[i]);
    }
}

std::vector<NativeFr> interpolate_recursive(
    const ProductTreeNode& node,
    const std::span<const NativeFr> weights)
{
    if (node.end - node.begin == 1) {
        return { weights[node.begin] };
    }

    std::vector<NativeFr> left = interpolate_recursive(*node.left, weights);
    std::vector<NativeFr> right = interpolate_recursive(*node.right, weights);
    left = multiply_polynomials(left, node.right->polynomial);
    right = multiply_polynomials(right, node.left->polynomial);
    return add_polynomials(left, right);
}

void interpolate_polynomial(
    const std::span<const NativeFr> points,
    const std::span<const NativeFr> values,
    const std::span<NativeFr> out)
{
    if (points.empty() || points.size() != values.size() || points.size() != out.size()) {
        throw std::invalid_argument("interpolation input shape mismatch");
    }

    const std::unique_ptr<ProductTreeNode> tree = build_product_tree(points, 0, points.size());
    const std::vector<NativeFr> derivative = derivative_polynomial(tree->polynomial);
    const std::vector<NativeFr> denominators = evaluate_many_with_tree(*tree, derivative);
    const std::vector<NativeFr> denominator_inverses = batch_invert_nonzero(denominators);

    std::vector<NativeFr> weights(points.size(), fr_zero());
    for (std::size_t i = 0; i < points.size(); ++i) {
        weights[i] = fr_mul(values[i], denominator_inverses[i]);
    }

    std::vector<NativeFr> coefficients = interpolate_recursive(*tree, weights);
    coefficients.resize(points.size(), fr_zero());
    std::copy(coefficients.begin(), coefficients.end(), out.begin());
}

std::size_t interpolate_fft_domain_suffix_scratch_len(
    const std::size_t n,
    const std::size_t fft_len)
{
    if (n == 0 || !valid_fft_size(fft_len)) {
        return 0;
    }

    return 3 * fft_len + 3 * n + 2;
}

bool interpolate_fft_domain_suffix_with_scratch(
    NativeFr* const out_suffix,
    const std::size_t* const indices,
    const NativeFr* const values,
    const std::size_t n,
    const std::size_t ell,
    const std::size_t fft_len,
    NativeFr* const scratch,
    const std::size_t scratch_len)
{
    if (out_suffix == nullptr || indices == nullptr || values == nullptr ||
        scratch == nullptr || n == 0 || ell > n || !valid_fft_size(fft_len)) {
        return false;
    }

    const std::size_t required_scratch =
        interpolate_fft_domain_suffix_scratch_len(n, fft_len);
    if (scratch_len < required_scratch) {
        return false;
    }

    try {
        NativeFr* cursor = scratch;
        const std::span<NativeFr> roots(cursor, fft_len);
        cursor += fft_len;
        const std::span<NativeFr> inv_roots(cursor, fft_len);
        cursor += fft_len;
        const std::span<NativeFr> product(cursor, n + 1);
        cursor += n + 1;
        const std::span<NativeFr> domain_work(cursor, fft_len);
        cursor += fft_len;
        const std::span<NativeFr> denominators(cursor, n);
        cursor += n;
        const std::span<NativeFr> prefixes(cursor, n + 1);
        (void)required_scratch;

        const NativeFr primitive = blst::primRootOfUnity<NativeFr>(fft_len);
        NativeFr inverse_primitive = primitive;
        NativeFr::reciprocal(inverse_primitive, inverse_primitive);
        NativeFr current = NativeFr::one();
        NativeFr inverse_current = NativeFr::one();
        for (std::size_t i = 0; i < fft_len; ++i) {
            roots[i] = current;
            inv_roots[i] = inverse_current;
            NativeFr::mul(current, current, primitive);
            NativeFr::mul(inverse_current, inverse_current, inverse_primitive);
        }

        std::fill(product.begin(), product.end(), fr_zero());
        product[0] = NativeFr::one();
        std::size_t degree = 0;
        for (std::size_t sample = 0; sample < n; ++sample) {
            if (indices[sample] >= fft_len) {
                return false;
            }
            const NativeFr& point = roots[indices[sample]];
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
                fr_mul(product[coeff], NativeFr::from_uint64(coeff));
        }
        fft<NativeFr, NativeFr, blst::FrCoeffCtx>(
            domain_work,
            std::span<const NativeFr>(roots),
            false);

        for (std::size_t sample = 0; sample < n; ++sample) {
            const NativeFr denominator = domain_work[indices[sample]];
            if (fr_is_zero(denominator)) {
                return false;
            }
            denominators[sample] = denominator;
        }
        batch_invert_nonzero_into(denominators, prefixes, prefixes.subspan(0, n));

        std::fill(domain_work.begin(), domain_work.end(), fr_zero());
        const NativeFr fft_len_scalar = NativeFr::from_uint64(fft_len);
        for (std::size_t sample = 0; sample < n; ++sample) {
            const std::size_t index = indices[sample];
            const NativeFr& inverse_root = roots[(fft_len - index) % fft_len];
            const NativeFr domain_derivative = fr_mul(fft_len_scalar, inverse_root);
            const NativeFr weighted_value = fr_mul(values[sample], prefixes[sample]);
            domain_work[index] = fr_mul(weighted_value, domain_derivative);
        }
        fft<NativeFr, NativeFr, blst::FrCoeffCtx>(
            domain_work,
            std::span<const NativeFr>(inv_roots),
            true);

        const std::size_t suffix_len = n - ell;
        for (std::size_t suffix_index = 0; suffix_index < suffix_len; ++suffix_index) {
            const std::size_t target_degree = fft_len + ell + suffix_index;
            const std::size_t min_product_index = ell + suffix_index + 1;
            NativeFr acc = fr_zero();
            for (std::size_t product_index = min_product_index; product_index <= n;
                 ++product_index) {
                const std::size_t quotient_index = target_degree - product_index;
                acc = fr_add(acc, fr_mul(product[product_index], domain_work[quotient_index]));
            }
            out_suffix[suffix_index] = acc;
        }

        return true;
    } catch (...) {
        return false;
    }
}

} // namespace

extern "C" {

bool blst_fr_fft(blst_fr* values, size_t n, bool inverse)
{
    if (values == nullptr || !valid_fft_size(n)) {
        return false;
    }

    try {
        std::vector<NativeFr> roots = make_roots(n, inverse);
        auto* native_values = reinterpret_cast<NativeFr*>(values);
        fft<NativeFr, NativeFr, blst::FrCoeffCtx>(
            std::span<NativeFr>(native_values, n),
            std::span<const NativeFr>(roots),
            inverse);
        return true;
    } catch (...) {
        return false;
    }
}

bool blst_fr_fft_transpose_forward(blst_fr* values, size_t n)
{
    // The exported Fr FFT is the natural-order DFT over the subgroup:
    // F[i,j] = omega^(i*j). That matrix is symmetric, so the full forward
    // transform is also its transpose. Keep this in C++/BLST so Rust protocol
    // code dispatches to the owned FFT implementation instead of spelling out
    // a Vandermonde contraction.
    return blst_fr_fft(values, n, false);
}

bool blst_fr_interpolate(
    blst_fr* out_coefficients,
    const blst_fr* points,
    const blst_fr* values,
    size_t n)
{
    if (out_coefficients == nullptr || points == nullptr || values == nullptr || n == 0) {
        return false;
    }

    try {
        auto* native_out = reinterpret_cast<NativeFr*>(out_coefficients);
        const auto* native_points = reinterpret_cast<const NativeFr*>(points);
        const auto* native_values = reinterpret_cast<const NativeFr*>(values);
        interpolate_polynomial(
            std::span<const NativeFr>(native_points, n),
            std::span<const NativeFr>(native_values, n),
            std::span<NativeFr>(native_out, n));
        return true;
    } catch (...) {
        return false;
    }
}

bool blst_fr_interpolate_fft_domain(
    blst_fr* out_coefficients,
    const size_t* indices,
    const blst_fr* values,
    size_t n,
    size_t fft_len)
{
    if (out_coefficients == nullptr || indices == nullptr || values == nullptr ||
        n == 0 || !valid_fft_size(fft_len)) {
        return false;
    }

    try {
        const std::vector<NativeFr> roots = make_roots(fft_len, false);
        std::vector<NativeFr> points;
        points.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            if (indices[i] >= fft_len) {
                return false;
            }
            points.push_back(roots[indices[i]]);
        }

        auto* native_out = reinterpret_cast<NativeFr*>(out_coefficients);
        const auto* native_values = reinterpret_cast<const NativeFr*>(values);
        interpolate_polynomial(
            std::span<const NativeFr>(points.data(), points.size()),
            std::span<const NativeFr>(native_values, n),
            std::span<NativeFr>(native_out, n));
        return true;
    } catch (...) {
        return false;
    }
}

bool blst_fr_interpolate_fft_domain_suffix(
    blst_fr* out_suffix,
    const size_t* indices,
    const blst_fr* values,
    size_t n,
    size_t ell,
    size_t fft_len)
{
    const std::size_t scratch_len =
        interpolate_fft_domain_suffix_scratch_len(n, fft_len);
    if (scratch_len == 0) {
        return false;
    }

    try {
        std::vector<NativeFr> scratch(scratch_len);
        return interpolate_fft_domain_suffix_with_scratch(
            reinterpret_cast<NativeFr*>(out_suffix),
            indices,
            reinterpret_cast<const NativeFr*>(values),
            n,
            ell,
            fft_len,
            scratch.data(),
            scratch.size());
    } catch (...) {
        return false;
    }
}

size_t blst_fr_interpolate_fft_domain_suffix_scratch_len(size_t n, size_t fft_len)
{
    return interpolate_fft_domain_suffix_scratch_len(n, fft_len);
}

bool blst_fr_interpolate_fft_domain_suffix_with_scratch(
    blst_fr* out_suffix,
    const size_t* indices,
    const blst_fr* values,
    size_t n,
    size_t ell,
    size_t fft_len,
    blst_fr* scratch,
    size_t scratch_len)
{
    return interpolate_fft_domain_suffix_with_scratch(
        reinterpret_cast<NativeFr*>(out_suffix),
        indices,
        reinterpret_cast<const NativeFr*>(values),
        n,
        ell,
        fft_len,
        reinterpret_cast<NativeFr*>(scratch),
        scratch_len);
}

}
