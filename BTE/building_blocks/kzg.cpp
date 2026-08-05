#include "BTE/building_blocks/kzg.h"

#include <array>
#include <cstring>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>

#include "blst/fft/ctFft.h"
#include "blst/fft/primitveRoot.h"

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
            const blst::Fr zero = fr_zero();
            return blst::Fr::equal(value, zero);
        }

        std::array<std::uint8_t, 32> fr_to_scalar_bytes(const blst::Fr& scalar)
        {
            std::array<std::uint8_t, 32> bytes{};
            blst::Fr::to_scaler(bytes.data(), scalar);
            return bytes;
        }

        blst::G1Point scale_g1(const blst::G1Point& point, const blst::Fr& scalar)
        {
            if (fr_is_zero(scalar)) {
                return blst::G1Point::infinity();
            }
            const std::array<std::uint8_t, 32> scalar_bytes = fr_to_scalar_bytes(scalar);
            blst::G1Point out{};
            blst::G1Point::mult(out, point, scalar_bytes.data(), 256);
            return out;
        }

        blst::G2Point scale_g2(const blst::G2Point& point, const blst::Fr& scalar)
        {
            if (fr_is_zero(scalar)) {
                return blst::G2Point::infinity();
            }
            const std::array<std::uint8_t, 32> scalar_bytes = fr_to_scalar_bytes(scalar);
            blst::G2Point out{};
            blst::G2Point::mult(out, point, scalar_bytes.data(), 256);
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

        blst::G2Point sub_g2(const blst::G2Point& lhs, const blst::G2Point& rhs)
        {
            blst::G2Point neg_rhs{};
            blst::G2Point::neg(neg_rhs, rhs);
            return add_g2(lhs, neg_rhs);
        }

        blst::G1Point sub_g1(const blst::G1Point& lhs, const blst::G1Point& rhs)
        {
            blst::G1Point neg_rhs{};
            blst::G1Point::neg(neg_rhs, rhs);
            return add_g1(lhs, neg_rhs);
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

        std::vector<blst::Fr> multiply_polynomials(
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
            return out;
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

        std::uint64_t pack_pair_key(const std::size_t a, const std::size_t b)
        {
            return (static_cast<std::uint64_t>(a) << 32U) | static_cast<std::uint64_t>(b);
        }

        struct TripleKey
        {
            std::size_t a = 0;
            std::size_t b = 0;
            std::size_t c = 0;

            bool operator==(const TripleKey& other) const = default;
        };

        struct TripleKeyHash
        {
            std::size_t operator()(const TripleKey& key) const noexcept
            {
                std::size_t seed = key.a;
                seed ^= key.b + 0x9e3779b9U + (seed << 6U) + (seed >> 2U);
                seed ^= key.c + 0x9e3779b9U + (seed << 6U) + (seed >> 2U);
                return seed;
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
            const std::span<const blst::Fr> roots,
            const std::size_t begin,
            const std::size_t end)
        {
            auto node = std::make_unique<ProductTreeNode>();
            node->begin = begin;
            node->end = end;
            if (end - begin == 1) {
                node->polynomial = { fr_neg(roots[begin]), blst::Fr::one() };
                return node;
            }

            const std::size_t mid = begin + ((end - begin) >> 1);
            node->left = build_product_tree(roots, begin, mid);
            node->right = build_product_tree(roots, mid, end);
            node->polynomial = multiply_polynomials(
                std::span<const blst::Fr>(node->left->polynomial.data(), node->left->polynomial.size()),
                std::span<const blst::Fr>(node->right->polynomial.data(), node->right->polynomial.size()));
            return node;
        }

        struct G1FrCoeffCtx
        {
            static blst::G1Point plus(const blst::G1Point& lhs, const blst::G1Point& rhs)
            {
                return add_g1(lhs, rhs);
            }

            static blst::G1Point minus(const blst::G1Point& lhs, const blst::G1Point& rhs)
            {
                return sub_g1(lhs, rhs);
            }

            static blst::G1Point mult(const blst::G1Point& point, const blst::Fr& scalar)
            {
                return scale_g1(point, scalar);
            }

            static blst::Fr plus(const blst::Fr& lhs, const blst::Fr& rhs)
            {
                return fr_add(lhs, rhs);
            }

            static blst::Fr minus(const blst::Fr& lhs, const blst::Fr& rhs)
            {
                return fr_sub(lhs, rhs);
            }

            static blst::Fr mult(const blst::Fr& lhs, const blst::Fr& rhs)
            {
                return fr_mul(lhs, rhs);
            }

            static blst::Fr inverse(const blst::Fr& value)
            {
                return fr_inverse(value);
            }

            template<typename F>
            static F from_u64(const std::uint64_t value)
            {
                static_assert(std::is_same_v<F, blst::Fr>);
                return blst::Fr::from_uint64(value);
            }
        };

        struct FkFftCache
        {
            std::unordered_map<std::size_t, std::vector<blst::Fr>> roots_by_size;
            std::unordered_map<std::size_t, std::vector<blst::Fr>> inv_roots_by_size;
            std::unordered_map<std::uint64_t, FftTruncatedForwardPlan> forward_plans;
            std::unordered_map<TripleKey, FftInverseTruncatedRangePlan, TripleKeyHash> inverse_plans;

            const std::vector<blst::Fr>& roots(const std::size_t n)
            {
                auto [it, inserted] = roots_by_size.try_emplace(n);
                if (inserted) {
                    it->second = make_fft_roots(n);
                }
                return it->second;
            }

            const std::vector<blst::Fr>& inv_roots(const std::size_t n)
            {
                auto [it, inserted] = inv_roots_by_size.try_emplace(n);
                if (inserted) {
                    it->second = make_ifft_roots(roots(n));
                }
                return it->second;
            }

            const FftTruncatedForwardPlan& forward_plan(const std::size_t n, const std::size_t xn)
            {
                const std::uint64_t key = pack_pair_key(n, xn);
                auto [it, inserted] = forward_plans.try_emplace(key);
                if (inserted) {
                    it->second = make_fft_truncated_forward_plan(n, xn);
                }
                return it->second;
            }

            const FftInverseTruncatedRangePlan& inverse_plan(
                const std::size_t n,
                const std::size_t start,
                const std::size_t length)
            {
                const TripleKey key{ n, start, length };
                auto [it, inserted] = inverse_plans.try_emplace(key);
                if (inserted) {
                    it->second = make_fft_inverse_truncated_range_plan(n, start, length);
                }
                return it->second;
            }
        };

        std::vector<blst::G1Point> middle_window_g1_fr(
            const std::span<const blst::G1Point> group_coefficients,
            const std::span<const blst::Fr> scalar_coefficients,
            const std::size_t window_start,
            const std::size_t window_length,
            FkFftCache& cache)
        {
            if (window_length == 0) {
                return {};
            }
            if (group_coefficients.empty() || scalar_coefficients.empty()) {
                throw std::invalid_argument("Middle product requires non-empty inputs");
            }

            const std::size_t convolution_size = group_coefficients.size() + scalar_coefficients.size() - 1;
            if (window_start + window_length > convolution_size) {
                throw std::invalid_argument("Requested middle window is out of range");
            }

            const std::size_t fft_size = next_power_of_two(convolution_size);
            std::vector<blst::G1Point> lhs(fft_size, blst::G1Point::infinity());
            std::copy(group_coefficients.begin(), group_coefficients.end(), lhs.begin());

            std::vector<blst::Fr> rhs(fft_size, fr_zero());
            for (std::size_t i = 0; i < scalar_coefficients.size(); ++i) {
                rhs[i] = scalar_coefficients[scalar_coefficients.size() - 1 - i];
            }

            const std::vector<blst::Fr>& roots = cache.roots(fft_size);
            const std::vector<blst::Fr>& inv_roots = cache.inv_roots(fft_size);
            const FftTruncatedForwardPlan& lhs_plan =
                cache.forward_plan(fft_size, group_coefficients.size());
            const FftTruncatedForwardPlan& rhs_plan =
                cache.forward_plan(fft_size, scalar_coefficients.size());

            fft_forward_truncated_prefix<blst::G1Point, blst::Fr, G1FrCoeffCtx>(
                std::span<blst::G1Point>(lhs.data(), lhs.size()),
                std::span<const blst::Fr>(roots.data(), roots.size()),
                group_coefficients.size(),
                lhs_plan);
            fft_forward_truncated_prefix<blst::Fr, blst::Fr, G1FrCoeffCtx>(
                std::span<blst::Fr>(rhs.data(), rhs.size()),
                std::span<const blst::Fr>(roots.data(), roots.size()),
                scalar_coefficients.size(),
                rhs_plan);

            std::vector<blst::G1Point> product(fft_size, blst::G1Point::infinity());
            for (std::size_t i = 0; i < fft_size; ++i) {
                product[i] = scale_g1(lhs[i], rhs[i]);
            }

            const FftInverseTruncatedRangePlan& inverse_plan =
                cache.inverse_plan(fft_size, window_start, window_length);
            fft_inverse_truncated_range<blst::G1Point, blst::Fr, G1FrCoeffCtx>(
                std::span<blst::G1Point>(product.data(), product.size()),
                std::span<const blst::Fr>(inv_roots.data(), inv_roots.size()),
                window_start,
                window_length,
                inverse_plan);

            std::vector<blst::G1Point> out(window_length, blst::G1Point::infinity());
            for (std::size_t i = 0; i < window_length; ++i) {
                out[i] = product[window_start + i];
            }
            return out;
        }

        std::vector<blst::G1Point> derive_child_crs(
            const std::span<const blst::G1Point> parent_crs,
            const std::span<const blst::Fr> sibling_polynomial,
            const std::size_t output_size,
            FkFftCache& cache)
        {
            if (output_size == 0) {
                return {};
            }
            if (sibling_polynomial.empty()) {
                throw std::invalid_argument("Sibling polynomial must be non-empty");
            }
            if (parent_crs.size() < output_size + sibling_polynomial.size() - 1) {
                throw std::invalid_argument("Parent CRS window is too short for FK derivation");
            }
            return middle_window_g1_fr(
                parent_crs,
                sibling_polynomial,
                sibling_polynomial.size() - 1,
                output_size,
                cache);
        }

        void compute_fk_eval_proofs_recursive(
            const ProductTreeNode& node,
            const std::span<const blst::G1Point> node_crs,
            std::vector<blst::G1Point>& out,
            FkFftCache& cache)
        {
            const std::size_t node_size = node.end - node.begin;
            if (node_crs.size() != node_size) {
                throw std::invalid_argument("FK node CRS length must match the node size");
            }

            if (node_size == 1) {
                out[node.begin] = node_crs.front();
                return;
            }

            const std::size_t left_size = node.left->end - node.left->begin;
            const std::size_t right_size = node.right->end - node.right->begin;

            const std::vector<blst::G1Point> left_crs = derive_child_crs(
                node_crs,
                std::span<const blst::Fr>(node.right->polynomial.data(), node.right->polynomial.size()),
                left_size,
                cache);
            const std::vector<blst::G1Point> right_crs = derive_child_crs(
                node_crs,
                std::span<const blst::Fr>(node.left->polynomial.data(), node.left->polynomial.size()),
                right_size,
                cache);

            compute_fk_eval_proofs_recursive(
                *node.left,
                std::span<const blst::G1Point>(left_crs.data(), left_crs.size()),
                out,
                cache);
            compute_fk_eval_proofs_recursive(
                *node.right,
                std::span<const blst::G1Point>(right_crs.data(), right_crs.size()),
                out,
                cache);
        }

        blst::GT pairing(const blst::G1Point& p, const blst::G2Point& q)
        {
            blst::G1AffinePoint p_aff{};
            blst::G2AffinePoint q_aff{};
            blst::G1Point::to_affine(p_aff, p);
            blst::G2Point::to_affine(q_aff, q);

            blst::MillerResult miller{};
            blst::GT::miller_loop(miller, &q_aff, &p_aff);
            blst::GT out{};
            blst::GT::final_exp(out, miller);
            return out;
        }
    }

    blst::G1Point commit_polynomial(
        const std::span<const blst::G1Point> powers,
        const std::vector<blst::Fr>& coefficients)
    {
        if (coefficients.empty()) {
            return blst::G1Point::infinity();
        }
        if (powers.size() < coefficients.size()) {
            throw std::invalid_argument("Not enough CRS powers for polynomial commitment");
        }

        blst::G1Point accumulator = blst::G1Point::infinity();
        for (std::size_t i = 0; i < coefficients.size(); ++i) {
            if (fr_is_zero(coefficients[i])) {
                continue;
            }
            accumulator = add_g1(accumulator, scale_g1(powers[i], coefficients[i]));
        }
        return accumulator;
    }

    std::vector<blst::G1Point> compute_fk_eval_proofs(
        const std::span<const blst::G1Point> powers,
        const std::span<const blst::Fr> roots)
    {
        if (roots.empty()) {
            return {};
        }
        if (powers.size() < roots.size()) {
            throw std::invalid_argument("Need at least batch-many CRS powers for FK evaluation proofs");
        }

        std::vector<blst::G1Point> out(roots.size(), blst::G1Point::infinity());
        const std::unique_ptr<ProductTreeNode> product_tree = build_product_tree(roots, 0, roots.size());
        FkFftCache cache{};
        compute_fk_eval_proofs_recursive(
            *product_tree,
            powers.first(roots.size()),
            out,
            cache);
        return out;
    }

    bool verify_zero_opening(
        const blst::G1Point& commitment,
        const blst::G1Point& proof,
        const blst::Fr& root,
        const blst::G2Point& h_tau)
    {
        const blst::G2Point scaled_h = scale_g2(blst::G2Point::generator(), root);
        const blst::G2Point opening_base = sub_g2(h_tau, scaled_h);

        const blst::GT lhs = pairing(proof, opening_base);
        const blst::GT rhs = pairing(commitment, blst::G2Point::generator());
        return blst::GT::equal(lhs, rhs) != 0;
    }
}
