#include "BTE/building_blocks/lagrange.h"

#include <cstring>
#include <stdexcept>
#include <unordered_set>

#include "BTE/building_blocks/polynomial.h"

namespace bte::building_blocks
{
    namespace
    {
        blst::Fr fr_from_u64(const std::uint64_t value)
        {
            return blst::Fr::from_uint64(value);
        }

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
    }

    std::vector<blst::Fr> lagrange_at_zero(const std::span<const std::uint32_t> indices)
    {
        if (indices.empty()) {
            throw std::invalid_argument("Need at least one index for Lagrange interpolation");
        }

        std::unordered_set<std::uint32_t> dedupe;
        dedupe.reserve(indices.size());
        for (const std::uint32_t index : indices) {
            if (index == 0) {
                throw std::invalid_argument("Validator indices must be nonzero");
            }
            if (!dedupe.insert(index).second) {
                throw std::invalid_argument("Duplicate validator index");
            }
        }

        std::vector<blst::Fr> points;
        points.reserve(indices.size());
        for (const std::uint32_t index : indices) {
            points.push_back(fr_from_u64(index));
        }

        const std::vector<blst::Fr> vanishing = build_vanishing_polynomial(points);
        const std::vector<blst::Fr> derivative = derivative_polynomial(vanishing);
        const std::vector<blst::Fr> derivative_values =
            evaluate_many(
                std::span<const blst::Fr>(derivative.data(), derivative.size()),
                std::span<const blst::Fr>(points.data(), points.size()));

        std::vector<blst::Fr> denominators(indices.size(), blst::Fr::one());
        for (std::size_t i = 0; i < indices.size(); ++i) {
            denominators[i] = fr_mul(points[i], derivative_values[i]);
            if (fr_is_zero(denominators[i])) {
                throw std::runtime_error("Lagrange denominator vanished unexpectedly");
            }
        }

        std::vector<blst::Fr> prefixes(indices.size() + 1, blst::Fr::one());
        for (std::size_t i = 0; i < indices.size(); ++i) {
            prefixes[i + 1] = fr_mul(prefixes[i], denominators[i]);
        }

        blst::Fr all_inverse{};
        blst::Fr::reciprocal(all_inverse, prefixes.back());

        std::vector<blst::Fr> denominator_inverses(indices.size(), blst::Fr::one());
        for (std::size_t i = indices.size(); i-- > 0;) {
            denominator_inverses[i] = fr_mul(all_inverse, prefixes[i]);
            all_inverse = fr_mul(all_inverse, denominators[i]);
        }

        const blst::Fr minus_m0 = fr_neg(vanishing.front());
        std::vector<blst::Fr> lambdas(indices.size(), blst::Fr::one());
        for (std::size_t i = 0; i < indices.size(); ++i) {
            lambdas[i] = fr_mul(minus_m0, denominator_inverses[i]);
        }

        return lambdas;
    }
}
