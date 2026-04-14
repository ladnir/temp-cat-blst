#pragma once

#include <span>
#include <utility>
#include <vector>

#include "blst/fields/Fr.h"

namespace bte::building_blocks
{
    std::vector<blst::Fr> multiply_polynomials(
        std::span<const blst::Fr> lhs,
        std::span<const blst::Fr> rhs);

    std::vector<blst::Fr> derivative_polynomial(
        std::span<const blst::Fr> polynomial);

    std::pair<std::vector<blst::Fr>, std::vector<blst::Fr>> divide_remainder(
        std::span<const blst::Fr> dividend,
        std::span<const blst::Fr> divisor);

    std::vector<blst::Fr> evaluate_many(
        std::span<const blst::Fr> polynomial,
        std::span<const blst::Fr> points);

    std::vector<std::vector<blst::Fr>> exclusion_quotients(
        std::span<const blst::Fr> roots);

    std::vector<blst::Fr> build_vanishing_polynomial(const std::vector<blst::Fr>& roots);

    std::vector<blst::Fr> synthetic_divide(
        const std::vector<blst::Fr>& polynomial,
        const blst::Fr& root,
        blst::Fr* remainder = nullptr);
}
