#pragma once

#include <span>
#include <vector>

#include "blst/curve/GT.h"
#include "blst/curve/point.h"
#include "blst/fields/Fr.h"

namespace bte::building_blocks
{
    blst::G1Point commit_polynomial(
        std::span<const blst::G1Point> powers,
        const std::vector<blst::Fr>& coefficients);

    std::vector<blst::G1Point> compute_fk_eval_proofs(
        std::span<const blst::G1Point> powers,
        std::span<const blst::Fr> roots);

    bool verify_zero_opening(
        const blst::G1Point& commitment,
        const blst::G1Point& proof,
        const blst::Fr& root,
        const blst::G2Point& h_tau);
}
