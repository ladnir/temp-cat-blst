// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#include "fuzz/internal/adapter.h"

#include <array>
#include <cstddef>
#include <span>
#include <vector>

int main()
{
    using namespace catblst::fuzz;
    namespace adapter = catblst::fuzz::adapter;

    adapter::run_abi_regressions();

    const std::array<std::byte, 1> one{std::byte{1}};
    const std::span<const std::byte> empty;

    const ReducedFpValue fp_zero = adapter::blst_fp_from_le_bytes(empty);
    const ReducedFpValue fp_one = adapter::blst_fp_from_le_bytes(one);
    if (adapter::refine_nonzero(fp_zero) || !adapter::refine_nonzero(fp_one))
        return 1;

    const ReducedFp2Value fp2_zero = adapter::reduce_fp2(empty);
    const ReducedFp2Value fp2_one = adapter::reduce_fp2(one);
    if (adapter::refine_nonzero(fp2_zero) || !adapter::refine_nonzero(fp2_one))
        return 2;

    const ReducedFp12Value fp12_zero = adapter::reduce_fp12(empty);
    const ReducedFp12Value fp12_one = adapter::reduce_fp12(one);
    if (adapter::refine_nonzero(fp12_zero) || !adapter::refine_nonzero(fp12_one))
        return 3;

    const SubgroupG1AffineValue g1 = adapter::blst_p1_affine_generator();
    const OnCurveG1AffineValue g1_curve = adapter::as_on_curve(g1);
    const auto g1_decoded = adapter::blst_p1_deserialize(
        adapter::blst_p1_affine_compress(g1_curve));
    if (!g1_decoded || !adapter::blst_p1_affine_in_g1(*g1_decoded))
        return 4;
    if (adapter::blst_p1_affine_in_g1(adapter::retained_non_subgroup_g1()))
        return 12;
    if (!adapter::blst_p1_affine_on_curve(adapter::as_candidate(g1_curve)))
        return 5;

    const SubgroupG2AffineValue g2 = adapter::blst_p2_affine_generator();
    const OnCurveG2AffineValue g2_curve = adapter::as_on_curve(g2);
    const auto g2_decoded = adapter::blst_p2_deserialize(
        adapter::blst_p2_affine_compress(g2_curve));
    if (!g2_decoded || !adapter::blst_p2_affine_in_g2(*g2_decoded))
        return 6;
    if (adapter::blst_p2_affine_in_g2(adapter::retained_non_subgroup_g2()))
        return 13;
    if (!adapter::blst_p2_affine_on_curve(adapter::as_candidate(g2_curve)))
        return 7;

    if (adapter::blst_p1_affine_on_curve(
            adapter::make_g1_candidate(fp_zero, fp_one)))
        return 8;
    if (adapter::blst_p2_affine_on_curve(
            adapter::make_g2_candidate(fp2_zero, fp2_one)))
        return 9;

    const OnCurveG1Value g1_projective = adapter::blst_p1_from_affine(g1_curve);
    const OnCurveG1Value g1_double = adapter::blst_p1_double(g1_projective);
    if (!adapter::incomplete_add_pair(g1_projective, g1_double) ||
        adapter::incomplete_add_pair(g1_projective, g1_projective) ||
        adapter::incomplete_add_pair(
            g1_projective,
            adapter::blst_p1_cneg(g1_projective)))
        return 10;

    const OnCurveG2Value g2_projective = adapter::blst_p2_from_affine(g2_curve);
    const OnCurveG2Value g2_double = adapter::blst_p2_double(g2_projective);
    if (!adapter::incomplete_add_pair(g2_projective, g2_double) ||
        adapter::incomplete_add_pair(g2_projective, g2_projective) ||
        adapter::incomplete_add_pair(
            g2_projective,
            adapter::blst_p2_cneg(g2_projective)))
        return 11;

    const MillerValueStorage miller = adapter::blst_miller_loop(g2, g1);
    [[maybe_unused]] const GtElementValue gt = adapter::blst_final_exp(miller);

    if (adapter::prepare_fft_domain(0) || adapter::prepare_fft_domain(3) ||
        !adapter::prepare_fft_domain(4))
        return 14;

    const SubgroupG2Value g2_subgroup = adapter::blst_p2_from_affine(g2);
    const std::array<SubgroupG2Value, 2> b_values{g2_subgroup, g2_subgroup};
    if (adapter::precompute_middle_product_b(b_values, 0) ||
        adapter::precompute_middle_product_b(b_values, 2) ||
        !adapter::precompute_middle_product_b(b_values, 1))
        return 15;

    if (!adapter::make_pairing_state(false, empty))
        return 16;
    auto pending = adapter::make_pairing_state(true, one);
    if (!pending)
        return 17;
    [[maybe_unused]] PairingUnsetCommitted committed =
        adapter::pairing_commit(std::move(*pending));
    return 0;
}
