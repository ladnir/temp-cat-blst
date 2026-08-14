// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "fuzz/support/storage.h"

#include <optional>
#include <array>
#include <span>
#include <vector>

namespace catblst::fuzz::adapter
{
    [[nodiscard]] ReducedFrValue blst_scalar_from_le_bytes(
        std::span<const std::byte> input) noexcept;
    [[nodiscard]] ReducedFpValue blst_fp_from_le_bytes(
        std::span<const std::byte> input) noexcept;
    [[nodiscard]] ReducedFp2Value reduce_fp2(
        std::span<const std::byte> input) noexcept;
    [[nodiscard]] ReducedFp12Value reduce_fp12(
        std::span<const std::byte> input) noexcept;

    [[nodiscard]] std::optional<NonzeroFpValue> refine_nonzero(
        const ReducedFpValue& input) noexcept;
    [[nodiscard]] std::optional<NonzeroFp2Value> refine_nonzero(
        const ReducedFp2Value& input) noexcept;
    [[nodiscard]] std::optional<NonzeroFp12Value> refine_nonzero(
        const ReducedFp12Value& input) noexcept;

    [[nodiscard]] FieldValidG1CandidateValue make_g1_candidate(
        const ReducedFpValue& x,
        const ReducedFpValue& y) noexcept;
    [[nodiscard]] FieldValidG2CandidateValue make_g2_candidate(
        const ReducedFp2Value& x,
        const ReducedFp2Value& y) noexcept;

    [[nodiscard]] std::optional<OnCurveG1AffineValue> blst_p1_deserialize(
        std::span<const std::byte> input) noexcept;
    [[nodiscard]] std::optional<OnCurveG2AffineValue> blst_p2_deserialize(
        std::span<const std::byte> input) noexcept;
    [[nodiscard]] std::optional<OnCurveG1AffineValue> blst_p1_affine_on_curve(
        const FieldValidG1CandidateValue& input) noexcept;
    [[nodiscard]] std::optional<OnCurveG2AffineValue> blst_p2_affine_on_curve(
        const FieldValidG2CandidateValue& input) noexcept;
    [[nodiscard]] std::optional<SubgroupG1AffineValue> blst_p1_affine_in_g1(
        const OnCurveG1AffineValue& input) noexcept;
    [[nodiscard]] std::optional<SubgroupG2AffineValue> blst_p2_affine_in_g2(
        const OnCurveG2AffineValue& input) noexcept;

    [[nodiscard]] SubgroupG1AffineValue blst_p1_affine_generator() noexcept;
    [[nodiscard]] SubgroupG2AffineValue blst_p2_affine_generator() noexcept;
    [[nodiscard]] OnCurveG1AffineValue as_on_curve(
        const SubgroupG1AffineValue& input) noexcept;
    [[nodiscard]] OnCurveG2AffineValue as_on_curve(
        const SubgroupG2AffineValue& input) noexcept;
    [[nodiscard]] FieldValidG1CandidateValue as_candidate(
        const OnCurveG1AffineValue& input) noexcept;
    [[nodiscard]] FieldValidG2CandidateValue as_candidate(
        const OnCurveG2AffineValue& input) noexcept;

    [[nodiscard]] OnCurveG1Value blst_p1_from_affine(
        const OnCurveG1AffineValue& input) noexcept;
    [[nodiscard]] OnCurveG2Value blst_p2_from_affine(
        const OnCurveG2AffineValue& input) noexcept;
    [[nodiscard]] SubgroupG1Value blst_p1_from_affine(
        const SubgroupG1AffineValue& input) noexcept;
    [[nodiscard]] SubgroupG2Value blst_p2_from_affine(
        const SubgroupG2AffineValue& input) noexcept;
    [[nodiscard]] OnCurveG1Value blst_p1_double(
        const OnCurveG1Value& input) noexcept;
    [[nodiscard]] OnCurveG2Value blst_p2_double(
        const OnCurveG2Value& input) noexcept;
    [[nodiscard]] OnCurveG1Value blst_p1_cneg(
        const OnCurveG1Value& input) noexcept;
    [[nodiscard]] OnCurveG2Value blst_p2_cneg(
        const OnCurveG2Value& input) noexcept;

    [[nodiscard]] std::array<std::byte, 48> blst_p1_affine_compress(
        const OnCurveG1AffineValue& input) noexcept;
    [[nodiscard]] std::array<std::byte, 96> blst_p2_affine_compress(
        const OnCurveG2AffineValue& input) noexcept;

    [[nodiscard]] std::optional<IncompleteOnCurveG1Pair> incomplete_add_pair(
        const OnCurveG1Value& left,
        const OnCurveG1Value& right) noexcept;
    [[nodiscard]] std::optional<IncompleteOnCurveG2Pair> incomplete_add_pair(
        const OnCurveG2Value& left,
        const OnCurveG2Value& right) noexcept;

    [[nodiscard]] GtElementValue blst_final_exp(
        const MillerValueStorage& input) noexcept;
    [[nodiscard]] MillerValueStorage blst_miller_loop(
        const SubgroupG2AffineValue& q,
        const SubgroupG1AffineValue& p) noexcept;

    // Deterministic retained fixtures for refinement rejection tests.
    [[nodiscard]] OnCurveG1AffineValue retained_non_subgroup_g1() noexcept;
    [[nodiscard]] OnCurveG2AffineValue retained_non_subgroup_g2() noexcept;

    [[nodiscard]] std::optional<PreparedFftDomainValue> prepare_fft_domain(
        std::size_t n) noexcept;
    [[nodiscard]] std::optional<LiveMiddleProductHandleValue>
    precompute_middle_product_b(
        std::span<const SubgroupG2Value> values,
        std::size_t n) noexcept;

    [[nodiscard]] std::optional<PairingUnsetPending> make_pairing_state(
        bool hash_or_encode,
        std::span<const std::byte> dst) noexcept;
    [[nodiscard]] PairingUnsetCommitted pairing_commit(
        PairingUnsetPending&& state) noexcept;
    [[nodiscard]] PairingG1Committed pairing_commit(
        PairingG1Pending&& state) noexcept;
    [[nodiscard]] PairingG2Committed pairing_commit(
        PairingG2Pending&& state) noexcept;

    /** Direct byte-boundary campaigns. A failed property aborts the process. */
    void run_decode_case(std::span<const std::byte> input);
    void run_bytes_case(std::span<const std::byte> input);
    void run_fields_case(std::span<const std::byte> input);
    void run_groups_case(std::span<const std::byte> input);
    void run_protocol_case(std::span<const std::byte> input);
    void run_pairing_state_case(std::span<const std::byte> input);
    void run_fft_middle_case(std::span<const std::byte> input);
    void run_msm_case(std::span<const std::byte> input);
    void run_uniq_case(std::span<const std::byte> input);
    void run_avx_case(std::span<const std::byte> input);
}
