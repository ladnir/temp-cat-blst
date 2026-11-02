// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs
//
// Experimental scalar G1 DFT for middle-product crossover measurements.
// This is deliberately not part of the public API or production dispatch.

#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <span>
#include <stdexcept>
#include <vector>

#include "blst/fft/middle_product.h"
#include "bindings/blst_scalar.h"

namespace blst::experimental {

/**
 * Scalar-backend dense DFT using CatBLST's reusable fixed-table MSM.
 *
 * The input-dependent affine multiple table is built once per transform and
 * reused for every output row. Fixed DFT scalar rows are prepared with the
 * immutable plan. Table construction and batch normalization are included in
 * the online operation; calls allocate nothing after construction.
 */
template<std::size_t InputSize, std::size_t Wbits>
class ScalarG1MiddleDftFixedTable {
    static_assert(InputSize != 0);
    static_assert((InputSize & (InputSize - 1)) == 0);
    static_assert(Wbits >= 2 && Wbits <= 14);

    static constexpr std::size_t kTransformSize = 2 * InputSize;
    static constexpr std::size_t kScalarBytes = 32;
    static constexpr std::size_t kTableEntries =
        InputSize * (std::size_t{1} << (Wbits - 1));

public:
    explicit ScalarG1MiddleDftFixedTable(std::span<const Fr> roots)
        : coefficient_bytes_(kTransformSize * InputSize * kScalarBytes),
          affine_points_(InputSize),
          table_(kTableEntries),
          scratch_(scratch_words()),
          point_descriptors_(InputSize, nullptr),
          affine_descriptors_(InputSize, nullptr),
          scalar_descriptors_(InputSize, nullptr)
    {
        if (roots.size() != kTransformSize) {
            throw std::invalid_argument(
                "experimental DFT requires exactly 2*n roots");
        }
        if (blst_p1s_mult_wbits_precompute_sizeof(Wbits, InputSize) !=
                table_.size() * sizeof(table_[0]) ||
            scratch_.empty()) {
            throw std::length_error(
                "experimental fixed-table DFT scratch size is invalid");
        }

        std::vector<std::array<byte, kScalarBytes>> canonical_roots(
            kTransformSize);
        for (std::size_t root = 0; root < kTransformSize; ++root) {
            Fr canonical_fr{};
            Fr::from_mont(canonical_fr, roots[root]);
            static_assert(sizeof(canonical_fr) == kScalarBytes);
            std::memcpy(
                canonical_roots[root].data(),
                &canonical_fr,
                sizeof(canonical_fr));
        }
        for (std::size_t output = 0; output < kTransformSize; ++output) {
            for (std::size_t input = 0; input < InputSize; ++input) {
                const std::size_t root =
                    (output * input) & (kTransformSize - 1);
                std::memcpy(
                    coefficient_bytes_.data() +
                        (output * InputSize + input) * kScalarBytes,
                    canonical_roots[root].data(),
                    kScalarBytes);
            }
        }
    }

    void transform(
        std::span<scalar::G1Point> out,
        std::span<const scalar::G1Point> input)
    {
        if (input.size() != InputSize || out.size() != kTransformSize) {
            throw std::invalid_argument(
                "experimental DFT input or output has the wrong length");
        }
        const std::uintptr_t out_begin =
            reinterpret_cast<std::uintptr_t>(out.data());
        const std::uintptr_t input_begin =
            reinterpret_cast<std::uintptr_t>(input.data());
        const std::uintptr_t out_end = out_begin + out.size_bytes();
        const std::uintptr_t input_end = input_begin + input.size_bytes();
        if (out_begin < input_end && input_begin < out_end) {
            throw std::invalid_argument(
                "experimental DFT output must not overlap its input");
        }
        if constexpr (InputSize == 1) {
            out[0] = input[0];
            out[1] = input[0];
            return;
        }

        point_descriptors_[0] = input.data();
        blst_p1s_to_affine(
            affine_points_.data(), point_descriptors_.data(), InputSize);
        affine_descriptors_[0] = affine_points_.data();
        blst_p1s_mult_wbits_precompute(
            table_.data(), Wbits, affine_descriptors_.data(), InputSize);

        for (std::size_t output = 0; output < kTransformSize; ++output) {
            scalar_descriptors_[0] = coefficient_bytes_.data() +
                output * InputSize * kScalarBytes;
            blst_p1s_mult_wbits(
                &out[output],
                table_.data(),
                Wbits,
                InputSize,
                scalar_descriptors_.data(),
                256,
                scratch_.data());
        }
    }

private:
    static std::size_t scratch_words()
    {
        const std::size_t bytes =
            blst_p1s_mult_wbits_scratch_sizeof(InputSize);
        if (bytes == 0 || bytes % sizeof(limb_t) != 0) {
            throw std::length_error(
                "experimental fixed-table DFT scratch size is invalid");
        }
        return bytes / sizeof(limb_t);
    }

    std::vector<byte> coefficient_bytes_;
    std::vector<scalar::G1AffinePoint> affine_points_;
    std::vector<scalar::G1AffinePoint> table_;
    std::vector<limb_t> scratch_;
    std::vector<const scalar::G1Point*> point_descriptors_;
    std::vector<const scalar::G1AffinePoint*> affine_descriptors_;
    std::vector<const byte*> scalar_descriptors_;
};

} // namespace blst::experimental
