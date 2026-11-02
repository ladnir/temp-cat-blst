// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs
//
// Experimental dense G1 DFT for middle-product crossover measurements.
// This is deliberately not part of the public API or production dispatch.

#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <span>
#include <stdexcept>
#include <vector>

#include "blst/avx/curve/point.h"
#include "blst/avx/curve/multi_scalar.h"
#include "blst/avx/curve/scalar_mult.h"
#include "blst/fft/experimental_scalar_g1_dft_msm.h"
#include "blst/fft/middle_product.h"
#include "bindings/blst_scalar.h"

namespace blst::experimental {

/**
 * Dense forward DFT specialized to the middle-product input shape.
 *
 * The input contains `InputSize` native AVX G1 points. The mathematical input
 * is padded with `InputSize` identities, so the transform size is
 * `2 * InputSize`. Every output is therefore an MSM over the same
 * `InputSize` non-padding bases:
 *
 * `out[u] = sum_v roots[(u * v) mod (2 * InputSize)] * input[v]`.
 *
 * The root matrix is immutable and constructed once. Every call rebuilds the
 * input-dependent odd-multiple tables exactly once and reuses them across all
 * output rows. That table construction is part of the online operation and is
 * intentionally included in benchmarks. Calls allocate nothing after object
 * construction.
 *
 * This kernel is variable-time in its public root scalars. Inputs must already
 * be on-curve prime-order points. `out` and `input` must not overlap. The class
 * is mutable workspace and must not be used concurrently.
 */
template<std::size_t InputSize, int Width = 5>
class AvxG1MiddleDftMsm {
    static_assert(InputSize != 0);
    static_assert((InputSize & (InputSize - 1)) == 0);
    static_assert(Width >= 4 && Width <= 8);

    static constexpr std::size_t kTransformSize = 2 * InputSize;
    static constexpr int kChunkSize =
        InputSize >= 32 ? 32 : static_cast<int>(InputSize);
    static constexpr std::size_t kChunkCount = InputSize / kChunkSize;
    static constexpr int kWidth = Width;
    static constexpr std::size_t kTableSize = std::size_t{1} << (kWidth - 2);

    using Curve = avx::avx_curve_detail::CurveG1;
    using Ring = avx::avx_curve_detail::CurveRing;
    using Glv = avx::avx_curve_detail::CurveG1GLV;
    using RawPoint = Curve::ProjPoint;
    using Tables = typename Glv::template BatchWnafTables<kChunkSize, kWidth>;
    using ScalarBytes = std::array<byte, 32>;
    using ComponentDigits = std::array<std::int8_t, 129>;

    struct RootRecoding {
        std::array<ComponentDigits, 2> digits{};
        std::size_t used = 0;
    };

public:
    explicit AvxG1MiddleDftMsm(std::span<const Fr> roots)
        : coefficient_indices_(kTransformSize * InputSize),
          root_recodings_(kTransformSize),
          tables_(kChunkCount)
    {
        if (roots.size() != kTransformSize) {
            throw std::invalid_argument(
                "experimental DFT requires exactly 2*n roots");
        }
        for (std::size_t root = 0; root < kTransformSize; ++root) {
            ScalarBytes canonical{};
            Fr canonical_fr{};
            Fr::from_mont(canonical_fr, roots[root]);
            static_assert(sizeof(canonical_fr) == sizeof(canonical));
            std::memcpy(
                canonical.data(), &canonical_fr, sizeof(canonical_fr));
            const avx::div_by_zz_result split = avx::div_by_zz(canonical);
            RootRecoding& recoding = root_recodings_[root];
            const std::size_t low = recode_component(
                recoding.digits[0].data(), split.first.data());
            const std::size_t high = recode_component(
                recoding.digits[1].data(), split.second.data());
            recoding.used = std::max(low, high);
        }
        for (std::size_t output = 0; output < kTransformSize; ++output) {
            for (std::size_t input = 0; input < InputSize; ++input) {
                coefficient_indices_[output * InputSize + input] =
                    static_cast<std::uint32_t>(
                        (output * input) & (kTransformSize - 1));
            }
        }
    }

    static constexpr std::size_t input_size() noexcept
    {
        return InputSize;
    }

    static constexpr std::size_t output_size() noexcept
    {
        return kTransformSize;
    }

    void transform(
        std::span<avx::G1Point> out,
        std::span<const avx::G1Point> input)
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

        prepare_tables(input);
        transform_prepared(out);
    }

    /** Experimental benchmark hook; the caller must preserve these bases. */
    void prepare_bases(std::span<const avx::G1Point> input)
    {
        if (input.size() != InputSize) {
            throw std::invalid_argument(
                "experimental DFT input has the wrong length");
        }
        if constexpr (InputSize != 1) {
            prepare_tables(input);
        }
    }

    /** Evaluate after `prepare_bases`; used only to bound table-build cost. */
    void transform_prepared(std::span<avx::G1Point> out)
    {
        if (out.size() != kTransformSize) {
            throw std::invalid_argument(
                "experimental DFT output has the wrong length");
        }
        if constexpr (InputSize == 1) {
            throw std::logic_error(
                "the one-point experimental DFT needs no prepared tables");
        }
        for (std::size_t output = 0; output < kTransformSize; ++output) {
            const std::uint32_t* const row =
                coefficient_indices_.data() + output * InputSize;
            avx::G1Point sum = avx::G1Point::infinity();
            for (std::size_t chunk = 0; chunk < kChunkCount; ++chunk) {
                const RawPoint partial = evaluate_chunk(
                    tables_[chunk], row + chunk * kChunkSize);
                const avx::G1Point native_partial{
                    partial.x, partial.y, partial.z};
                if (chunk == 0) {
                    sum = native_partial;
                } else {
                    avx::G1Point::dadd(sum, sum, native_partial);
                }
            }
            out[output] = sum;
        }
    }

private:
    static std::size_t recode_component(
        std::int8_t digits[129],
        const byte scalar_bytes[16])
    {
        if constexpr (Width <= 7) {
            return scalar::recode_wnaf_vartime<128, Width>(
                digits, scalar_bytes);
        } else {
            static_assert(Width == 8);
            constexpr std::size_t kScalarLimbs =
                (128 + BLST_LIMB_SIZE - 1) / BLST_LIMB_SIZE;
            constexpr std::size_t kLimbCount = kScalarLimbs + 1;
            limb_t value[kLimbCount] = {};
            limbs_from_le_bytes(value, scalar_bytes, 16);
            std::fill_n(digits, 129, std::int8_t{0});

            std::size_t used = 0;
            while (used <= 128) {
                limb_t nonzero = 0;
                for (const limb_t limb : value) {
                    nonzero |= limb;
                }
                if (nonzero == 0) {
                    break;
                }
                if ((value[0] & 1U) != 0) {
                    const unsigned raw =
                        static_cast<unsigned>(value[0] & 0xffU);
                    const int digit = raw >= 128U
                        ? static_cast<int>(raw) - 256
                        : static_cast<int>(raw);
                    digits[used] = static_cast<std::int8_t>(digit);
                    if (digit > 0) {
                        value[0] -= static_cast<limb_t>(digit);
                    } else {
                        limb_t add = static_cast<limb_t>(-digit);
                        for (std::size_t limb = 0;
                             limb < kLimbCount && add != 0;
                             ++limb) {
                            const limb_t before = value[limb];
                            value[limb] += add;
                            add = value[limb] < before ? 1 : 0;
                        }
                    }
                }
                for (std::size_t limb = 0;
                     limb + 1 < kLimbCount;
                     ++limb) {
                    value[limb] =
                        (value[limb] >> 1) |
                        (value[limb + 1] << (BLST_LIMB_SIZE - 1));
                }
                value[kLimbCount - 1] >>= 1;
                ++used;
            }
            return used;
        }
    }

    void prepare_tables(std::span<const avx::G1Point> input)
    {
        const Ring& ring = avx::avx_curve_detail::curve_ring();
        const Glv& glv = avx::avx_curve_detail::curve_g1_glv();
        for (std::size_t chunk = 0; chunk < kChunkCount; ++chunk) {
            Tables& tables = tables_[chunk];
            for (int lane = 0; lane < kChunkSize; ++lane) {
                const std::size_t input_index =
                    chunk * kChunkSize + static_cast<std::size_t>(lane);
                const avx::G1Point& point =
                    input[input_index];
                const RawPoint raw{point.X, point.Y, point.Z};
                tables[2 * lane][0] = raw;
                // The first DFT column is identically one, so its input never
                // needs an odd-multiple or endomorphism table.
                if (input_index == 0) {
                    continue;
                }
                const RawPoint twice = glv.curve.double_point(raw, ring);
                for (std::size_t entry = 1; entry < kTableSize; ++entry) {
                    tables[2 * lane][entry] = glv.curve.add_point(
                        tables[2 * lane][entry - 1], twice, ring);
                }
                for (std::size_t entry = 0; entry < kTableSize; ++entry) {
                    const RawPoint& value = tables[2 * lane][entry];
                    tables[2 * lane + 1][entry] = {
                        ring.modmul(value.x, glv.beta),
                        ring.standard_negate(value.y),
                        value.z,
                    };
                }
            }
        }
    }

    RawPoint evaluate_chunk(
        const Tables& tables,
        const std::uint32_t* coefficient_row) const
    {
        constexpr int kComponentCount = 2 * kChunkSize;
        const Ring& ring = avx::avx_curve_detail::curve_ring();
        const Glv& glv = avx::avx_curve_detail::curve_g1_glv();
        std::size_t used = 0;
        RawPoint simple = glv.curve.zero(ring);
        bool simple_started = false;
        for (int lane = 0; lane < kChunkSize; ++lane) {
            const std::uint32_t root = coefficient_row[lane];
            if (root == 0 || root == InputSize) {
                RawPoint addend = tables[2 * lane][0];
                if (root == InputSize) {
                    addend.y = ring.standard_negate(addend.y);
                }
                if (!simple_started) {
                    simple = addend;
                    simple_started = true;
                } else {
                    simple = glv.curve.add_point(simple, addend, ring);
                }
                continue;
            }
            used = std::max(
                used,
                root_recodings_[root].used);
        }

        RawPoint result = glv.curve.zero(ring);
        bool started = false;
        for (std::size_t position = used; position != 0;) {
            --position;
            if (started) {
                result = glv.curve.double_point(result, ring);
            }
            for (int component = 0; component < kComponentCount; ++component) {
                const std::size_t lane =
                    static_cast<std::size_t>(component >> 1);
                const std::uint32_t root = coefficient_row[lane];
                if (root == 0 || root == InputSize) {
                    continue;
                }
                const std::size_t glv_component =
                    static_cast<std::size_t>(component & 1);
                const int digit = root_recodings_[root]
                    .digits[glv_component][position];
                if (digit == 0) {
                    continue;
                }
                const unsigned magnitude = digit < 0
                    ? static_cast<unsigned>(-digit)
                    : static_cast<unsigned>(digit);
                BLST_ASSERT((magnitude & 1U) != 0);
                BLST_ASSERT((magnitude >> 1) < kTableSize);
                RawPoint addend = tables[component][magnitude >> 1];
                if (digit < 0) {
                    addend.y = ring.standard_negate(addend.y);
                }
                if (!started) {
                    result = addend;
                    started = true;
                } else {
                    result = glv.curve.add_point(result, addend, ring);
                }
            }
        }
        if (!started) {
            return simple;
        }
        if (simple_started) {
            result = glv.curve.add_point(result, simple, ring);
        }
        return result;
    }

    std::vector<std::uint32_t> coefficient_indices_;
    std::vector<RootRecoding> root_recodings_;
    std::vector<Tables> tables_;
};

/**
 * Dense middle-product DFT using the production AVX GLV Pippenger kernel.
 *
 * The immutable plan stores pre-split scalar rows. Each transform batch-
 * normalizes the input bases and constructs their two GLV affine images once,
 * then reuses those images for all `2 * InputSize` Pippenger evaluations.
 * Input-dependent normalization and point preparation are included in the
 * measured online operation. Calls allocate nothing after construction.
 */
template<std::size_t InputSize>
class AvxG1MiddleDftPippenger {
    static_assert(InputSize != 0);
    static_assert((InputSize & (InputSize - 1)) == 0);

    static constexpr std::size_t kTransformSize = 2 * InputSize;
    static constexpr std::size_t kScalarBytes = 32;

    using ScalarBytes = std::array<byte, kScalarBytes>;

public:
    explicit AvxG1MiddleDftPippenger(std::span<const Fr> roots)
        : split_scalar_rows_(kTransformSize * InputSize * kScalarBytes),
          canonical_first_row_(InputSize * kScalarBytes),
          affine_points_(InputSize),
          affine_prefix_(InputSize),
          split_points_(2 * InputSize),
          prepare_split_scalars_(InputSize * kScalarBytes),
          buckets_(
              avx::avx_g1_mult_pippenger_glv_scratch_sizeof(InputSize) /
              sizeof(avx::AvxG1PippengerBucket)),
          bucket_used_(
              avx::avx_g1_mult_pippenger_glv_bucket_used_scratch_sizeof(
                  InputSize))
    {
        if (roots.size() != kTransformSize) {
            throw std::invalid_argument(
                "experimental DFT requires exactly 2*n roots");
        }
        if (buckets_.empty() || bucket_used_.empty()) {
            throw std::length_error(
                "experimental Pippenger DFT scratch size is invalid");
        }

        std::vector<ScalarBytes> canonical_roots(kTransformSize);
        std::vector<ScalarBytes> split_roots(kTransformSize);
        for (std::size_t root = 0; root < kTransformSize; ++root) {
            Fr canonical_fr{};
            Fr::from_mont(canonical_fr, roots[root]);
            static_assert(sizeof(canonical_fr) == kScalarBytes);
            std::memcpy(
                canonical_roots[root].data(),
                &canonical_fr,
                sizeof(canonical_fr));
            const avx::div_by_zz_result split =
                avx::div_by_zz(canonical_roots[root]);
            std::memcpy(
                split_roots[root].data(),
                split.first.data(),
                split.first.size());
            std::memcpy(
                split_roots[root].data() + split.first.size(),
                split.second.data(),
                split.second.size());
        }

        for (std::size_t input = 0; input < InputSize; ++input) {
            std::memcpy(
                canonical_first_row_.data() + input * kScalarBytes,
                canonical_roots[0].data(),
                kScalarBytes);
        }
        for (std::size_t output = 0; output < kTransformSize; ++output) {
            for (std::size_t input = 0; input < InputSize; ++input) {
                const std::size_t root =
                    (output * input) & (kTransformSize - 1);
                std::memcpy(
                    split_scalar_rows_.data() +
                        (output * InputSize + input) * kScalarBytes,
                    split_roots[root].data(),
                    kScalarBytes);
            }
        }
    }

    void transform(
        std::span<avx::G1Point> out,
        std::span<const avx::G1Point> input)
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

        avx::G1Point::to_affine_batch(
            affine_points_, input, affine_prefix_);
        avx::avx_g1_mult_pippenger_glv_prepare(
            split_points_.data(),
            prepare_split_scalars_.data(),
            affine_points_.data(),
            canonical_first_row_.data(),
            InputSize);

        for (std::size_t output = 0; output < kTransformSize; ++output) {
            avx::avx_g1_mult_pippenger_glv_prepared(
                out[output],
                split_points_.data(),
                split_scalar_rows_.data() +
                    output * InputSize * kScalarBytes,
                InputSize,
                buckets_.data(),
                0,
                bucket_used_.data());
        }
    }

private:
    std::vector<byte> split_scalar_rows_;
    std::vector<byte> canonical_first_row_;
    std::vector<avx::G1AffinePoint> affine_points_;
    std::vector<avx::G1Point::Field_type> affine_prefix_;
    std::vector<avx::G1AffinePoint> split_points_;
    std::vector<byte> prepare_split_scalars_;
    std::vector<avx::AvxG1PippengerBucket> buckets_;
    std::vector<std::uint8_t> bucket_used_;
};

} // namespace blst::experimental
