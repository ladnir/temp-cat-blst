// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#include "fuzz/internal/adapter.h"
#include "fuzz/internal/adapter_access.h"

#include "bindings/blst_scalar.h"
#include "bindings/blst_fft.h"
#include "bindings/blst_middle_product.h"
#include "blst/scalar/curve/map_to_G.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstring>
#include <memory>
#include <new>
#include <vector>

namespace catblst::fuzz::adapter
{
    namespace
    {
        template<class Native, std::size_t Size, class Tag>
        [[nodiscard]] Native native(
            const DomainValue<OpaqueStorage<Size>, Tag>& input) noexcept
        {
            return internal::AdapterAccess::unwrap<Native>(input);
        }

        template<class Native, std::size_t Size, class Tag>
        [[nodiscard]] DomainValue<OpaqueStorage<Size>, Tag> domain(
            const Native& input) noexcept
        {
            return internal::AdapterAccess::wrap<Native, Size, Tag>(input);
        }

        [[nodiscard]] blst_fp reduce_fp_component(
            std::span<const std::byte> input) noexcept
        {
            std::array<blst_byte, 1024> readable{};
            const std::size_t length = std::min(input.size(), readable.size());
            if (length != 0)
                std::memcpy(readable.data(), input.data(), length);
            blst_fp output;
            ::blst_fp_from_le_bytes(&output, readable.data(), length);
            return output;
        }

        template<std::size_t Components>
        [[nodiscard]] std::array<blst_fp, Components> reduce_components(
            std::span<const std::byte> input) noexcept
        {
            std::array<blst_fp, Components> output{};
            const std::size_t chunk_size =
                input.empty() ? 0 : (input.size() + Components - 1) / Components;
            for (std::size_t index = 0; index < Components; ++index)
            {
                const std::size_t offset = std::min(input.size(), index * chunk_size);
                const std::size_t length =
                    std::min(chunk_size, input.size() - offset);
                output[index] = reduce_fp_component(input.subspan(offset, length));
            }
            return output;
        }
    }

    ReducedFrValue blst_scalar_from_le_bytes(
        std::span<const std::byte> input) noexcept
    {
        std::array<blst_byte, 1024> readable{};
        const std::size_t length = std::min(input.size(), readable.size());
        if (length != 0)
            std::memcpy(readable.data(), input.data(), length);
        blst_scalar scalar{};
        ::blst_scalar_from_le_bytes(&scalar, readable.data(), length);
        blst_fr output;
        ::blst_fr_from_scalar(&output, &scalar);
        return domain<blst_fr, sizeof(blst_fr), ReducedFrTag>(output);
    }

    ReducedFpValue blst_fp_from_le_bytes(
        std::span<const std::byte> input) noexcept
    {
        const blst_fp output = reduce_fp_component(input);
        return domain<blst_fp, sizeof(blst_fp), ReducedFpTag>(output);
    }

    ReducedFp2Value reduce_fp2(std::span<const std::byte> input) noexcept
    {
        const auto components = reduce_components<2>(input);
        blst_fp2 output{components[0], components[1]};
        return domain<blst_fp2, sizeof(blst_fp2), ReducedFp2Tag>(output);
    }

    ReducedFp12Value reduce_fp12(std::span<const std::byte> input) noexcept
    {
        const auto components = reduce_components<12>(input);
        blst_fp12 output{};
        std::size_t index = 0;
        for (auto& fp6 : output.fp6)
            for (auto& fp2 : fp6.fp2)
                for (auto& fp : fp2.fp)
                    fp = components[index++];
        return domain<blst_fp12, sizeof(blst_fp12), ReducedFp12Tag>(output);
    }

    std::optional<NonzeroFpValue> refine_nonzero(
        const ReducedFpValue& input) noexcept
    {
        const blst_fp value = native<blst_fp>(input);
        const blst_fp zero{};
        if (std::memcmp(&value, &zero, sizeof(value)) == 0)
            return std::nullopt;
        return domain<blst_fp, sizeof(blst_fp), NonzeroTag<ReducedFpTag>>(value);
    }

    std::optional<NonzeroFp2Value> refine_nonzero(
        const ReducedFp2Value& input) noexcept
    {
        const blst_fp2 value = native<blst_fp2>(input);
        const blst_fp2 zero{};
        if (std::memcmp(&value, &zero, sizeof(value)) == 0)
            return std::nullopt;
        return domain<blst_fp2, sizeof(blst_fp2), NonzeroTag<ReducedFp2Tag>>(value);
    }

    std::optional<NonzeroFp12Value> refine_nonzero(
        const ReducedFp12Value& input) noexcept
    {
        const blst_fp12 value = native<blst_fp12>(input);
        const blst_fp12 zero{};
        if (std::memcmp(&value, &zero, sizeof(value)) == 0)
            return std::nullopt;
        return domain<blst_fp12, sizeof(blst_fp12), NonzeroTag<ReducedFp12Tag>>(value);
    }

    FieldValidG1CandidateValue make_g1_candidate(
        const ReducedFpValue& x,
        const ReducedFpValue& y) noexcept
    {
        blst_p1_affine point{native<blst_fp>(x), native<blst_fp>(y)};
        return domain<blst_p1_affine, sizeof(blst_p1_affine), FieldValidG1CandidateTag>(point);
    }

    FieldValidG2CandidateValue make_g2_candidate(
        const ReducedFp2Value& x,
        const ReducedFp2Value& y) noexcept
    {
        blst_p2_affine point{native<blst_fp2>(x), native<blst_fp2>(y)};
        return domain<blst_p2_affine, sizeof(blst_p2_affine), FieldValidG2CandidateTag>(point);
    }

    std::optional<OnCurveG1AffineValue> blst_p1_deserialize(
        std::span<const std::byte> input) noexcept
    {
        std::array<blst_byte, 96> readable{};
        const std::size_t length = std::min(input.size(), readable.size());
        if (length != 0)
            std::memcpy(readable.data(), input.data(), length);
        blst_p1_affine output{};
        if (::blst_p1_deserialize(&output, readable.data()) != BLST_SUCCESS)
            return std::nullopt;
        return domain<blst_p1_affine, sizeof(blst_p1_affine), OnCurveG1AffineTag>(output);
    }

    std::optional<OnCurveG2AffineValue> blst_p2_deserialize(
        std::span<const std::byte> input) noexcept
    {
        std::array<blst_byte, 192> readable{};
        const std::size_t length = std::min(input.size(), readable.size());
        if (length != 0)
            std::memcpy(readable.data(), input.data(), length);
        blst_p2_affine output{};
        if (::blst_p2_deserialize(&output, readable.data()) != BLST_SUCCESS)
            return std::nullopt;
        return domain<blst_p2_affine, sizeof(blst_p2_affine), OnCurveG2AffineTag>(output);
    }

    std::optional<OnCurveG1AffineValue> blst_p1_affine_on_curve(
        const FieldValidG1CandidateValue& input) noexcept
    {
        const blst_p1_affine point = native<blst_p1_affine>(input);
        if (!::blst_p1_affine_on_curve(&point))
            return std::nullopt;
        return domain<blst_p1_affine, sizeof(blst_p1_affine), OnCurveG1AffineTag>(point);
    }

    std::optional<OnCurveG2AffineValue> blst_p2_affine_on_curve(
        const FieldValidG2CandidateValue& input) noexcept
    {
        const blst_p2_affine point = native<blst_p2_affine>(input);
        if (!::blst_p2_affine_on_curve(&point))
            return std::nullopt;
        return domain<blst_p2_affine, sizeof(blst_p2_affine), OnCurveG2AffineTag>(point);
    }

    std::optional<SubgroupG1AffineValue> blst_p1_affine_in_g1(
        const OnCurveG1AffineValue& input) noexcept
    {
        const blst_p1_affine point = native<blst_p1_affine>(input);
        if (!::blst_p1_affine_in_g1(&point))
            return std::nullopt;
        return domain<blst_p1_affine, sizeof(blst_p1_affine), SubgroupG1AffineTag>(point);
    }

    std::optional<SubgroupG2AffineValue> blst_p2_affine_in_g2(
        const OnCurveG2AffineValue& input) noexcept
    {
        const blst_p2_affine point = native<blst_p2_affine>(input);
        if (!::blst_p2_affine_in_g2(&point))
            return std::nullopt;
        return domain<blst_p2_affine, sizeof(blst_p2_affine), SubgroupG2AffineTag>(point);
    }

    SubgroupG1AffineValue blst_p1_affine_generator() noexcept
    {
        const blst_p1_affine point = *::blst_p1_affine_generator();
        return domain<blst_p1_affine, sizeof(blst_p1_affine), SubgroupG1AffineTag>(point);
    }

    SubgroupG2AffineValue blst_p2_affine_generator() noexcept
    {
        const blst_p2_affine point = *::blst_p2_affine_generator();
        return domain<blst_p2_affine, sizeof(blst_p2_affine), SubgroupG2AffineTag>(point);
    }

    OnCurveG1AffineValue as_on_curve(const SubgroupG1AffineValue& input) noexcept
    {
        const blst_p1_affine point = native<blst_p1_affine>(input);
        return domain<blst_p1_affine, sizeof(blst_p1_affine), OnCurveG1AffineTag>(point);
    }

    OnCurveG2AffineValue as_on_curve(const SubgroupG2AffineValue& input) noexcept
    {
        const blst_p2_affine point = native<blst_p2_affine>(input);
        return domain<blst_p2_affine, sizeof(blst_p2_affine), OnCurveG2AffineTag>(point);
    }

    FieldValidG1CandidateValue as_candidate(const OnCurveG1AffineValue& input) noexcept
    {
        const blst_p1_affine point = native<blst_p1_affine>(input);
        return domain<blst_p1_affine, sizeof(blst_p1_affine), FieldValidG1CandidateTag>(point);
    }

    FieldValidG2CandidateValue as_candidate(const OnCurveG2AffineValue& input) noexcept
    {
        const blst_p2_affine point = native<blst_p2_affine>(input);
        return domain<blst_p2_affine, sizeof(blst_p2_affine), FieldValidG2CandidateTag>(point);
    }

    OnCurveG1Value blst_p1_from_affine(const OnCurveG1AffineValue& input) noexcept
    {
        const blst_p1_affine affine = native<blst_p1_affine>(input);
        blst_p1 output;
        ::blst_p1_from_affine(&output, &affine);
        return domain<blst_p1, sizeof(blst_p1), OnCurveG1Tag>(output);
    }

    OnCurveG2Value blst_p2_from_affine(const OnCurveG2AffineValue& input) noexcept
    {
        const blst_p2_affine affine = native<blst_p2_affine>(input);
        blst_p2 output;
        ::blst_p2_from_affine(&output, &affine);
        return domain<blst_p2, sizeof(blst_p2), OnCurveG2Tag>(output);
    }

    SubgroupG1Value blst_p1_from_affine(const SubgroupG1AffineValue& input) noexcept
    {
        const blst_p1_affine affine = native<blst_p1_affine>(input);
        blst_p1 output;
        ::blst_p1_from_affine(&output, &affine);
        return domain<blst_p1, sizeof(blst_p1), SubgroupG1Tag>(output);
    }

    SubgroupG2Value blst_p2_from_affine(const SubgroupG2AffineValue& input) noexcept
    {
        const blst_p2_affine affine = native<blst_p2_affine>(input);
        blst_p2 output;
        ::blst_p2_from_affine(&output, &affine);
        return domain<blst_p2, sizeof(blst_p2), SubgroupG2Tag>(output);
    }

    OnCurveG1Value blst_p1_double(const OnCurveG1Value& input) noexcept
    {
        const blst_p1 point = native<blst_p1>(input);
        blst_p1 output;
        ::blst_p1_double(&output, &point);
        return domain<blst_p1, sizeof(blst_p1), OnCurveG1Tag>(output);
    }

    OnCurveG2Value blst_p2_double(const OnCurveG2Value& input) noexcept
    {
        const blst_p2 point = native<blst_p2>(input);
        blst_p2 output;
        ::blst_p2_double(&output, &point);
        return domain<blst_p2, sizeof(blst_p2), OnCurveG2Tag>(output);
    }

    OnCurveG1Value blst_p1_cneg(const OnCurveG1Value& input) noexcept
    {
        blst_p1 output = native<blst_p1>(input);
        ::blst_p1_cneg(&output, true);
        return domain<blst_p1, sizeof(blst_p1), OnCurveG1Tag>(output);
    }

    OnCurveG2Value blst_p2_cneg(const OnCurveG2Value& input) noexcept
    {
        blst_p2 output = native<blst_p2>(input);
        ::blst_p2_cneg(&output, true);
        return domain<blst_p2, sizeof(blst_p2), OnCurveG2Tag>(output);
    }

    std::array<std::byte, 48> blst_p1_affine_compress(
        const OnCurveG1AffineValue& input) noexcept
    {
        const blst_p1_affine point = native<blst_p1_affine>(input);
        std::array<std::byte, 48> output{};
        std::array<blst_byte, 48> encoded{};
        ::blst_p1_affine_compress(encoded.data(), &point);
        std::memcpy(output.data(), encoded.data(), output.size());
        return output;
    }

    std::array<std::byte, 96> blst_p2_affine_compress(
        const OnCurveG2AffineValue& input) noexcept
    {
        const blst_p2_affine point = native<blst_p2_affine>(input);
        std::array<std::byte, 96> output{};
        std::array<blst_byte, 96> encoded{};
        ::blst_p2_affine_compress(encoded.data(), &point);
        std::memcpy(output.data(), encoded.data(), output.size());
        return output;
    }

    std::optional<IncompleteOnCurveG1Pair> incomplete_add_pair(
        const OnCurveG1Value& left,
        const OnCurveG1Value& right) noexcept
    {
        const blst_p1 a = native<blst_p1>(left);
        blst_p1 b = native<blst_p1>(right);
        if (::blst_p1_is_inf(&a) || ::blst_p1_is_inf(&b) ||
            ::blst_p1_is_equal(&a, &b))
            return std::nullopt;
        ::blst_p1_cneg(&b, true);
        if (::blst_p1_is_equal(&a, &b))
            return std::nullopt;
        return internal::AdapterAccess::make<
            std::array<OnCurveG1Value, 2>,
            IncompleteAddPairTag<OnCurveG1Tag>>({left, right});
    }

    std::optional<IncompleteOnCurveG2Pair> incomplete_add_pair(
        const OnCurveG2Value& left,
        const OnCurveG2Value& right) noexcept
    {
        const blst_p2 a = native<blst_p2>(left);
        blst_p2 b = native<blst_p2>(right);
        if (::blst_p2_is_inf(&a) || ::blst_p2_is_inf(&b) ||
            ::blst_p2_is_equal(&a, &b))
            return std::nullopt;
        ::blst_p2_cneg(&b, true);
        if (::blst_p2_is_equal(&a, &b))
            return std::nullopt;
        return internal::AdapterAccess::make<
            std::array<OnCurveG2Value, 2>,
            IncompleteAddPairTag<OnCurveG2Tag>>({left, right});
    }

    GtElementValue blst_final_exp(const MillerValueStorage& input) noexcept
    {
        const blst_fp12 value = native<blst_fp12>(input);
        blst_fp12 output;
        ::blst_final_exp(&output, &value);
        return domain<blst_fp12, sizeof(blst_fp12), GtElementTag>(output);
    }

    MillerValueStorage blst_miller_loop(
        const SubgroupG2AffineValue& q,
        const SubgroupG1AffineValue& p) noexcept
    {
        const blst_p2_affine q_native = native<blst_p2_affine>(q);
        const blst_p1_affine p_native = native<blst_p1_affine>(p);
        blst_fp12 output;
        ::blst_miller_loop(&output, &q_native, &p_native);
        return domain<blst_fp12, sizeof(blst_fp12), MillerValueTag>(output);
    }

    OnCurveG1AffineValue retained_non_subgroup_g1() noexcept
    {
        using blst::scalar::G1AffinePoint;
        using blst::scalar::G1Point;
        G1Point isogenous{};
        G1Point point{};
        blst::scalar::map_to_isogenous<G1Point>(
            isogenous, G1AffinePoint::generator().x);
        blst::scalar::isogeny_map_to<G1Point>(point, isogenous);
        G1AffinePoint affine;
        G1Point::to_affine(affine, point);
        return domain<blst_p1_affine, sizeof(blst_p1_affine), OnCurveG1AffineTag>(affine);
    }

    OnCurveG2AffineValue retained_non_subgroup_g2() noexcept
    {
        using blst::scalar::G2AffinePoint;
        using blst::scalar::G2Point;
        G2Point isogenous{};
        G2Point point{};
        blst::scalar::map_to_isogenous<G2Point>(
            isogenous, G2AffinePoint::generator().x);
        blst::scalar::isogeny_map_to<G2Point>(point, isogenous);
        G2AffinePoint affine;
        G2Point::to_affine(affine, point);
        return domain<blst_p2_affine, sizeof(blst_p2_affine), OnCurveG2AffineTag>(affine);
    }

    std::optional<PreparedFftDomainValue> prepare_fft_domain(
        std::size_t n) noexcept
    {
        const std::size_t length = ::blst_fr_fft_domain_len(n);
        if (length == 0)
            return std::nullopt;
        try
        {
            std::vector<blst_fr> native_domain(length);
            if (!::blst_fr_fft_domain_init(native_domain.data(), length, n))
                return std::nullopt;

            PreparedFftStorage storage;
            storage.elements.resize(length);
            storage.n = n;
            for (std::size_t index = 0; index < length; ++index)
            {
                std::memcpy(
                    storage.elements[index].bytes.data(),
                    &native_domain[index],
                    sizeof(blst_fr));
            }
            return internal::AdapterAccess::make<
                PreparedFftStorage,
                PreparedFftDomainTag>(std::move(storage));
        }
        catch (...)
        {
            return std::nullopt;
        }
    }

    std::optional<LiveMiddleProductHandleValue> precompute_middle_product_b(
        std::span<const SubgroupG2Value> values,
        std::size_t n) noexcept
    {
        if (n == 0 || n > values.size() || values.size() - n != n)
            return std::nullopt;
        try
        {
            std::vector<blst_p2> native_values;
            native_values.reserve(values.size());
            for (const auto& value : values)
                native_values.push_back(native<blst_p2>(value));
            blst_middle_product_precomputed_b* raw =
                ::blst_middle_product_precompute_b(native_values.data(), n);
            if (raw == nullptr)
                return std::nullopt;

            MiddleProductHandleStorage storage{
                std::shared_ptr<void>(
                    raw,
                    [](void* pointer) noexcept {
                        ::blst_middle_product_precomputed_b_free(
                            static_cast<blst_middle_product_precomputed_b*>(pointer));
                    }),
                n};
            return internal::AdapterAccess::make<
                MiddleProductHandleStorage,
                LiveMiddleProductHandleTag>(std::move(storage));
        }
        catch (...)
        {
            return std::nullopt;
        }
    }

    std::optional<PairingUnsetPending> make_pairing_state(
        bool hash_or_encode,
        std::span<const std::byte> dst) noexcept
    {
        const std::size_t size = ::blst_pairing_sizeof();
        const std::size_t alignment = ::blst_pairing_alignof();
        if (size == 0 || alignment == 0)
            return std::nullopt;
        try
        {
            std::vector<std::uint8_t> owned_dst(dst.size());
            if (!dst.empty())
                std::memcpy(owned_dst.data(), dst.data(), dst.size());
            void* memory = ::operator new(size, std::align_val_t(alignment));
            PairingStorage storage(
                memory,
                alignment,
                std::move(owned_dst),
                size);
            const auto* dst_pointer = storage.dst.empty() ? nullptr : storage.dst.data();
            ::blst_pairing_init(
                static_cast<blst_pairing*>(storage.context.get()),
                hash_or_encode,
                dst_pointer,
                storage.dst.size());
            return internal::AdapterAccess::make<
                PairingStorage,
                PairingStateTag<
                    PairingOrientation::unset,
                    PairingPendingState::pending>>(std::move(storage));
        }
        catch (...)
        {
            return std::nullopt;
        }
    }

    namespace
    {
        template<PairingOrientation Orientation>
        [[nodiscard]] PairingStateValue<
            Orientation,
            PairingPendingState::committed>
        pairing_commit_impl(PairingStateValue<
            Orientation,
            PairingPendingState::pending>&& state) noexcept
        {
            PairingStorage storage = internal::AdapterAccess::take(std::move(state));
            ::blst_pairing_commit(
                static_cast<blst_pairing*>(storage.context.get()));
            return internal::AdapterAccess::make<
                PairingStorage,
                PairingStateTag<Orientation, PairingPendingState::committed>>(
                    std::move(storage));
        }
    }

    PairingUnsetCommitted pairing_commit(PairingUnsetPending&& state) noexcept
    {
        return pairing_commit_impl(std::move(state));
    }

    PairingG1Committed pairing_commit(PairingG1Pending&& state) noexcept
    {
        return pairing_commit_impl(std::move(state));
    }

    PairingG2Committed pairing_commit(PairingG2Pending&& state) noexcept
    {
        return pairing_commit_impl(std::move(state));
    }
}
