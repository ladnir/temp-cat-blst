// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#include "fuzz/internal/adapter.h"
#include "fuzz/internal/adapter_access.h"

#include "bindings/blst_scalar.h"
#include "bindings/blst_fft.h"
#include "bindings/blst_middle_product.h"
#include "blst/scalar/curve/map_to_G.h"
#include "fuzz/support/modulus_oracle.h"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <cstdio>
#include <cstddef>
#include <cstring>
#include <memory>
#include <new>
#include <source_location>
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

        void require(
            bool condition,
            const std::source_location location =
                std::source_location::current()) noexcept
        {
            if (!condition)
            {
                std::fprintf(
                    stderr,
                    "fuzz property failed at %s:%u\n",
                    location.file_name(),
                    location.line());
                std::abort();
            }
        }

        template<class T>
        [[nodiscard]] bool object_equal(const T& left, const T& right) noexcept
        {
            static_assert(std::is_trivially_copyable_v<T>);
            return std::memcmp(&left, &right, sizeof(T)) == 0;
        }

        template<std::size_t Size>
        [[nodiscard]] std::array<blst_byte, Size> readable_bytes(
            std::span<const std::byte> input) noexcept
        {
            std::array<blst_byte, Size> output{};
            const std::size_t copied = std::min(input.size(), output.size());
            if (copied != 0)
                std::memcpy(output.data(), input.data(), copied);
            return output;
        }

        struct ParsedHex final
        {
            bool syntax = false;
            std::size_t digits = 0;
            std::array<std::byte, 48> little_endian{};
        };

        [[nodiscard]] ParsedHex parse_hex(
            const std::array<blst_byte, 99>& text) noexcept
        {
            ParsedHex result;
            std::size_t length = 0;
            while (length < text.size() && text[length] != 0)
                ++length;
            std::size_t offset = 0;
            if (length >= 2 && text[0] == '0' &&
                (text[1] == 'x' || text[1] == 'X'))
                offset = 2;
            if (offset == length)
                return result;

            result.digits = length - offset;
            if (result.digits > result.little_endian.size() * 2)
                return result;
            for (std::size_t index = 0; index < result.digits; ++index)
            {
                const unsigned char character = text[offset + index];
                unsigned char nibble = 0;
                if (character >= '0' && character <= '9')
                    nibble = static_cast<unsigned char>(character - '0');
                else if (character >= 'a' && character <= 'f')
                    nibble = static_cast<unsigned char>(character - 'a' + 10);
                else if (character >= 'A' && character <= 'F')
                    nibble = static_cast<unsigned char>(character - 'A' + 10);
                else
                    return result;

                const std::size_t reverse = result.digits - 1 - index;
                const std::size_t byte = reverse / 2;
                if ((reverse & 1U) == 0)
                    result.little_endian[byte] |= std::byte{nibble};
                else
                    result.little_endian[byte] |= std::byte{
                        static_cast<unsigned char>(nibble << 4)};
            }
            result.syntax = true;
            return result;
        }

        template<std::size_t Size>
        [[nodiscard]] bool equal_prefix(
            const std::array<blst_byte, Size>& actual,
            const std::byte* expected) noexcept
        {
            return std::memcmp(actual.data(), expected, Size) == 0;
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

    void run_decode_case(std::span<const std::byte> input)
    {
        const auto g1_wire = readable_bytes<96>(input.first(
            std::min<std::size_t>(input.size(), 96)));
        const auto g2_wire = readable_bytes<192>(input.first(
            std::min<std::size_t>(input.size(), 192)));

        blst_p1_affine g1_sentinel;
        std::memset(&g1_sentinel, 0xa5, sizeof(g1_sentinel));
        blst_p1_affine g1_deserialized = g1_sentinel;
        const BLST_ERROR g1_deserialize_status =
            ::blst_p1_deserialize(&g1_deserialized, g1_wire.data());
        if (g1_deserialize_status != BLST_SUCCESS)
            require(object_equal(g1_deserialized, g1_sentinel));
        else
        {
            require(::blst_p1_affine_on_curve(&g1_deserialized));
            std::array<blst_byte, 96> serialized{};
            std::array<blst_byte, 48> compressed{};
            std::array<blst_byte, 96> readable_compressed{};
            ::blst_p1_affine_serialize(serialized.data(), &g1_deserialized);
            ::blst_p1_affine_compress(compressed.data(), &g1_deserialized);
            std::copy(compressed.begin(), compressed.end(), readable_compressed.begin());

            blst_p1_affine roundtrip{};
            require(::blst_p1_deserialize(&roundtrip, serialized.data()) == BLST_SUCCESS);
            require(::blst_p1_affine_is_equal(&g1_deserialized, &roundtrip));
            require(::blst_p1_uncompress(&roundtrip, compressed.data()) == BLST_SUCCESS);
            require(::blst_p1_affine_is_equal(&g1_deserialized, &roundtrip));
            require(::blst_p1_deserialize(&roundtrip, readable_compressed.data()) == BLST_SUCCESS);
            require(::blst_p1_affine_is_equal(&g1_deserialized, &roundtrip));

            blst_p1 projective{};
            ::blst_p1_from_affine(&projective, &g1_deserialized);
            require(::blst_p1_in_g1(&projective) ==
                    ::blst_p1_affine_in_g1(&g1_deserialized));
        }

        blst_p1_affine g1_uncompressed = g1_sentinel;
        const BLST_ERROR g1_uncompress_status =
            ::blst_p1_uncompress(&g1_uncompressed, g1_wire.data());
        if (g1_uncompress_status != BLST_SUCCESS)
            require(object_equal(g1_uncompressed, g1_sentinel));
        else
        {
            blst_p1_affine via_deserialize{};
            require(::blst_p1_deserialize(&via_deserialize, g1_wire.data()) == BLST_SUCCESS);
            require(::blst_p1_affine_is_equal(&g1_uncompressed, &via_deserialize));
        }

        blst_p1 aggregate_sentinel;
        std::memset(&aggregate_sentinel, 0x5a, sizeof(aggregate_sentinel));
        blst_p1 aggregate = aggregate_sentinel;
        const BLST_ERROR aggregate_g1_status =
            ::blst_aggregate_in_g1(&aggregate, nullptr, g1_wire.data());
        if (aggregate_g1_status != BLST_SUCCESS)
            require(object_equal(aggregate, aggregate_sentinel));
        else
        {
            require(::blst_p1_in_g1(&aggregate));
            require(g1_deserialize_status == BLST_SUCCESS);
            require(::blst_p1_affine_in_g1(&g1_deserialized));
            blst_p1 expected{};
            ::blst_p1_from_affine(&expected, &g1_deserialized);
            require(::blst_p1_is_equal(&aggregate, &expected));
        }

        blst_p2_affine g2_sentinel;
        std::memset(&g2_sentinel, 0xa5, sizeof(g2_sentinel));
        blst_p2_affine g2_deserialized = g2_sentinel;
        const BLST_ERROR g2_deserialize_status =
            ::blst_p2_deserialize(&g2_deserialized, g2_wire.data());
        if (g2_deserialize_status != BLST_SUCCESS)
            require(object_equal(g2_deserialized, g2_sentinel));
        else
        {
            require(::blst_p2_affine_on_curve(&g2_deserialized));
            std::array<blst_byte, 192> serialized{};
            std::array<blst_byte, 96> compressed{};
            std::array<blst_byte, 192> readable_compressed{};
            ::blst_p2_affine_serialize(serialized.data(), &g2_deserialized);
            ::blst_p2_affine_compress(compressed.data(), &g2_deserialized);
            std::copy(compressed.begin(), compressed.end(), readable_compressed.begin());

            blst_p2_affine roundtrip{};
            require(::blst_p2_deserialize(&roundtrip, serialized.data()) == BLST_SUCCESS);
            require(::blst_p2_affine_is_equal(&g2_deserialized, &roundtrip));
            require(::blst_p2_uncompress(&roundtrip, compressed.data()) == BLST_SUCCESS);
            require(::blst_p2_affine_is_equal(&g2_deserialized, &roundtrip));
            require(::blst_p2_deserialize(&roundtrip, readable_compressed.data()) == BLST_SUCCESS);
            require(::blst_p2_affine_is_equal(&g2_deserialized, &roundtrip));

            blst_p2 projective{};
            ::blst_p2_from_affine(&projective, &g2_deserialized);
            require(::blst_p2_in_g2(&projective) ==
                    ::blst_p2_affine_in_g2(&g2_deserialized));
        }

        blst_p2_affine g2_uncompressed = g2_sentinel;
        const BLST_ERROR g2_uncompress_status =
            ::blst_p2_uncompress(&g2_uncompressed, g2_wire.data());
        if (g2_uncompress_status != BLST_SUCCESS)
            require(object_equal(g2_uncompressed, g2_sentinel));
        else
        {
            blst_p2_affine via_deserialize{};
            require(::blst_p2_deserialize(&via_deserialize, g2_wire.data()) == BLST_SUCCESS);
            require(::blst_p2_affine_is_equal(&g2_uncompressed, &via_deserialize));
        }

        blst_p2 aggregate_g2_sentinel;
        std::memset(&aggregate_g2_sentinel, 0x5a, sizeof(aggregate_g2_sentinel));
        blst_p2 aggregate_g2 = aggregate_g2_sentinel;
        const BLST_ERROR aggregate_g2_status =
            ::blst_aggregate_in_g2(&aggregate_g2, nullptr, g2_wire.data());
        if (aggregate_g2_status != BLST_SUCCESS)
            require(object_equal(aggregate_g2, aggregate_g2_sentinel));
        else
        {
            require(::blst_p2_in_g2(&aggregate_g2));
            require(g2_deserialize_status == BLST_SUCCESS);
            require(::blst_p2_affine_in_g2(&g2_deserialized));
            blst_p2 expected{};
            ::blst_p2_from_affine(&expected, &g2_deserialized);
            require(::blst_p2_is_equal(&aggregate_g2, &expected));
        }
    }

    void run_bytes_case(std::span<const std::byte> input)
    {
        input = input.first(std::min<std::size_t>(input.size(), 1024));
        const auto readable = readable_bytes<1024>(input);
        const blst_byte* bytes = input.empty() ? nullptr : readable.data();

        const ModulusOracle<8> scalar_oracle(scalar_modulus32);
        const auto scalar_le_expected =
            scalar_oracle.reduce(input, ByteOrder::little);
        const auto scalar_be_expected =
            scalar_oracle.reduce(input, ByteOrder::big);
        const auto scalar_le_bytes =
            ModulusOracle<8>::to_le_bytes(scalar_le_expected);
        const auto scalar_be_bytes =
            ModulusOracle<8>::to_le_bytes(scalar_be_expected);

        blst_scalar scalar_le{};
        blst_scalar scalar_be{};
        const bool le_nonzero =
            ::blst_scalar_from_le_bytes(&scalar_le, bytes, input.size());
        const bool be_nonzero =
            ::blst_scalar_from_be_bytes(&scalar_be, bytes, input.size());
        std::array<blst_byte, 32> scalar_export{};
        ::blst_lendian_from_scalar(scalar_export.data(), &scalar_le);
        require(equal_prefix(scalar_export, scalar_le_bytes.data()));
        require(le_nonzero == !ModulusOracle<8>::is_zero(scalar_le_expected));
        require(::blst_scalar_fr_check(&scalar_le));
        require(::blst_sk_check(&scalar_le) == le_nonzero);
        ::blst_lendian_from_scalar(scalar_export.data(), &scalar_be);
        require(equal_prefix(scalar_export, scalar_be_bytes.data()));
        require(be_nonzero == !ModulusOracle<8>::is_zero(scalar_be_expected));
        require(::blst_scalar_fr_check(&scalar_be));
        require(::blst_sk_check(&scalar_be) == be_nonzero);

        blst_scalar empty_scalar;
        std::memset(&empty_scalar, 0xa5, sizeof(empty_scalar));
        require(!::blst_scalar_from_le_bytes(&empty_scalar, nullptr, 0));
        std::array<blst_byte, 32> zero_export{};
        ::blst_lendian_from_scalar(scalar_export.data(), &empty_scalar);
        require(scalar_export == zero_export);
        std::memset(&empty_scalar, 0xa5, sizeof(empty_scalar));
        require(!::blst_scalar_from_be_bytes(&empty_scalar, nullptr, 0));
        ::blst_lendian_from_scalar(scalar_export.data(), &empty_scalar);
        require(scalar_export == zero_export);

        const auto fixed_scalar = readable_bytes<32>(input);
        blst_scalar fixed{};
        ::blst_scalar_from_lendian(&fixed, fixed_scalar.data());
        ::blst_lendian_from_scalar(scalar_export.data(), &fixed);
        require(scalar_export == fixed_scalar);
        ::blst_scalar_from_bendian(&fixed, fixed_scalar.data());
        std::array<blst_byte, 32> big_export{};
        ::blst_bendian_from_scalar(big_export.data(), &fixed);
        require(big_export == fixed_scalar);

        std::array<std::uint32_t, 8> words32{};
        std::array<std::uint64_t, 4> words64{};
        for (std::size_t index = 0; index < fixed_scalar.size(); ++index)
        {
            words32[index / 4] |= std::to_integer<std::uint32_t>(
                std::byte{fixed_scalar[index]})
                << ((index % 4) * 8);
            words64[index / 8] |= std::to_integer<std::uint64_t>(
                std::byte{fixed_scalar[index]})
                << ((index % 8) * 8);
        }
        ::blst_scalar_from_uint32(&fixed, words32.data());
        std::array<std::uint32_t, 8> words32_out{};
        ::blst_uint32_from_scalar(words32_out.data(), &fixed);
        require(words32_out == words32);
        ::blst_scalar_from_uint64(&fixed, words64.data());
        std::array<std::uint64_t, 4> words64_out{};
        ::blst_uint64_from_scalar(words64_out.data(), &fixed);
        require(words64_out == words64);

        blst_fr fr{};
        ::blst_fr_from_scalar(&fr, &scalar_le);
        blst_scalar scalar_roundtrip{};
        ::blst_scalar_from_fr(&scalar_roundtrip, &fr);
        require(object_equal(scalar_roundtrip, scalar_le));
        const std::array<blst_fr, 2> fr_batch{fr, fr};
        std::array<blst_scalar, 2> scalar_batch{};
        require(::blst_scalars_from_fr(
            scalar_batch.data(), fr_batch.data(), fr_batch.size()));
        require(object_equal(scalar_batch[0], scalar_le));
        require(object_equal(scalar_batch[1], scalar_le));
        require(::blst_scalars_from_fr(nullptr, nullptr, 0));

        const ModulusOracle<12> fp_oracle(fp_modulus32);
        const auto fp_le_expected = fp_oracle.reduce(input, ByteOrder::little);
        const auto fp_be_expected = fp_oracle.reduce(input, ByteOrder::big);
        const auto fp_le_bytes = ModulusOracle<12>::to_le_bytes(fp_le_expected);
        const auto fp_be_bytes = ModulusOracle<12>::to_le_bytes(fp_be_expected);
        blst_fp fp_le{};
        blst_fp fp_be{};
        ::blst_fp_from_le_bytes(&fp_le, bytes, input.size());
        ::blst_fp_from_be_bytes(&fp_be, bytes, input.size());
        std::array<blst_byte, 48> fp_export{};
        ::blst_lendian_from_fp(fp_export.data(), &fp_le);
        require(equal_prefix(fp_export, fp_le_bytes.data()));
        ::blst_lendian_from_fp(fp_export.data(), &fp_be);
        require(equal_prefix(fp_export, fp_be_bytes.data()));

        blst_fp empty_fp;
        std::memset(&empty_fp, 0xa5, sizeof(empty_fp));
        ::blst_fp_from_le_bytes(&empty_fp, nullptr, 0);
        ::blst_lendian_from_fp(fp_export.data(), &empty_fp);
        require(fp_export == std::array<blst_byte, 48>{});
        std::memset(&empty_fp, 0xa5, sizeof(empty_fp));
        ::blst_fp_from_be_bytes(&empty_fp, nullptr, 0);
        ::blst_lendian_from_fp(fp_export.data(), &empty_fp);
        require(fp_export == std::array<blst_byte, 48>{});

        const auto fixed_fp = readable_bytes<48>(input);
        ::blst_fp_from_lendian(&fp_le, fixed_fp.data());
        ::blst_lendian_from_fp(fp_export.data(), &fp_le);
        const auto fixed_fp_le_expected = fp_oracle.reduce(
            std::as_bytes(std::span(fixed_fp)), ByteOrder::little);
        const auto fixed_fp_le_bytes =
            ModulusOracle<12>::to_le_bytes(fixed_fp_le_expected);
        require(equal_prefix(fp_export, fixed_fp_le_bytes.data()));
        ::blst_fp_from_bendian(&fp_be, fixed_fp.data());
        ::blst_lendian_from_fp(fp_export.data(), &fp_be);
        const auto fixed_fp_be_expected = fp_oracle.reduce(
            std::as_bytes(std::span(fixed_fp)), ByteOrder::big);
        const auto fixed_fp_be_bytes =
            ModulusOracle<12>::to_le_bytes(fixed_fp_be_expected);
        require(equal_prefix(fp_export, fixed_fp_be_bytes.data()));

        std::array<std::uint32_t, 12> fp_words32{};
        std::array<std::uint64_t, 6> fp_words64{};
        for (std::size_t index = 0; index < fixed_fp.size(); ++index)
        {
            fp_words32[index / 4] |= std::to_integer<std::uint32_t>(
                std::byte{fixed_fp[index]})
                << ((index % 4) * 8);
            fp_words64[index / 8] |= std::to_integer<std::uint64_t>(
                std::byte{fixed_fp[index]})
                << ((index % 8) * 8);
        }
        ::blst_fp_from_uint32(&fp_le, fp_words32.data());
        std::array<std::uint32_t, 12> fp_words32_out{};
        ::blst_uint32_from_fp(fp_words32_out.data(), &fp_le);
        require(fp_words32_out == fixed_fp_le_expected);
        ::blst_fp_from_uint64(&fp_le, fp_words64.data());
        std::array<std::uint64_t, 6> fp_words64_out{};
        ::blst_uint64_from_fp(fp_words64_out.data(), &fp_le);
        for (std::size_t index = 0; index < fp_words64_out.size(); ++index)
        {
            const std::uint64_t expected =
                static_cast<std::uint64_t>(fixed_fp_le_expected[index * 2]) |
                (static_cast<std::uint64_t>(fixed_fp_le_expected[index * 2 + 1]) << 32);
            require(fp_words64_out[index] == expected);
        }

        const blst_fp12 fp12 = native<blst_fp12>(reduce_fp12(input));
        std::array<blst_byte, 48 * 12> fp12_export{};
        std::array<blst_byte, 48 * 12> fp12_repeat{};
        ::blst_bendian_from_fp12(fp12_export.data(), &fp12);
        ::blst_bendian_from_fp12(fp12_repeat.data(), &fp12);
        require(fp12_export == fp12_repeat);

        std::array<blst_byte, 99> hex{};
        const std::size_t hex_size = std::min<std::size_t>(input.size(), 98);
        if (hex_size != 0)
            std::memcpy(hex.data(), input.data(), hex_size);
        const ParsedHex parsed = parse_hex(hex);

        blst_scalar scalar_hex;
        std::memset(&scalar_hex, 0xa5, sizeof(scalar_hex));
        const blst_scalar scalar_hex_before = scalar_hex;
        const bool scalar_hex_ok =
            ::blst_scalar_from_hexascii(&scalar_hex, hex.data());
        const bool scalar_hex_expected = parsed.syntax && parsed.digits <= 64;
        require(scalar_hex_ok == scalar_hex_expected);
        if (scalar_hex_ok)
        {
            ::blst_lendian_from_scalar(scalar_export.data(), &scalar_hex);
            require(std::memcmp(
                scalar_export.data(), parsed.little_endian.data(), 32) == 0);
        }
        else
            require(object_equal(scalar_hex, scalar_hex_before));

        const auto parsed_scalar_span = std::span(parsed.little_endian).first(32);
        const auto parsed_scalar_reduced = scalar_oracle.reduce(
            parsed_scalar_span, ByteOrder::little);
        const auto parsed_scalar_reduced_bytes =
            ModulusOracle<8>::to_le_bytes(parsed_scalar_reduced);
        blst_fr fr_hex;
        std::memset(&fr_hex, 0xa5, sizeof(fr_hex));
        const blst_fr fr_hex_before = fr_hex;
        const bool fr_hex_ok = ::blst_fr_from_hexascii(&fr_hex, hex.data());
        require(fr_hex_ok == scalar_hex_expected);
        if (fr_hex_ok)
        {
            ::blst_scalar_from_fr(&scalar_roundtrip, &fr_hex);
            ::blst_lendian_from_scalar(scalar_export.data(), &scalar_roundtrip);
            require(std::memcmp(
                scalar_export.data(), parsed_scalar_reduced_bytes.data(), 32) == 0);
        }
        else
            require(object_equal(fr_hex, fr_hex_before));

        const auto parsed_fp_reduced = fp_oracle.reduce(
            std::span(parsed.little_endian), ByteOrder::little);
        const auto parsed_fp_reduced_bytes =
            ModulusOracle<12>::to_le_bytes(parsed_fp_reduced);
        blst_fp fp_hex;
        std::memset(&fp_hex, 0xa5, sizeof(fp_hex));
        const blst_fp fp_hex_before = fp_hex;
        const bool fp_hex_ok = ::blst_fp_from_hexascii(&fp_hex, hex.data());
        require(fp_hex_ok == parsed.syntax);
        if (fp_hex_ok)
        {
            ::blst_lendian_from_fp(fp_export.data(), &fp_hex);
            require(std::memcmp(
                fp_export.data(), parsed_fp_reduced_bytes.data(), 48) == 0);
        }
        else
            require(object_equal(fp_hex, fp_hex_before));
    }
}
