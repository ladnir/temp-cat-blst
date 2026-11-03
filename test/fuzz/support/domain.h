// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <array>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <utility>

namespace catblst::fuzz
{
    namespace internal
    {
        struct AdapterAccess;
        struct FactoryAccess;
    }

    /**
     * A fuzz-only value whose tag records a proven library precondition.
     *
     * The storage representation is deliberately not observable here. Only
     * the reviewed factory and adapter layers can construct or unwrap values.
     * This prevents a harness from fabricating a narrow domain with raw bytes.
     */
    template<class Storage, class Tag>
    class DomainValue final
    {
      public:
        using storage_type = Storage;
        using tag_type = Tag;

        DomainValue(const DomainValue&) = default;
        DomainValue(DomainValue&&) noexcept(
            std::is_nothrow_move_constructible_v<Storage>) = default;
        DomainValue& operator=(const DomainValue&) = default;
        DomainValue& operator=(DomainValue&&) noexcept(
            std::is_nothrow_move_assignable_v<Storage>) = default;

      private:
        explicit DomainValue(Storage storage) : storage_(std::move(storage)) {}

        Storage storage_;

        friend struct internal::AdapterAccess;
        friend struct internal::FactoryAccess;
    };

    struct ScalarBitsTag final {};
    struct CanonicalScalarTag final {};
    struct SecretScalarTag final {};
    struct ReducedFrTag final {};
    struct ReducedFpTag final {};
    struct ReducedFp2Tag final {};
    struct ReducedFp12Tag final {};
    template<class InnerTag>
    struct NonzeroTag final {};

    struct FieldValidG1CandidateTag final {};
    struct FieldValidG2CandidateTag final {};
    struct OnCurveG1Tag final {};
    struct OnCurveG1AffineTag final {};
    struct OnCurveG2Tag final {};
    struct OnCurveG2AffineTag final {};
    struct SubgroupG1Tag final {};
    struct SubgroupG1AffineTag final {};
    struct SubgroupG2Tag final {};
    struct SubgroupG2AffineTag final {};
    template<class PointTag>
    struct IncompleteAddPairTag final {};

    struct MillerValueTag final {};
    struct GtElementTag final {};
    struct PreparedFftDomainTag final {};
    struct LiveMiddleProductHandleTag final {};

    enum class PairingOrientation : std::uint8_t
    {
        unset,
        pk_in_g1,
        pk_in_g2,
    };

    enum class PairingPendingState : std::uint8_t
    {
        pending,
        committed,
    };

    template<PairingOrientation Orientation, PairingPendingState Pending>
    struct PairingStateTag final
    {
        static constexpr PairingOrientation orientation = Orientation;
        static constexpr PairingPendingState pending = Pending;
    };

    template<class Storage>
    using ScalarBits = DomainValue<Storage, ScalarBitsTag>;
    template<class Storage>
    using CanonicalScalar = DomainValue<Storage, CanonicalScalarTag>;
    template<class Storage>
    using SecretScalar = DomainValue<Storage, SecretScalarTag>;
    template<class Storage>
    using ReducedFr = DomainValue<Storage, ReducedFrTag>;
    template<class Storage>
    using ReducedFp = DomainValue<Storage, ReducedFpTag>;
    template<class Storage>
    using ReducedFp2 = DomainValue<Storage, ReducedFp2Tag>;
    template<class Storage>
    using ReducedFp12 = DomainValue<Storage, ReducedFp12Tag>;

    template<class Domain>
    using Nonzero = DomainValue<
        typename Domain::storage_type,
        NonzeroTag<typename Domain::tag_type>>;

    template<class Storage>
    using FieldValidG1Candidate = DomainValue<Storage, FieldValidG1CandidateTag>;
    template<class Storage>
    using FieldValidG2Candidate = DomainValue<Storage, FieldValidG2CandidateTag>;
    template<class Storage>
    using OnCurveG1 = DomainValue<Storage, OnCurveG1Tag>;
    template<class Storage>
    using OnCurveG1Affine = DomainValue<Storage, OnCurveG1AffineTag>;
    template<class Storage>
    using OnCurveG2 = DomainValue<Storage, OnCurveG2Tag>;
    template<class Storage>
    using OnCurveG2Affine = DomainValue<Storage, OnCurveG2AffineTag>;
    template<class Storage>
    using SubgroupG1 = DomainValue<Storage, SubgroupG1Tag>;
    template<class Storage>
    using SubgroupG1Affine = DomainValue<Storage, SubgroupG1AffineTag>;
    template<class Storage>
    using SubgroupG2 = DomainValue<Storage, SubgroupG2Tag>;
    template<class Storage>
    using SubgroupG2Affine = DomainValue<Storage, SubgroupG2AffineTag>;

    template<class Domain>
    using IncompleteAddPair = DomainValue<
        std::array<Domain, 2>,
        IncompleteAddPairTag<typename Domain::tag_type>>;

    template<class Storage>
    using MillerValue = DomainValue<Storage, MillerValueTag>;
    template<class Storage>
    using GtElement = DomainValue<Storage, GtElementTag>;
    template<class Storage>
    using PreparedFftDomain = DomainValue<Storage, PreparedFftDomainTag>;
    template<class Storage>
    using LiveMiddleProductHandle =
        DomainValue<Storage, LiveMiddleProductHandleTag>;

    template<
        class Storage,
        PairingOrientation Orientation,
        PairingPendingState Pending>
    using PairingState = DomainValue<
        Storage,
        PairingStateTag<Orientation, Pending>>;

    template<class T, class Tag>
    concept DomainWithTag = requires {
        typename T::tag_type;
        requires std::same_as<typename T::tag_type, Tag>;
    };

    static_assert(!std::is_aggregate_v<DomainValue<std::array<std::byte, 32>, ScalarBitsTag>>);
    static_assert(!std::is_constructible_v<
                  DomainValue<std::array<std::byte, 32>, SecretScalarTag>,
                  std::array<std::byte, 32>>);
}
