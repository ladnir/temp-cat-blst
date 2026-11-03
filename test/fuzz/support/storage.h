// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "test/fuzz/support/domain.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <new>
#include <vector>

namespace catblst::fuzz
{
    template<std::size_t Size>
    struct OpaqueStorage final
    {
        std::array<std::byte, Size> bytes{};
    };

    using FrStorage = OpaqueStorage<32>;
    using FpStorage = OpaqueStorage<48>;
    using Fp2Storage = OpaqueStorage<96>;
    using Fp12Storage = OpaqueStorage<576>;
    using G1Storage = OpaqueStorage<144>;
    using G1AffineStorage = OpaqueStorage<96>;
    using G2Storage = OpaqueStorage<288>;
    using G2AffineStorage = OpaqueStorage<192>;

    struct PreparedFftStorage final
    {
        std::vector<FrStorage> elements;
        std::size_t n = 0;
    };

    struct MiddleProductHandleStorage final
    {
        std::shared_ptr<void> handle;
        std::size_t n = 0;
    };

    struct AlignedDelete final
    {
        std::size_t alignment = alignof(std::max_align_t);

        void operator()(void* pointer) const noexcept
        {
            if (pointer != nullptr)
                ::operator delete(pointer, std::align_val_t(alignment));
        }
    };

    struct PairingStorage final
    {
        std::unique_ptr<void, AlignedDelete> context;
        std::vector<std::uint8_t> dst;
        std::size_t size = 0;

        PairingStorage(
            void* pointer,
            std::size_t alignment,
            std::vector<std::uint8_t> owned_dst,
            std::size_t context_size) noexcept
            : context(pointer, AlignedDelete{alignment}),
              dst(std::move(owned_dst)),
              size(context_size)
        {}

        PairingStorage(PairingStorage&&) noexcept = default;
        PairingStorage& operator=(PairingStorage&&) noexcept = default;
        PairingStorage(const PairingStorage&) = delete;
        PairingStorage& operator=(const PairingStorage&) = delete;
    };

    using ReducedFrValue = ReducedFr<FrStorage>;
    using ReducedFpValue = ReducedFp<FpStorage>;
    using ReducedFp2Value = ReducedFp2<Fp2Storage>;
    using ReducedFp12Value = ReducedFp12<Fp12Storage>;
    using NonzeroFpValue = Nonzero<ReducedFpValue>;
    using NonzeroFp2Value = Nonzero<ReducedFp2Value>;
    using NonzeroFp12Value = Nonzero<ReducedFp12Value>;

    using FieldValidG1CandidateValue =
        FieldValidG1Candidate<G1AffineStorage>;
    using FieldValidG2CandidateValue =
        FieldValidG2Candidate<G2AffineStorage>;
    using OnCurveG1Value = OnCurveG1<G1Storage>;
    using OnCurveG1AffineValue = OnCurveG1Affine<G1AffineStorage>;
    using OnCurveG2Value = OnCurveG2<G2Storage>;
    using OnCurveG2AffineValue = OnCurveG2Affine<G2AffineStorage>;
    using SubgroupG1Value = SubgroupG1<G1Storage>;
    using SubgroupG1AffineValue = SubgroupG1Affine<G1AffineStorage>;
    using SubgroupG2Value = SubgroupG2<G2Storage>;
    using SubgroupG2AffineValue = SubgroupG2Affine<G2AffineStorage>;
    using IncompleteOnCurveG1Pair = IncompleteAddPair<OnCurveG1Value>;
    using IncompleteOnCurveG2Pair = IncompleteAddPair<OnCurveG2Value>;

    using MillerValueStorage = MillerValue<Fp12Storage>;
    using GtElementValue = GtElement<Fp12Storage>;
    using PreparedFftDomainValue = PreparedFftDomain<PreparedFftStorage>;
    using LiveMiddleProductHandleValue =
        LiveMiddleProductHandle<MiddleProductHandleStorage>;

    template<PairingOrientation Orientation, PairingPendingState Pending>
    using PairingStateValue = PairingState<PairingStorage, Orientation, Pending>;

    using PairingUnsetPending = PairingStateValue<
        PairingOrientation::unset,
        PairingPendingState::pending>;
    using PairingG1Pending = PairingStateValue<
        PairingOrientation::pk_in_g1,
        PairingPendingState::pending>;
    using PairingG2Pending = PairingStateValue<
        PairingOrientation::pk_in_g2,
        PairingPendingState::pending>;
    using PairingUnsetCommitted = PairingStateValue<
        PairingOrientation::unset,
        PairingPendingState::committed>;
    using PairingG1Committed = PairingStateValue<
        PairingOrientation::pk_in_g1,
        PairingPendingState::committed>;
    using PairingG2Committed = PairingStateValue<
        PairingOrientation::pk_in_g2,
        PairingPendingState::committed>;

    static_assert(!std::is_copy_constructible_v<PairingUnsetPending>);
}
