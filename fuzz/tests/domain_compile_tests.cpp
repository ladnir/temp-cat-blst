// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#include "fuzz/support/domain.h"

#include <array>
#include <concepts>
#include <cstddef>

namespace
{
    using Bytes32 = std::array<std::byte, 32>;
    using Scalar = catblst::fuzz::ScalarBits<Bytes32>;
    using OnCurve = catblst::fuzz::OnCurveG1<Bytes32>;
    using Subgroup = catblst::fuzz::SubgroupG1<Bytes32>;
    using Fp = catblst::fuzz::ReducedFp<Bytes32>;
    using NonzeroFp = catblst::fuzz::Nonzero<Fp>;

    template<class T>
    concept CanEnterSubgroupMult = catblst::fuzz::DomainWithTag<
        T,
        catblst::fuzz::SubgroupG1Tag>;

    template<class T>
    concept CanEnterUncheckedMult = catblst::fuzz::DomainWithTag<
        T,
        catblst::fuzz::OnCurveG1Tag>;

    template<class T>
    concept CanEnterFpInverse = catblst::fuzz::DomainWithTag<
        T,
        catblst::fuzz::NonzeroTag<catblst::fuzz::ReducedFpTag>>;

    static_assert(CanEnterSubgroupMult<Subgroup>);
    static_assert(!CanEnterSubgroupMult<OnCurve>);
    static_assert(CanEnterUncheckedMult<OnCurve>);
    static_assert(!CanEnterUncheckedMult<Subgroup>);
    static_assert(CanEnterFpInverse<NonzeroFp>);
    static_assert(!CanEnterFpInverse<Fp>);
    static_assert(std::same_as<typename Scalar::tag_type, catblst::fuzz::ScalarBitsTag>);
}
