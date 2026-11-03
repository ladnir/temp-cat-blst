// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#include "test/fuzz/support/call_frame.h"

#include <cstddef>
#include <limits>

int main()
{
    using catblst::fuzz::CallFrame;
    using catblst::fuzz::StorageShape;

    std::size_t bytes = 1;
    if (!CallFrame::byte_size(3, 7, bytes) || bytes != 21)
        return 1;
    if (CallFrame::byte_size(
            std::numeric_limits<std::size_t>::max(), 2, bytes) ||
        bytes != 0)
        return 2;

    CallFrame frame = CallFrame::disjoint();
    auto first = frame.add_region(8);
    first[0] = std::byte{0x5a};
    auto second = frame.add_region(16);
    second[0] = std::byte{0xa5};
    if (first[0] != std::byte{0x5a} || second[0] != std::byte{0xa5})
        return 3;
    if (frame.shape() != StorageShape::disjoint ||
        CallFrame::exact_alias().shape() != StorageShape::exact_alias ||
        CallFrame::null_zero_length().shape() != StorageShape::null_zero_length ||
        CallFrame::invalid_size().shape() != StorageShape::invalid_size ||
        CallFrame::short_scratch().shape() != StorageShape::short_scratch ||
        CallFrame::rejected_overlap().shape() != StorageShape::rejected_overlap)
        return 4;

    auto disjoint_frame = CallFrame::disjoint();
    auto disjoint = disjoint_frame.two_ranges(8, 8);
    if (disjoint.first.data() == disjoint.second.data())
        return 5;
    auto alias_frame = CallFrame::exact_alias();
    auto aliased = alias_frame.two_ranges(8, 8);
    auto bad_alias_frame = CallFrame::exact_alias();
    if (aliased.first.data() != aliased.second.data() ||
        !bad_alias_frame.two_ranges(8, 7).first.empty())
        return 6;
    auto null_frame = CallFrame::null_zero_length();
    auto null_zero = null_frame.two_ranges(0, 0);
    if (null_zero.first.data() != nullptr || null_zero.second.data() != nullptr)
        return 7;
    auto short_frame = CallFrame::short_scratch();
    auto short_scratch = short_frame.two_ranges(8, 8);
    if (short_scratch.first.size() != 8 || short_scratch.second.size() != 7)
        return 8;
    auto overlap_frame = CallFrame::rejected_overlap();
    auto overlap = overlap_frame.two_ranges(8, 8);
    if (overlap.second.data() != overlap.first.data() + 7)
        return 9;
    return 0;
}
