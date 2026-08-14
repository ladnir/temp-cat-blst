// Copyright 2026 Category Labs
// SPDX-License-Identifier: Apache-2.0

#include "fuzz/support/call_frame.h"

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
    return 0;
}
