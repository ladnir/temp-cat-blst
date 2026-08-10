// SPDX-License-Identifier: MIT
// Copyright 2026 Simon Langowski
// Copyright 2026 Category Labs
#pragma once

namespace blst::avx {

// Statically bounded product awaiting the selected MatrixNoK basis change.
// Both members own ordinary values; construction allocates and throws nothing.
template<class WideM2, class ReducedM1>
class ReadyToReduce {
    public:
    WideM2 acc;
    ReducedM1 input;
    
    // Constructor for aggregate initialization compatibility
    ReadyToReduce() = default;
    ReadyToReduce(const WideM2 &a, const ReducedM1 &i) : acc(a), input(i) {}
};

}  // namespace blst::avx
