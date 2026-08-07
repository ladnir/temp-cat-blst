// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs
#pragma once

#include "../common/preprocess.h"
#include "../cpu/vector/vector_impl.hpp"
#include "../generated_matrixnok_constants.hpp"

#include <array>
#include <cstddef>
#include <utility>

namespace blst::avx::matrixnok
{
    inline constexpr std::array<uint64_t, 8> r1_cyclic_shift_indices = {
        7, 0, 1, 2, 3, 4, 5, 6};

    namespace detail
    {
        template<std::size_t diagonal, std::size_t batch, std::size_t... index>
        INLINE void r1_change_base_step(
            std::array<AVXVector<8>, batch>& residues,
            std::array<AVXVector<8>, batch>& out_hi,
            std::array<AVXVector<8>, batch>& out_lo,
            std::index_sequence<index...>)
        {
            static_assert(diagonal < 8);
            const AVXVector<8> row(generated::r1_perm_mat_fixed[diagonal]);
            if constexpr (diagonal != 0) {
                const AVXVector<8> shift(r1_cyclic_shift_indices);
                ((residues[index] = residues[index].permute(shift)), ...);
            }
            ((out_hi[index] = out_hi[index].mulhi(row, residues[index])), ...);
            ((out_lo[index] = out_lo[index].mullo(row, residues[index])), ...);

            if constexpr (diagonal + 1 < 8) {
                r1_change_base_step<diagonal + 1>(
                    residues, out_hi, out_lo, std::index_sequence<index...>{});
            }
        }
    } // namespace detail

    // Accumulate the fixed r1 base-change matrix into independent IFMA high
    // and low accumulators. The MatrixNoK QR representation intentionally
    // omits the conventional quotient-correction term. Callers must satisfy
    // the separately proved input and redundant-range bounds.
    template<std::size_t batch>
    INLINE void r1_change_base_raw_batch(
        std::array<AVXVector<8>, batch>& residues,
        std::array<AVXVector<8>, batch>& out_hi,
        std::array<AVXVector<8>, batch>& out_lo)
    {
        static_assert(batch > 0, "MatrixNoK change of base requires a nonempty fixed batch");
        detail::r1_change_base_step<0>(
            residues, out_hi, out_lo, std::make_index_sequence<batch>{});
    }
} // namespace blst::avx::matrixnok
