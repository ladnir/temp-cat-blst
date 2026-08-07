// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs
#pragma once

#include "../common/preprocess.h"
#include "../cpu/vector/vector_impl.hpp"
#include "../generated_matrixnok_constants.hpp"

#include <array>
#include <cstddef>

namespace blst::avx::matrixnok
{
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
        const AVXVector<8> shift(_vv512_const_u64(7, 0, 1, 2, 3, 4, 5, 6));
        BLST_UNROLL
        for (std::size_t diagonal = 0; diagonal < 8; ++diagonal) {
            const AVXVector<8> row(generated::r1_perm_mat_fixed[diagonal]);
            if (diagonal != 0) {
                BLST_UNROLL
                for (std::size_t index = 0; index < batch; ++index) {
                    residues[index] = residues[index].permute(shift);
                }
            }
            BLST_UNROLL
            for (std::size_t index = 0; index < batch; ++index) {
                out_hi[index] = out_hi[index].mulhi(row, residues[index]);
                out_lo[index] = out_lo[index].mullo(row, residues[index]);
            }
        }
    }


    // Register-only latency path. The input is permuted in place.
    INLINE void r1_change_base_raw(
        AVXVector<8>& residues,
        AVXVector<8>& out_hi,
        AVXVector<8>& out_lo)
    {
        const AVXVector<8> shift(_vv512_const_u64(7, 0, 1, 2, 3, 4, 5, 6));
        BLST_UNROLL
        for (std::size_t diagonal = 0; diagonal < 8; ++diagonal) {
            const AVXVector<8> row(generated::r1_perm_mat_fixed[diagonal]);
            if (diagonal != 0) {
                residues = residues.permute(shift);
            }
            out_hi = out_hi.mulhi(row, residues);
            out_lo = out_lo.mullo(row, residues);
        }
    }
} // namespace blst::avx::matrixnok
