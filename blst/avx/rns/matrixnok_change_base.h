// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs
//
// Fixed BLS12-381 MatrixNoK construction adapted from the algorithm and
// parameters in Simon Langowski's VROOM commit d6ecd534. This implementation
// is locally composed around CatBLST's selected-vector and proof interfaces.
#pragma once

#include "../common/preprocess.h"
#include "../cpu/vector/vector_impl.hpp"
#include "../generated_matrixnok_constants.hpp"

#include <array>
#include <cstddef>

namespace blst::avx::matrixnok
{
    inline constexpr std::size_t r1_diagonal_count =
        generated::BLS12_381_RNS_DIM;
    static_assert(r1_diagonal_count == 8);
    static_assert(generated::r1_skip_correction);

    namespace detail
    {
        INLINE void debug_assert_r1_input(const AVXVector<8>& residues)
        {
#ifndef NDEBUG
            const auto lanes = residues.to_array();
            for (std::size_t lane = 0; lane < lanes.size(); ++lane) {
                BLST_ASSERT(lanes[lane] < 2 * generated::moduli1[lane]);
            }
#else
            (void)residues;
#endif
        }
    } // namespace detail

    // Add R1(residues) to the caller's IFMA high/low accumulator. Row d holds
    // the d-th cyclic diagonal of the fixed 8x8 matrix. Rotating residues by
    // one lane before each later row therefore reconstructs the ordinary
    // matrix-vector product without horizontal reductions.
    //
    // MatrixNoK's QR representation intentionally omits the conventional
    // approximate-quotient correction. Correctness consequently depends on
    // the generated r1_skip_correction claim and the separately proved input
    // and redundant-range bounds. The input array is scratch and is rotated.
    // This batch schedule is kept separate from the latency-oriented scalar
    // schedule below so each retains its intended register/codegen shape.
    template<std::size_t batch>
    INLINE void r1_change_base_raw_batch(
        std::array<AVXVector<8>, batch>& residues,
        std::array<AVXVector<8>, batch>& out_hi,
        std::array<AVXVector<8>, batch>& out_lo)
    {
        static_assert(batch > 0, "MatrixNoK change of base requires a nonempty fixed batch");
        BLST_ASSERT(&residues != &out_hi);
        BLST_ASSERT(&residues != &out_lo);
        BLST_ASSERT(&out_hi != &out_lo);
#ifndef NDEBUG
        for (std::size_t index = 0; index < batch; ++index) {
            detail::debug_assert_r1_input(residues[index]);
        }
#endif
        const AVXVector<8> shift(_vv512_const_u64(7, 0, 1, 2, 3, 4, 5, 6));
        BLST_UNROLL
        for (std::size_t diagonal = 0; diagonal < r1_diagonal_count; ++diagonal) {
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
    // Register-only latency path with the same recurrence as the batch path.
    // residues is scratch; out_hi/out_lo are in/out accumulators.
    INLINE void r1_change_base_raw(
        AVXVector<8>& residues,
        AVXVector<8>& out_hi,
        AVXVector<8>& out_lo)
    {
        BLST_ASSERT(&residues != &out_hi);
        BLST_ASSERT(&residues != &out_lo);
        BLST_ASSERT(&out_hi != &out_lo);
        detail::debug_assert_r1_input(residues);
        const AVXVector<8> shift(_vv512_const_u64(7, 0, 1, 2, 3, 4, 5, 6));
        BLST_UNROLL
        for (std::size_t diagonal = 0; diagonal < r1_diagonal_count; ++diagonal) {
            const AVXVector<8> row(generated::r1_perm_mat_fixed[diagonal]);
            if (diagonal != 0) {
                residues = residues.permute(shift);
            }
            out_hi = out_hi.mulhi(row, residues);
            out_lo = out_lo.mullo(row, residues);
        }
    }
} // namespace blst::avx::matrixnok
