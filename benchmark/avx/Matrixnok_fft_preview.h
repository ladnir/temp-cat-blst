// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Category Labs
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>

struct MatrixNoKFftGlv4Schedule
{
    uint8_t sign_mask = 0;
    uint8_t sched_len = 0;
    std::array<uint8_t, 70> schedule = {};
};

extern "C" bool catblst_matrixnok_fft_roots(
    unsigned char* roots_le,
    MatrixNoKFftGlv4Schedule* schedules,
    std::size_t n,
    bool inverse,
    unsigned char scale_le[32],
    MatrixNoKFftGlv4Schedule* scale_schedule);

extern "C" void catblst_matrixnok_g1_generator_coordinates(
    unsigned char coordinates_be[2 * 48]);

extern "C" void catblst_matrixnok_gt_generator_coefficients(
    unsigned char coefficients_be[12 * 48]);
