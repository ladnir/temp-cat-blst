#ifndef BLST_CURVE_PAIRING_AVX_H
#define BLST_CURVE_PAIRING_AVX_H

#include <cstddef>

#include "blst/curve/GT.h"
#include "blst/curve/point.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

namespace blst::avx
{
void miller_loop_n(
    GT& ret,
    const G2AffinePoint Q[],
    const G1AffinePoint P[],
    std::size_t n);

void final_exp(GT& ret, const GT& f);

void optimal_ate_pairing(
    GT& ret,
    const G2AffinePoint Q[],
    const G1AffinePoint P[],
    std::size_t n);
} // namespace blst::avx

#endif

#endif // BLST_CURVE_PAIRING_AVX_H
