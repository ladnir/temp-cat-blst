#ifndef BLST_CURVE_PAIRING_AVX_H
#define BLST_CURVE_PAIRING_AVX_H

/*
 * WARNING: This file is intentionally excluded from the active build during the
 * VROOM backend reset. It is retained only as reference and must not be wired
 * back into active code until the post-reset AVX API is defined.
 */

#include <cstddef>
#include <span>

#include "blst/curve/GT.h"
#include "blst/curve/point.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

namespace blst::avx
{
void miller_loop_n(
    GT& ret,
    const G2AffinePoint* BLST_restrict q,
    const G1AffinePoint* BLST_restrict p,
    std::size_t n);

void miller_loop_n(
    GT& ret,
    std::span<const G2AffinePoint> q,
    std::span<const G1AffinePoint> p);

void final_exp(GT& ret, const GT& f);

void optimal_ate_pairing(
    GT& ret,
    const G2AffinePoint* BLST_restrict q,
    const G1AffinePoint* BLST_restrict p,
    std::size_t n);

void optimal_ate_pairing(
    GT& ret,
    std::span<const G2AffinePoint> q,
    std::span<const G1AffinePoint> p);
} // namespace blst::avx

#endif

#endif // BLST_CURVE_PAIRING_AVX_H
