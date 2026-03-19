#pragma once

#include "blst/backend/backend.h"
#include "blst/avx/fp12/fp12.h"

#include <array>
#include <cstdint>

#if defined(_MSC_VER) && defined(_M_X64)
#include <intrin.h>
#endif

struct VroomGLV4Exp {
    uint8_t sign_mask = 0;
    uint8_t sched_len = 0;
    std::array<uint8_t, 70> schedule = {};
};

namespace vroom_gt_glv4_detail {

static constexpr uint64_t kBlsZ = 0xd201000000010000ULL;

inline uint64_t udiv128_64(uint64_t hi, uint64_t lo, uint64_t div, uint64_t* rem)
{
#if defined(_MSC_VER) && defined(_M_X64)
    unsigned __int64 r = 0;
    const unsigned __int64 q = _udiv128(hi, lo, div, &r);
    *rem = r;
    return q;
#else
    const unsigned __int128 n = (static_cast<unsigned __int128>(hi) << 64) | lo;
    *rem = static_cast<uint64_t>(n % div);
    return static_cast<uint64_t>(n / div);
#endif
}

inline void div_rem_u256_u64(uint64_t q[4], uint64_t& r, const uint64_t n[4], uint64_t d)
{
    uint64_t rem = 0;
    for (int i = 3; i >= 0; --i) {
        q[i] = udiv128_64(rem, n[i], d, &rem);
    }
    r = rem;
}

inline void gt_split_scalar_4d(int64_t k_out[4], const uint64_t k_in[4])
{
    uint64_t temp_k[4] = {k_in[0], k_in[1], k_in[2], k_in[3]};

    for (int i = 0; i < 4; ++i) {
        uint64_t remainder = 0;
        uint64_t quotient[4] = {};

        if ((temp_k[0] | temp_k[1] | temp_k[2] | temp_k[3]) == 0) {
            k_out[i] = 0;
            continue;
        }

        div_rem_u256_u64(quotient, remainder, temp_k, kBlsZ);

        bool carry = false;
        if (remainder > (kBlsZ / 2)) {
            carry = true;
            k_out[i] = static_cast<int64_t>(remainder) - static_cast<int64_t>(kBlsZ);
        } else {
            k_out[i] = static_cast<int64_t>(remainder);
        }

        if (carry) {
            for (int j = 0; j < 4; ++j) {
                ++quotient[j];
                if (quotient[j] != 0) {
                    break;
                }
            }
        }

        for (int j = 0; j < 4; ++j) {
            temp_k[j] = quotient[j];
        }
    }

    if (temp_k[0] != 0) {
        k_out[2] += static_cast<int64_t>(temp_k[0]);
        k_out[0] -= static_cast<int64_t>(temp_k[0]);
    }

    k_out[1] = -k_out[1];
    k_out[3] = -k_out[3];
}

}  // namespace vroom_gt_glv4_detail

inline VroomGLV4Exp vroom_glv4_compile(const uint64_t* scalar)
{
    VroomGLV4Exp props = {};
    int64_t k[4] = {};
    vroom_gt_glv4_detail::gt_split_scalar_4d(k, scalar);

    uint64_t abs_k[4] = {};
    for (int i = 0; i < 4; ++i) {
        if (k[i] < 0) {
            props.sign_mask |= static_cast<uint8_t>(1u << i);
            abs_k[i] = static_cast<uint64_t>(-k[i]);
        } else {
            abs_k[i] = static_cast<uint64_t>(k[i]);
        }
    }

    const uint64_t combined = abs_k[0] | abs_k[1] | abs_k[2] | abs_k[3];
    int msb = 63;
    while (msb >= 0 && ((combined >> msb) & 1u) == 0) {
        --msb;
    }

    if (msb < 0) {
        return props;
    }

    auto push = [&](uint8_t op) {
        if (props.sched_len < props.schedule.size()) {
            props.schedule[props.sched_len++] = op;
        }
    };

    unsigned idx = 0;
    for (int i = 0; i < 4; ++i) {
        if ((abs_k[i] >> msb) & 1u) {
            idx |= (1u << i);
        }
    }
    push(static_cast<uint8_t>(idx));

    int sqr_run = 1;
    for (int bit = msb - 1; bit >= 0; --bit) {
        idx = 0;
        for (int i = 0; i < 4; ++i) {
            if ((abs_k[i] >> bit) & 1u) {
                idx |= (1u << i);
            }
        }

        if (idx == 0) {
            ++sqr_run;
            continue;
        }

        while (sqr_run > 15) {
            push(static_cast<uint8_t>(15u << 4));
            sqr_run -= 15;
        }
        push(static_cast<uint8_t>((sqr_run << 4) | idx));
        sqr_run = 1;
    }

    --sqr_run;
    while (sqr_run > 0) {
        const int chunk = (sqr_run > 15) ? 15 : sqr_run;
        push(static_cast<uint8_t>(chunk << 4));
        sqr_run -= chunk;
    }

    return props;
}

template<class Ring>
inline typename FP12Ring<Ring>::StandardElement vroom_cyclotomic_exp_basic(
    const typename FP12Ring<Ring>::StandardElement& base,
    const uint64_t* scalar,
    const FP12Ring<Ring>& fp12ring,
    const Ring& ring)
{
    using Accumulator = typename FP12Ring<Ring>::template CommonElement<6>;

    Accumulator out = fp12ring.recast_common(fp12ring.one(ring));
    bool started = false;

    for (int limb_idx = 3; limb_idx >= 0; --limb_idx) {
        for (int bit_idx = 63; bit_idx >= 0; --bit_idx) {
            if (started) {
                out = fp12ring.recast_common(fp12ring.fp12_flat_cyclotomic_sqr(out, ring));
            }
            if (((scalar[limb_idx] >> bit_idx) & 1u) == 0) {
                continue;
            }
            if (!started) {
                out = fp12ring.recast_common(base);
                started = true;
            } else {
                out = fp12ring.recast_common(fp12ring.fp12_flat_mul(out, base, ring));
            }
        }
    }

    return fp12ring.fp12_flat_mul(fp12ring.one(ring), out, ring);
}

template<class Ring>
inline typename FP12Ring<Ring>::StandardElement vroom_cyclotomic_exp_glv4(
    const typename FP12Ring<Ring>::StandardElement& base,
    const VroomGLV4Exp& props,
    const FP12Ring<Ring>& fp12ring,
    const Ring& ring)
{
    using Accumulator = typename FP12Ring<Ring>::template CommonElement<6>;
    using StandardElement = typename FP12Ring<Ring>::StandardElement;

    if (props.sched_len == 0) {
        return fp12ring.one(ring);
    }

    std::array<StandardElement, 4> basis = {
        base,
        fp12ring.fp12_flat_frobenius_map1(base, ring),
        fp12ring.fp12_flat_frobenius_map2(base, ring),
        fp12ring.fp12_flat_frobenius_map3(base, ring),
    };

    for (int i = 0; i < 4; ++i) {
        if ((props.sign_mask >> i) & 1u) {
            basis[i] = fp12ring.fp12_flat_conjugate(basis[i], ring);
        }
    }

    std::array<StandardElement, 16> table = {};
    table[0] = fp12ring.one(ring);
    table[1] = basis[0];
    table[2] = basis[1];
    table[4] = basis[2];
    table[8] = basis[3];
    table[3] = fp12ring.fp12_flat_mul(table[1], table[2], ring);
    table[5] = fp12ring.fp12_flat_mul(table[1], table[4], ring);
    table[6] = fp12ring.fp12_flat_mul(table[2], table[4], ring);
    table[7] = fp12ring.fp12_flat_mul(table[3], table[4], ring);
    table[9] = fp12ring.fp12_flat_mul(table[1], table[8], ring);
    table[10] = fp12ring.fp12_flat_mul(table[2], table[8], ring);
    table[11] = fp12ring.fp12_flat_mul(table[3], table[8], ring);
    table[12] = fp12ring.fp12_flat_mul(table[4], table[8], ring);
    table[13] = fp12ring.fp12_flat_mul(table[5], table[8], ring);
    table[14] = fp12ring.fp12_flat_mul(table[6], table[8], ring);
    table[15] = fp12ring.fp12_flat_mul(table[7], table[8], ring);

    Accumulator out = fp12ring.recast_common(table[props.schedule[0] & 0x0f]);
    for (size_t op_idx = 1; op_idx < props.sched_len; ++op_idx) {
        const uint8_t op = props.schedule[op_idx];
        const int sqr_count = op >> 4;
        const int mul_idx = op & 0x0f;

        for (int i = 0; i < sqr_count; ++i) {
            out = fp12ring.recast_common(fp12ring.fp12_flat_cyclotomic_sqr(out, ring));
        }
        if (mul_idx != 0) {
            out = fp12ring.recast_common(fp12ring.fp12_flat_mul(out, table[mul_idx], ring));
        }
    }

    return fp12ring.fp12_flat_mul(fp12ring.one(ring), out, ring);
}

template<class Ring>
inline typename FP12Ring<Ring>::StandardElement vroom_cyclotomic_exp_glv4(
    const typename FP12Ring<Ring>::StandardElement& base,
    const uint64_t* scalar,
    const FP12Ring<Ring>& fp12ring,
    const Ring& ring)
{
    return vroom_cyclotomic_exp_glv4(base, vroom_glv4_compile(scalar), fp12ring, ring);
}
