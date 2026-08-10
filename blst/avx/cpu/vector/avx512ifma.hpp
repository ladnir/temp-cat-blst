#pragma once

#include "blst/backend/avx512/intrin_compat.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdio>

namespace blst::avx {

inline constexpr int ceil_div(int n, int d) {
    return (n + d - 1) / d;
}

using AVXScalar = vv512i;

// GCC diagnoses its native vector attributes when the vector type appears
// directly as a standard-library template argument. The attributes still
// govern the element representation; contain that template-identity warning
// at this alias rather than suppressing it throughout AVX translation units.
#if defined(__GNUC__) && !defined(__clang__) && defined(BLST_AVX512)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
#endif
template <std::size_t count>
using AVXScalarArray = std::array<AVXScalar, count>;
#if defined(__GNUC__) && !defined(__clang__) && defined(BLST_AVX512)
#pragma GCC diagnostic pop
#endif

static_assert(sizeof(AVXScalar) == 64);
static_assert(sizeof(AVXScalarArray<2>) == 2 * sizeof(AVXScalar));
static_assert(alignof(AVXScalarArray<2>) == alignof(AVXScalar));

inline AVXScalar Scalar(uint64_t v) {
    return _vv512_set1_epi64(static_cast<long long>(v));
}

inline AVXScalar get_scalar_mask(int bits) {
    const uint64_t value = (static_cast<uint64_t>(1) << bits) - 1;
    return Scalar(value);
}

inline AVXScalar scalar_shift(const AVXScalar &s, int bits) {
    return _vv512_srli_epi64(s, bits);
}

template <int limbs>
class AVXVector {
public:
    static constexpr int LIMBS_PER_VEC = 8;
    static constexpr int VEC_LIMBS = ceil_div(limbs, LIMBS_PER_VEC);

    AVXScalarArray<VEC_LIMBS> data;

    inline vv512i operator[] (const int idx) const {
        return data[idx];
    }

    inline vv512i& operator[] (const int idx) {
        return data[idx];
    }

    AVXVector() = default;

    constexpr explicit AVXVector(vv512i v) : data{v} {
        static_assert(VEC_LIMBS == 1, "A literal vector constructor requires exactly one AVX register");
    }

    inline AVXVector(uint64_t v) {
        for (int i = 0; i < VEC_LIMBS; i++) {
            data[i] = Scalar(v);
        }
    }

    inline AVXVector(const std::array<uint64_t, limbs> &array) {
        if constexpr (limbs % LIMBS_PER_VEC == 0) {
            load(array.data());
        } else {
            std::array<uint64_t, VEC_LIMBS * LIMBS_PER_VEC> values{};
            for (int i = 0; i < limbs; i++) {
                values[i] = array[i];
            }
            load(values.data());
        }
    }

    inline AVXVector(const AVXScalarArray<VEC_LIMBS> &array) {
        for (int i = 0; i < VEC_LIMBS; i++) {
            data[i] = array[i];
        }
    }

    inline AVXVector permute(AVXVector perm) const {
        AVXVector out;
        for (int i = 0; i < VEC_LIMBS; i++) {
            out[i] = _vv512_permutexvar_epi64(perm[i], data[i]);
        }
        return out;
    }

    inline vv512i getm512() const {
        return data[0];
    }

    inline AVXVector mulhi(const AVXVector &a, const AVXVector &b) const {
        AVXVector out;
        for (int i = 0; i < VEC_LIMBS; i++) {
            out[i] = _vv512_madd52hi_epu64(data[i], a[i], b[i]);
        }
        return out;
    }

    inline AVXVector mullo(const AVXVector &a, const AVXVector &b) const {
        AVXVector out;
        for (int i = 0; i < VEC_LIMBS; i++) {
            out[i] = _vv512_madd52lo_epu64(data[i], a[i], b[i]);
        }
        return out;
    }

    inline AVXVector mulhi_scalar(const AVXVector &a, const AVXScalar &s) const {
        AVXVector out;
        for (int i = 0; i < VEC_LIMBS; i++) {
            out[i] = _vv512_madd52hi_epu64(data[i], a[i], s);
        }
        return out;
    }

    inline AVXVector mullo_scalar(const AVXVector &a, const AVXScalar &s) const {
        AVXVector out;
        for (int i = 0; i < VEC_LIMBS; i++) {
            out[i] = _vv512_madd52lo_epu64(data[i], a[i], s);
        }
        return out;
    }

    inline AVXVector srli(const int bits) const {
        AVXVector out;
        for (int i = 0; i < VEC_LIMBS; i++) {
            out[i] = _vv512_srli_epi64(data[i], bits);
        }
        return out;
    }

    inline AVXVector slli(const int bits) const {
        AVXVector out;
        for (int i = 0; i < VEC_LIMBS; i++) {
            out[i] = _vv512_slli_epi64(data[i], bits);
        }
        return out;
    }

    inline AVXVector add(const AVXVector &other) const {
        AVXVector out;
        for (int i = 0; i < VEC_LIMBS; i++) {
            out[i] = _vv512_add_epi64(data[i], other[i]);
        }
        return out;
    }

    inline AVXVector sub_min(const AVXVector &other) const {
        AVXVector out;
        for (int i = 0; i < VEC_LIMBS; i++) {
            const vv512i diff = _vv512_sub_epi64(data[i], other[i]);
            out[i] = _vv512_min_epu64(data[i], diff);
        }
        return out;
    }

    inline AVXVector min(const AVXVector &other) const {
        AVXVector out;
        for (int i = 0; i < VEC_LIMBS; i++) {
            out[i] = _vv512_min_epu64(data[i], other[i]);
        }
        return out;
    }

    inline AVXVector add_min(const AVXVector &other) const {
        AVXVector out;
        for (int i = 0; i < VEC_LIMBS; i++) {
            const vv512i sum = _vv512_add_epi64(data[i], other[i]);
            out[i] = _vv512_min_epu64(data[i], sum);
        }
        return out;
    }

    inline AVXVector cmp_sub(const AVXVector &value) const {
        AVXVector out;
        for (int i = 0; i < VEC_LIMBS; i++) {
            const uint64_t mask = _vv512_cmpge_epi64_mask(data[i], value[i]);
            out[i] = _vv512_mask_sub_epi64(data[i], mask, data[i], value[i]);
        }
        return out;
    }

    inline AVXVector cmp_add(const AVXVector &v1, const AVXVector &v2, const AVXVector &val) {
        AVXVector out;
        for (int i = 0; i < VEC_LIMBS; i++) {
            const uint64_t mask = _vv512_cmpgt_epi64_mask(v1[i], v2[i]);
            out[i] = _vv512_mask_add_epi64(data[i], mask, data[i], val[i]);
        }
        return out;
    }

    inline AVXVector mask_add_cond(const AVXVector &v1, const AVXVector &v2, const AVXVector &val) const {
        AVXVector out;
        for (int i = 0; i < VEC_LIMBS; i++) {
            const uint64_t mask = _vv512_cmpgt_epi64_mask(v1[i], v2[i]);
            out[i] = _vv512_mask_add_epi64(data[i], mask, data[i], val[i]);
        }
        return out;
    }

    inline AVXVector and_scalar(const AVXScalar &other) const {
        AVXVector out;
        for (int i = 0; i < VEC_LIMBS; i++) {
            out[i] = _vv512_and_si512(data[i], other);
        }
        return out;
    }

    inline AVXVector bit_or(const AVXVector &other) const {
        AVXVector out;
        for (int i = 0; i < VEC_LIMBS; i++) {
            out[i] = _vv512_or_si512(data[i], other[i]);
        }
        return out;
    }

    inline AVXVector ctime_select(const bool selection, const AVXVector &b) const {
        const uint64_t mask = static_cast<uint64_t>(static_cast<uint8_t>(!selection) * static_cast<uint8_t>(~uint8_t(0)));
        AVXVector out;
        for (int i = 0; i < VEC_LIMBS; i++) {
            out[i] = _vv512_mask_blend_epi64(mask, data[i], b[i]);
        }
        return out;
    }

    inline bool_t is_zero() const {
        const vv512i zero = Scalar(0);
        uint64_t diff = 0;
        for (int i = 0; i < VEC_LIMBS; i++) {
            constexpr uint64_t full_mask = (uint64_t(1) << LIMBS_PER_VEC) - 1;
            const int remaining = limbs - i * LIMBS_PER_VEC;
            const int active = remaining < LIMBS_PER_VEC ? remaining : LIMBS_PER_VEC;
            const uint64_t active_mask = active == LIMBS_PER_VEC ? full_mask : ((uint64_t(1) << active) - 1);
            const uint64_t eq_mask =
                _vv512_cmpge_epi64_mask(data[i], zero) &
                _vv512_cmpge_epi64_mask(zero, data[i]) &
                active_mask;
            diff |= eq_mask ^ active_mask;
        }
        return diff == 0;
    }

    inline bool_t equal(const AVXVector &other) const {
        uint64_t diff = 0;
        for (int i = 0; i < VEC_LIMBS; i++) {
            constexpr uint64_t full_mask = (uint64_t(1) << LIMBS_PER_VEC) - 1;
            const int remaining = limbs - i * LIMBS_PER_VEC;
            const int active = remaining < LIMBS_PER_VEC ? remaining : LIMBS_PER_VEC;
            const uint64_t active_mask = active == LIMBS_PER_VEC ? full_mask : ((uint64_t(1) << active) - 1);
            const uint64_t eq_mask =
                _vv512_cmpge_epi64_mask(data[i], other[i]) &
                _vv512_cmpge_epi64_mask(other[i], data[i]) &
                active_mask;
            diff |= eq_mask ^ active_mask;
        }
        return diff == 0;
    }

    inline AVXVector sub(const AVXVector & other) const {
        AVXVector out;
        for (int i = 0; i < VEC_LIMBS; i++) {
            out[i] = _vv512_sub_epi64(data[i], other[i]);
        }
        return out;
    }

    template<int index>
    inline uint64_t extract() const {
        static_assert(index >= 0 && index < limbs, "AVX lane index is out of range");
        constexpr int vector_index = index / LIMBS_PER_VEC;
        constexpr int lane_index = index % LIMBS_PER_VEC;
        if constexpr (lane_index < 2) {
            return _vv128_extract_epi64<lane_index>(
                _vv512_castsi512_si128(data[vector_index]));
        } else {
            return _vv256_extract_epi64<lane_index % 4>(
                _vv512_extracti64x4_epi64<lane_index / 4>(data[vector_index]));
        }
    }

    inline uint64_t extract0() const {
        return extract<0>();
    }

    inline void store(uint64_t* dst) const {
        for (int i = 0; i < VEC_LIMBS; i++) {
            _vv512_storeu_si512(dst + LIMBS_PER_VEC * i, data[i]);
        }
    }

    inline void load(const uint64_t* src) {
        for (int i = 0; i < VEC_LIMBS; i++) {
            data[i] = _vv512_loadu_si512(src + LIMBS_PER_VEC * i);
        }
    }

    inline void print(const char* label) const {
        std::array<uint64_t, VEC_LIMBS * LIMBS_PER_VEC> values{};
        store(values.data());
        printf("%s: [", label);
        for (int i = 0; i < limbs; i++) {
            printf("%lu ", values[i]);
        }
        printf("]\n");
    }

    inline std::array<uint64_t, limbs> to_array() const {
        std::array<uint64_t, VEC_LIMBS * LIMBS_PER_VEC> values{};
        store(values.data());
        std::array<uint64_t, limbs> result{};
        std::copy(values.data(), values.data() + limbs, result.data());
        return result;
    }
};

class ShortVector {
    vv128i data;

public:
    ShortVector() = default;
    constexpr explicit ShortVector(vv128i value) : data(value) {}

    template <int limbs>
    static inline ShortVector from_avx_lo(const AVXVector<limbs>& value) {
        return ShortVector(_vv512_castsi512_si128(value[0]));
    }

    inline ShortVector add(const ShortVector& other) const {
        return ShortVector(_vv128_add_epi64(data, other.data));
    }

    inline ShortVector srli(int bits) const {
        return ShortVector(_vv128_srli_epi64(data, bits));
    }

    template <int limbs>
    inline AVXVector<limbs> broadcast_to_avx() const {
        AVXVector<limbs> out;
        const vv512i broadcast = _vv512_broadcastq_epi64(data);
        for (int i = 0; i < AVXVector<limbs>::VEC_LIMBS; ++i) {
            out[i] = broadcast;
        }
        return out;
    }
};

}  // namespace blst::avx
