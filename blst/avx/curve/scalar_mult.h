/*
 * Copyright Supranational LLC and CatBLST contributors
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 *
 * Imported with the VROOM AVX backend and adapted to CatBLST's scalar,
 * constant-time, and selected-vector contracts.
 */
#pragma once
#include "ec.h"
#include "../common/constexpr_utils.h"
#include "../fields/fp2.h"
#include "blst/backend/backend.h"
#include "blst/common/bytes.h"
#include "blst/common/vect.h"
#include "blst/scalar/curve/ec_mult.h"
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <utility>

#if defined(_MSC_VER)
# define BLST_AVX_SCALAR_MULT_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
# define BLST_AVX_SCALAR_MULT_NOINLINE __attribute__((noinline))
#else
# define BLST_AVX_SCALAR_MULT_NOINLINE
#endif

namespace blst::avx {

template<int table_log_size>
inline constexpr bool valid_precomputed_table_log_size =
    table_log_size > 0 && table_log_size < 8;

template<class Curve, class Ring, int table_log_size>
class PrecomputedTable {
    static_assert(valid_precomputed_table_log_size<table_log_size>,
                  "PrecomputedTable table_log_size must be between 1 and 7");

    public:
    using ProjectivePoint = typename Curve::ProjPoint;

    std::array<ProjectivePoint, (1 << table_log_size)> data;

    PrecomputedTable() = default;

    PrecomputedTable(const std::array<ProjectivePoint, 1 << table_log_size> &data) : data(data) {}

    PrecomputedTable(const ProjectivePoint &P, const Curve &curve, const Ring &ring) {
        build(P, curve, ring);
    }

    void build(const ProjectivePoint &P, const Curve &curve, const Ring &ring) {
        data[0] = P;
        // Build [P, 2P, ..., 2^table_log_size P] with one group operation
        // per new entry: double the half-indexed value or add P.
        for (int i = 1; i < (1 << table_log_size); i++) {
            if (i % 2 == 1) {
                data[i] = curve.double_point(data[i / 2], ring);
            } else {
                data[i] = curve.add_point(data[i - 1], P, ring);
            }
        }
    }

    std::pair<uint8_t, bool> booth_encode(uint64_t idx) const {
        BLST_ASSERT(idx <= (uint64_t{1} << (table_log_size + 1)));
        const uint64_t booth_sign = static_cast<uint64_t>((idx >> table_log_size) != 0);
        const uint64_t sign_mask = uint64_t{0} - booth_sign;
        const uint64_t negated = (uint64_t{1} << (table_log_size + 1)) - idx;
        const uint64_t booth_index = (idx & ~sign_mask) | (negated & sign_mask);
        return std::make_pair(static_cast<uint8_t>(booth_index), booth_sign != 0);
    }

    ProjectivePoint gather_booth(uint8_t booth_index, bool booth_sign, const Curve &curve, const Ring &ring) const {
        BLST_ASSERT(booth_index <= (1U << table_log_size));
        ProjectivePoint result = curve.zero(ring);
        for (uint8_t i = 1; i <= (1 << table_log_size); i++) {
            bool match = (i == booth_index);
            // Replace result with data[i-1] if we match
            result = result.ctime_replace_if(match, data[i - 1]);
        }
        // Negate the selected row for a negative Booth digit.
        result.y = result.y.ctime_replace_if(booth_sign, ring.standard_negate(result.y));
        return result;
    }

    ProjectivePoint gather_booth_vartime(uint8_t booth_index, bool booth_sign, const Curve &curve, const Ring &ring) const {
        BLST_ASSERT(booth_index <= (1U << table_log_size));
        if (booth_index == 0) {
            return curve.zero(ring);
        }
        ProjectivePoint result = data[booth_index - 1];
        if (booth_sign) {
            result.y = ring.standard_negate(result.y);
        }
        return result;
    }
};

template<int scalar_bits>
class ScalarForExponent {
    public:
    static_assert(scalar_bits > 0, "ScalarForExponent requires a positive bit width");
    static constexpr int limbs = static_cast<int>(blst::ceil_div(
        static_cast<std::size_t>(scalar_bits), std::size_t{8}));
    std::array<uint8_t, limbs> data;

    ScalarForExponent() = default;

    ScalarForExponent(const std::array<uint8_t, limbs> &raw_data) : data(raw_data) {}

    template<int num_bits>
    uint8_t get_bits(int start_bit) const  {
        static_assert(num_bits <= 8, "num_bits must be at most 8");
        
        // Check if request goes beyond scalar bits - return zero if so
        if (start_bit >= scalar_bits) {
            return 0;
        }
        
        // Clamp num_bits to not exceed scalar_bits
        int actual_num_bits = num_bits;
        if (start_bit + actual_num_bits > scalar_bits) {
            actual_num_bits = scalar_bits - start_bit;
        }
        
        int start_byte = start_bit / 8;
        int start_bit_in_byte = start_bit % 8;
        int bits_available_in_first_byte = 8 - start_bit_in_byte;
        
        // Extract bits from first byte
        uint8_t result = 0;
        if (start_byte < limbs) {
            if (actual_num_bits <= bits_available_in_first_byte) {
                // All bits fit in first byte
                uint8_t mask = (1 << actual_num_bits) - 1;
                result = (data[start_byte] >> start_bit_in_byte) & mask;
            } else {
                // Bits cross byte boundary
                uint8_t first_mask = (1 << bits_available_in_first_byte) - 1;
                uint8_t first_bits = (data[start_byte] >> start_bit_in_byte) & first_mask;
                result = first_bits;
                
                // Extract remaining bits from next byte
                int next_byte = start_byte + 1;
                int bits_needed_from_next = actual_num_bits - bits_available_in_first_byte;
                if (next_byte < limbs) {
                    uint8_t next_mask = (1 << bits_needed_from_next) - 1;
                    uint8_t next_bits = data[next_byte] & next_mask;
                    result |= (next_bits << bits_available_in_first_byte);
                }
                // If next_byte >= limbs, the remaining bits are zero (already set)
            }
        }
        // If start_byte >= limbs, result is already zero
        
        return result;
    }
};

template<class Curve, class Ring, int scalar_bits, int table_log_size, int num_tables>
typename Curve::ProjPoint table_msm(const std::array<ScalarForExponent<scalar_bits>, num_tables> &scalars, const std::array<PrecomputedTable<Curve, Ring, table_log_size>, num_tables> &tables, const Curve &curve, const Ring &ring) {
    static_assert(scalar_bits > 0);
    static_assert(num_tables > 0);
    constexpr int wbits = table_log_size;
    constexpr int num_windows = (scalar_bits + wbits) / (wbits + 1);
    constexpr int carry_possible = scalar_bits % (wbits + 1) == 0;
    std::array<std::array<std::pair<uint8_t, bool>, num_windows>, num_tables> booth_indexes;
    typename Curve::ProjPoint result = curve.zero(ring);
    
    for (int i = 0; i < num_tables; i++) {
        bool carry = false;
        for (int j = 0; j < num_windows; j++) {
            uint64_t window_val = scalars[i].template get_bits<wbits + 1>(j * (wbits + 1));
            window_val += carry ? 1 : 0;
            booth_indexes[i][j] = tables[i].booth_encode(window_val);
            carry = booth_indexes[i][j].second;
        }
        // A full final window may emit one additional carry digit.
        if constexpr (carry_possible) {
            typename Curve::ProjPoint carry_point = curve.zero(ring).ctime_replace_if(carry, tables[i].data[(1 << wbits) - 1]);
            result = curve.add_point(result, carry_point, ring);
        }
        else {
            BLST_ASSERT(!carry);
        }
    }
    // Align the extra carry digit with the fixed-width windows.
    if constexpr (carry_possible) {
        result = curve.double_point(result, ring);
    }
    for (int j = num_windows - 1; j >= 0; j--) {
        for (int i = 0; i < num_tables; i++) {
            result = curve.add_point(result, tables[i].gather_booth(booth_indexes[i][j].first, booth_indexes[i][j].second, curve, ring), ring);
        }
        if (j != 0) {
            for (int k = 0; k < wbits+1; k++) {
                result = curve.double_point(result, ring);
            }
        }
    }
    blst::vec_zero(booth_indexes.data(), sizeof(booth_indexes));
    return result;
}

template<class Curve, class Ring, int scalar_bits, int table_log_size, int max_tables>
BLST_AVX_SCALAR_MULT_NOINLINE typename Curve::ProjPoint table_msm_bounded(
    const std::array<ScalarForExponent<scalar_bits>, max_tables>& scalars,
    const std::array<PrecomputedTable<Curve, Ring, table_log_size>, max_tables>& tables,
    size_t num_tables,
    const Curve& curve,
    const Ring& ring)
{
    static_assert(scalar_bits > 0);
    static_assert(max_tables > 0);
    constexpr int wbits = table_log_size;
    constexpr int num_windows = (scalar_bits + wbits) / (wbits + 1);
    constexpr int carry_possible = scalar_bits % (wbits + 1) == 0;
    BLST_ASSERT(num_tables <= static_cast<size_t>(max_tables));
    std::array<std::array<std::pair<uint8_t, bool>, num_windows>, max_tables>
        booth_indexes;
    typename Curve::ProjPoint result = curve.zero(ring);

    for (size_t i = 0; i < num_tables; ++i) {
        bool carry = false;
        for (int j = 0; j < num_windows; ++j) {
            uint64_t window_val =
                scalars[i].template get_bits<wbits + 1>(j * (wbits + 1));
            window_val += carry ? 1 : 0;
            booth_indexes[i][j] = tables[i].booth_encode(window_val);
            carry = booth_indexes[i][j].second;
        }
        if constexpr (carry_possible) {
            const auto carry_point = curve.zero(ring).ctime_replace_if(
                carry, tables[i].data[(1 << wbits) - 1]);
            result = curve.add_point(result, carry_point, ring);
        }
        else {
            BLST_ASSERT(!carry);
        }
    }
    if constexpr (carry_possible) {
        result = curve.double_point(result, ring);
    }
    for (int j = num_windows - 1; j >= 0; --j) {
        for (size_t i = 0; i < num_tables; ++i) {
            result = curve.add_point(
                result,
                tables[i].gather_booth(
                    booth_indexes[i][j].first,
                    booth_indexes[i][j].second,
                    curve,
                    ring),
                ring);
        }
        if (j != 0) {
            for (int k = 0; k < wbits + 1; ++k) {
                result = curve.double_point(result, ring);
            }
        }
    }
    blst::vec_zero(booth_indexes.data(), sizeof(booth_indexes));
    return result;
}

template<class Curve, class Ring, int scalar_bits, int table_log_size, int num_tables>
typename Curve::ProjPoint table_msm_vartime(const std::array<ScalarForExponent<scalar_bits>, num_tables> &scalars, const std::array<PrecomputedTable<Curve, Ring, table_log_size>, num_tables> &tables, const Curve &curve, const Ring &ring) {
    static_assert(scalar_bits > 0);
    static_assert(num_tables > 0);
    constexpr int wbits = table_log_size;
    constexpr int num_windows = (scalar_bits + wbits) / (wbits + 1);
    constexpr int carry_possible = scalar_bits % (wbits + 1) == 0;
    std::array<std::array<std::pair<uint8_t, bool>, num_windows>, num_tables> booth_indexes;
    typename Curve::ProjPoint result = curve.zero(ring);

    for (int i = 0; i < num_tables; ++i) {
        bool carry = false;
        for (int j = 0; j < num_windows; ++j) {
            uint64_t window_val = scalars[i].template get_bits<wbits + 1>(j * (wbits + 1));
            window_val += static_cast<uint64_t>(carry);
            booth_indexes[i][j] = tables[i].booth_encode(window_val);
            carry = booth_indexes[i][j].second;
        }
        if constexpr (carry_possible) {
            if (carry) {
                result = curve.add_point(result, tables[i].data[(1 << wbits) - 1], ring);
            }
        }
        else {
            BLST_ASSERT(!carry);
        }
    }
    if constexpr (carry_possible) {
        result = curve.double_point(result, ring);
    }
    for (int j = num_windows - 1; j >= 0; --j) {
        for (int i = 0; i < num_tables; ++i) {
            const auto [index, sign] = booth_indexes[i][j];
            if (index != 0) {
                result = curve.add_point(
                    result,
                    tables[i].gather_booth_vartime(index, sign, curve, ring),
                    ring);
            }
        }
        if (j != 0) {
            for (int k = 0; k < wbits + 1; ++k) {
                result = curve.double_point(result, ring);
            }
        }
    }
    return result;
}




// Wrapper functions for blst division operations
// Returns {remainder, quotient}: {sk_zz_r, sk_zz}.
using div_by_zz_result =
    std::pair<std::array<uint8_t, 16>, std::array<uint8_t, 16>>;

template<bool clear_scalar = false>
inline void div_by_zz_into(div_by_zz_result& out, const std::array<uint8_t, 32>& a) {
    union { blst::vec_256 l; blst::pow256 s; } val;
    blst::limbs_from_le_bytes(val.l, a.data(), 32);
    blst::div_by_zz(val.l);
    blst::le_bytes_from_limbs(val.s, val.l, 32);

    // blst_div_by_zz stores: val.s is remainder, val.s+16 is quotient.
    // Preserve the established result ordering: {sk_zz_r, sk_zz}.
    std::memcpy(out.first.data(), val.s, 16);
    std::memcpy(out.second.data(), val.s + 16, 16);

    if constexpr (clear_scalar)
        blst::vec_zero(&val, sizeof(val));
}

inline div_by_zz_result div_by_zz(const std::array<uint8_t, 32>& a) {
    div_by_zz_result out;
    div_by_zz_into(out, a);
    return out;
}

// Returns {remainder, quotient}.
using div_by_z_result =
    std::pair<std::array<uint8_t, 8>, std::array<uint8_t, 8>>;

template<bool clear_scalar = false>
inline void div_by_z_into(div_by_z_result& out, const std::array<uint8_t, 16>& a) {
    typedef blst::byte pow128[16];
    typedef blst::limb_t vec128[blst::NLIMBS(128)];

    union { vec128 l; pow128 s; } val;
    blst::limbs_from_le_bytes(val.l, a.data(), 16);
    blst::div_by_z(val.l);
    blst::le_bytes_from_limbs(val.s, val.l, 16);

    // blst_div_by_z stores: val.s is remainder, val.s+8 is quotient.
    // Preserve the established result ordering.
    std::memcpy(out.first.data(), val.s, 8);
    std::memcpy(out.second.data(), val.s + 8, 8);

    if constexpr (clear_scalar)
        blst::vec_zero(&val, sizeof(val));
}

inline div_by_z_result div_by_z(const std::array<uint8_t, 16>& a) {
    div_by_z_result out;
    div_by_z_into(out, a);
    return out;
}

template<class Ring, class Curve>
class GLVMult {
    public:
    using G1Curve = Curve;
    using ProjectivePoint = typename G1Curve::ProjPoint;

    typename Ring::StandardElement beta;
    G1Curve curve;

    constexpr GLVMult(
        const G1Curve& curve,
        const typename Ring::StandardElement& beta_value) noexcept
        : beta(beta_value), curve(curve)
    {}

    template<int batch_size>
    std::array<ProjectivePoint, batch_size> batch_endomorphism(const std::array<ProjectivePoint, batch_size> &P, const Ring &ring) const {
        static_assert(batch_size > 0);
        std::array<ProjectivePoint, batch_size> result;
        for (int i = 0; i < batch_size; i++) {
            result[i].x = ring.modmul(P[i].x, beta);
            result[i].y = ring.standard_negate(P[i].y);
            result[i].z = P[i].z;
        }
        return result;
    }

    template<int table_log_size=4>
    ProjectivePoint multiply(const ProjectivePoint &P, const std::array<uint8_t, 32> &SK, const Ring &ring) const {
        constexpr int table_size = 1 << table_log_size;
        div_by_zz_result split;
        div_by_zz_into<true>(split, SK);
        std::array<ScalarForExponent<128>, 2> scalars = {
            split.first,
            split.second,
        };
        blst::vec_zero(&split, sizeof(split));
        PrecomputedTable<G1Curve, Ring, table_log_size> tableP(P, curve, ring);
        // Keep independently owned projective tables. A shared-z representation
        // requires separate code-generation and benchmark evidence.
        auto table_phi_data = batch_endomorphism<table_size>(tableP.data, ring);
        PrecomputedTable<G1Curve, Ring, table_log_size> tablePhi(table_phi_data);

        std::array<PrecomputedTable<G1Curve, Ring, table_log_size>, 2> tables = {tableP, tablePhi};
        ProjectivePoint result = table_msm<G1Curve, Ring, 128, table_log_size, 2>(scalars, tables, curve, ring);
        blst::vec_zero(&scalars, sizeof(scalars));
        return result;
    }

    /**
     * Constant-time MSM for one fixed-size batch of subgroup points.
     *
     * The batch shares one doubling schedule across every GLV component. The
     * point tables are public, while the split scalars and Booth indices are
     * secret and are cleared before return. `SK` contains `batch_size`
     * consecutive canonical 32-byte little-endian scalars.
     *
     * `scalars` and `tables` are caller-owned reusable storage. The function
     * allocates nothing. Keep `batch_size` bounded: increasing it expands the
     * live table set and can increase cache traffic enough to erase the saved
     * doublings.
     */
    template<int batch_size>
    using BatchScalars =
        std::array<ScalarForExponent<128>, 2 * batch_size>;

    template<int batch_size, int table_log_size=4>
    using BatchTables = std::array<
        PrecomputedTable<G1Curve, Ring, table_log_size>, 2 * batch_size>;

    template<int batch_size>
    using BatchWnafDigits =
        std::array<std::array<int8_t, 129>, 2 * batch_size>;

    template<int batch_size, int width=5>
    using BatchWnafTables = std::array<
        std::array<ProjectivePoint, size_t{1} << (width - 2)>,
        2 * batch_size>;

    template<int batch_size, int table_log_size=4>
    ProjectivePoint multiply_batch(
        const ProjectivePoint* P,
        const uint8_t* SK,
        const Ring& ring,
        BatchScalars<batch_size>& scalars,
        BatchTables<batch_size, table_log_size>& tables) const
    {
        static_assert(batch_size > 1);
        static_assert(table_log_size == 4,
            "the selected AVX G1 constant-time batch kernel uses width 5");
        constexpr int table_size = 1 << table_log_size;
        constexpr int component_count = 2 * batch_size;

        for (int i = 0; i < batch_size; ++i) {
            std::array<uint8_t, 32> scalar_bytes;
            std::memcpy(scalar_bytes.data(), SK + 32 * i, scalar_bytes.size());
            div_by_zz_result split;
            div_by_zz_into<true>(split, scalar_bytes);
            blst::vec_zero(scalar_bytes.data(), scalar_bytes.size());

            scalars[2 * i] = ScalarForExponent<128>(split.first);
            scalars[2 * i + 1] = ScalarForExponent<128>(split.second);
            blst::vec_zero(&split, sizeof(split));

            tables[2 * i].build(P[i], curve, ring);
            for (int j = 0; j < table_size; ++j) {
                const ProjectivePoint& point = tables[2 * i].data[j];
                tables[2 * i + 1].data[j] = {
                    ring.modmul(point.x, beta),
                    ring.standard_negate(point.y),
                    point.z,
                };
            }
        }

        ProjectivePoint result = table_msm_bounded<
            G1Curve, Ring, 128, table_log_size, component_count>(
                scalars, tables, component_count, curve, ring);
        blst::vec_zero(&scalars, sizeof(scalars));
        return result;
    }

    /**
     * Variable-time fixed-batch MSM for public scalars.
     *
     * This uses the same GLV split, point tables, and shared doubling schedule
     * as `multiply_batch`, but directly indexes public Booth digits and skips
     * zero digits. `SK` contains `batch_size` consecutive canonical 32-byte
     * little-endian scalars. `scalars` and `tables` are caller-owned reusable
     * storage, and the operation allocates nothing.
     */
    template<int batch_size, int width=5>
    ProjectivePoint multiply_batch_vartime(
        const ProjectivePoint* P,
        const uint8_t* SK,
        const Ring& ring,
        BatchWnafDigits<batch_size>& digits,
        BatchWnafTables<batch_size, width>& tables) const
    {
        static_assert(batch_size > 1);
        static_assert(width >= 4 && width <= 6);
        constexpr size_t table_size = size_t{1} << (width - 2);
        constexpr int component_count = 2 * batch_size;
        size_t used = 0;

        for (int i = 0; i < batch_size; ++i) {
            std::array<uint8_t, 32> scalar_bytes;
            std::memcpy(scalar_bytes.data(), SK + 32 * i, scalar_bytes.size());
            const div_by_zz_result split = div_by_zz(scalar_bytes);
            const size_t used_low =
                blst::scalar::recode_wnaf_vartime<128, width>(
                    digits[2 * i].data(), split.first.data());
            const size_t used_high =
                blst::scalar::recode_wnaf_vartime<128, width>(
                    digits[2 * i + 1].data(), split.second.data());
            used = std::max(used, std::max(used_low, used_high));

            const ProjectivePoint twice = curve.double_point(P[i], ring);
            tables[2 * i][0] = P[i];
            for (size_t j = 1; j < table_size; ++j) {
                tables[2 * i][j] =
                    curve.add_point(tables[2 * i][j - 1], twice, ring);
            }
            for (size_t j = 0; j < table_size; ++j) {
                const ProjectivePoint& point = tables[2 * i][j];
                tables[2 * i + 1][j] = {
                    ring.modmul(point.x, beta),
                    ring.standard_negate(point.y),
                    point.z,
                };
            }
        }

        ProjectivePoint result = curve.zero(ring);
        bool started = false;
        for (size_t position = used; position != 0;) {
            --position;
            if (started) {
                result = curve.double_point(result, ring);
            }
            for (int component = 0; component < component_count; ++component) {
                const int digit = digits[component][position];
                if (digit == 0) {
                    continue;
                }
                const unsigned magnitude = digit < 0
                    ? static_cast<unsigned>(-digit)
                    : static_cast<unsigned>(digit);
                BLST_ASSERT((magnitude & 1U) != 0);
                BLST_ASSERT((magnitude >> 1) < table_size);
                ProjectivePoint row = tables[component][magnitude >> 1];
                if (digit < 0) {
                    row.y = ring.standard_negate(row.y);
                }
                if (!started) {
                    result = row;
                    started = true;
                }
                else {
                    result = curve.add_point(result, row, ring);
                }
            }
        }
        return result;
    }

    template<int table_log_size=4>
    ProjectivePoint multiply_vartime(const ProjectivePoint &P, const std::array<uint8_t, 32> &SK, const Ring &ring) const {
        static_assert(table_log_size == 4,
            "the selected AVX G1 variable-time kernel uses wNAF-5");
        return multiply_wnaf_vartime<5>(P, SK, ring);
    }

    template<int width>
    ProjectivePoint multiply_wnaf_vartime(
        const ProjectivePoint& P, const std::array<uint8_t, 32>& SK,
        const Ring& ring) const
    {
        static_assert(width >= 4 && width <= 6);
        constexpr size_t table_size = size_t{1} << (width - 2);
        auto [low, high] = div_by_zz(SK);
        int8_t digits[2][129] = {};
        const size_t used0 = blst::scalar::recode_wnaf_vartime<128, width>(
            digits[0], low.data());
        const size_t used1 = blst::scalar::recode_wnaf_vartime<128, width>(
            digits[1], high.data());
        const size_t used = used0 > used1 ? used0 : used1;
        BLST_ASSERT(used <= 129);

        std::array<ProjectivePoint, table_size> table_p = {};
        std::array<ProjectivePoint, table_size> table_phi = {};
        const ProjectivePoint twice = curve.double_point(P, ring);
        table_p[0] = P;
        for (size_t i = 1; i < table_size; ++i)
            table_p[i] = curve.add_point(table_p[i - 1], twice, ring);
        table_phi = batch_endomorphism<table_size>(table_p, ring);

        ProjectivePoint result = curve.zero(ring);
        bool started = false;
        for (size_t position = used; position != 0;) {
            --position;
            if (started)
                result = curve.double_point(result, ring);
            for (size_t component = 0; component < 2; ++component) {
                const int digit = digits[component][position];
                if (digit == 0)
                    continue;
                const unsigned magnitude = digit < 0
                    ? static_cast<unsigned>(-digit)
                    : static_cast<unsigned>(digit);
                BLST_ASSERT((magnitude & 1U) != 0);
                BLST_ASSERT((magnitude >> 1) < table_size);
                ProjectivePoint row = component == 0
                    ? table_p[magnitude >> 1]
                    : table_phi[magnitude >> 1];
                if (digit < 0)
                    row.y = ring.standard_negate(row.y);
                if (!started) {
                    result = row;
                    started = true;
                }
                else {
                    result = curve.add_point(result, row, ring);
                }
            }
        }
        blst::vec_zero(digits, sizeof(digits));
        return result;
    }

};

template<class Ring, class Curve>
class GLSMult {
    public:
    using G2Curve = Curve;
    using FP2RingType = FP2Ring<Ring>;
    using ProjectivePoint = typename G2Curve::ProjPoint;

    typename Ring::StandardElement psi_x_i;
    typename Ring::StandardElement psi_y_r;
    typename Ring::StandardElement psi_y_i;

    constexpr GLSMult(
        const typename Ring::StandardElement& psi_x_i_value,
        const typename Ring::StandardElement& psi_y_r_value,
        const typename Ring::StandardElement& psi_y_i_value) noexcept
        : psi_x_i(psi_x_i_value),
          psi_y_r(psi_y_r_value),
          psi_y_i(psi_y_i_value)
    {}
   
    ProjectivePoint minus_psi(const ProjectivePoint &P, const Ring &ring) const {
        // Unroll the Fp2 products, skip multiplication by zero, and combine the
        // final negation. For this coefficient, psi_y_i == -psi_y_r (mod p).
        auto x_r = P.x.y() * psi_x_i;
        auto x_i = P.x.x() * psi_x_i;
        auto common = P.y.y() * psi_y_r;
        auto y_r = common + P.y.x() * psi_y_i;
        auto y_i = common + P.y.x() * psi_y_r;
        auto z_r = P.z.x();
        auto z_i = ring.standard_negate(P.z.y());
        auto [x_r_reduced, x_i_reduced, y_r_reduced, y_i_reduced] = ring.batch_reduce_expand(x_r, x_i, y_r, y_i);
        return ProjectivePoint(typename FP2RingType::StandardElement(x_r_reduced, x_i_reduced), typename FP2RingType::StandardElement(y_r_reduced, y_i_reduced), typename FP2RingType::StandardElement(z_r, z_i));
    }

    template<int table_log_size=4>
    ProjectivePoint multiply(const ProjectivePoint &P, const std::array<uint8_t, 32> &SK, const FP2RingType &fp2ring, const G2Curve &curve, const Ring &ring) const {
        constexpr int table_size = 1 << table_log_size;
        using Table = PrecomputedTable<G2Curve, FP2RingType, table_log_size>;
        div_by_zz_result split_zz;
        div_by_z_result split_low;
        div_by_z_result split_high;
        div_by_zz_into<true>(split_zz, SK);
        div_by_z_into<true>(split_low, split_zz.first);
        div_by_z_into<true>(split_high, split_zz.second);
        std::array<ScalarForExponent<64>, 4> scalars = {
            split_low.first,
            split_low.second,
            split_high.first,
            split_high.second,
        };
        blst::vec_zero(&split_zz, sizeof(split_zz));
        blst::vec_zero(&split_low, sizeof(split_low));
        blst::vec_zero(&split_high, sizeof(split_high));
        std::array<Table, 4> tables;
        tables[0].build(P, curve, fp2ring);
        for (int i = 0; i < table_size; i++) {
            tables[1].data[i] = minus_psi(tables[0].data[i], ring);
            tables[2].data[i] = minus_psi(tables[1].data[i], ring);
            tables[3].data[i] = minus_psi(tables[2].data[i], ring);
        }
        ProjectivePoint result = table_msm<G2Curve, FP2RingType, 64, table_log_size, 4>(scalars, tables, curve, fp2ring);
        blst::vec_zero(&scalars, sizeof(scalars));
        return result;
    }

    template<int table_log_size=4>
    ProjectivePoint multiply_vartime(const ProjectivePoint &P, const std::array<uint8_t, 32> &SK, const FP2RingType &fp2ring, const G2Curve &curve, const Ring &ring) const {
        static_assert(table_log_size == 4,
            "the selected AVX G2 variable-time entry uses wNAF-4");
        return multiply_wnaf_vartime<4>(P, SK, fp2ring, curve, ring);
    }

    template<int width>
    ProjectivePoint multiply_wnaf_vartime(
        const ProjectivePoint& P, const std::array<uint8_t, 32>& SK,
        const FP2RingType& fp2ring, const G2Curve& curve,
        const Ring& ring) const
    {
        static_assert(width >= 4 && width <= 6);
        constexpr size_t table_size = size_t{1} << (width - 2);
        auto [low128, high128] = div_by_zz(SK);
        auto [d0, d1] = div_by_z(low128);
        auto [d2, d3] = div_by_z(high128);
        const std::array<const uint8_t*, 4> scalar_digits = {
            d0.data(), d1.data(), d2.data(), d3.data()};
        int8_t digits[4][65] = {};
        size_t used = 0;
        for (size_t i = 0; i < 4; ++i) {
            const size_t count = blst::scalar::recode_wnaf_vartime<64, width>(
                digits[i], scalar_digits[i]);
            BLST_ASSERT(count <= 65);
            if (count > used)
                used = count;
        }

        std::array<std::array<ProjectivePoint, table_size>, 4> table = {};
        const ProjectivePoint twice = curve.double_point(P, fp2ring);
        table[0][0] = P;
        for (size_t i = 1; i < table_size; ++i)
            table[0][i] = curve.add_point(table[0][i - 1], twice, fp2ring);
        for (size_t i = 0; i < table_size; ++i) {
            table[1][i] = minus_psi(table[0][i], ring);
            table[2][i] = minus_psi(table[1][i], ring);
            table[3][i] = minus_psi(table[2][i], ring);
        }

        ProjectivePoint result = curve.zero(fp2ring);
        bool started = false;
        for (size_t position = used; position != 0;) {
            --position;
            if (started)
                result = curve.double_point(result, fp2ring);
            for (size_t component = 0; component < 4; ++component) {
                const int digit = digits[component][position];
                if (digit == 0)
                    continue;
                const unsigned magnitude = digit < 0
                    ? static_cast<unsigned>(-digit)
                    : static_cast<unsigned>(digit);
                BLST_ASSERT((magnitude & 1U) != 0);
                BLST_ASSERT((magnitude >> 1) < table_size);
                ProjectivePoint row = table[component][magnitude >> 1];
                if (digit < 0)
                    row.y = fp2ring.standard_negate(row.y);
                if (!started) {
                    result = row;
                    started = true;
                }
                else {
                    result = curve.add_point(result, row, fp2ring);
                }
            }
        }
        blst::vec_zero(digits, sizeof(digits));
        return result;
    }

};

}  // namespace blst::avx

#undef BLST_AVX_SCALAR_MULT_NOINLINE
