#pragma once
#include "ec.hpp"
#include "constexpr_utils.hpp"
#include "../blst/vect.h"
#include "../blst/bytes.h"
#include "../blst/ec_mult.h"
#include "gmp_wrapper.hpp"
#include <array>
#include <cstring>
#include <cstdint>
// Modified from blst

template<class Curve, class Ring, int table_log_size>
class PrecomputedTable {
    public:
    using ProjectivePoint = typename Curve::ProjPoint;

    std::array<ProjectivePoint, (1 << table_log_size)> data;

    PrecomputedTable(const std::array<ProjectivePoint, 1 << table_log_size> &data) : data(data) {}

    PrecomputedTable(const ProjectivePoint &P, const Curve &curve, const Ring &ring) {
        data[0] = P;
        // TODO: more efficient way
        for (int i = 1; i < (1 << table_log_size); i++) {
            if (i % 2 == 1) {
                data[i] = curve.double_point(data[i / 2], ring);
            } else {
                data[i] = curve.add_point(data[i - 1], P, ring);
            }
        }
    }

    std::pair<uint8_t, bool> booth_encode(uint64_t idx) const {
        bool booth_sign = (idx >> table_log_size);
        uint8_t booth_index = booth_sign ? (1 << (table_log_size+1)) - idx : idx;
        return std::make_pair(booth_index, booth_sign);
    }

    ProjectivePoint gather_booth(uint8_t booth_index, bool booth_sign, const Curve &curve, const Ring &ring) const {
        ProjectivePoint result = curve.zero(ring);
        for (uint8_t i = 1; i <= (1 << table_log_size); i++) {
            bool match = (i == booth_index);
            // Replace result with data[i-1] if we match
            result = result.ctime_replace_if(match, data[i - 1]);
        }
        // Replace y with negated y if booth_sign is true (matches Python: n = -val if sign else val)
        result.y = result.y.ctime_replace_if(booth_sign, ring.standard_negate(result.y));
        return result;
    }
};

template<int scalar_bits>
class ScalarForExponent {
    public:
    static constexpr int limbs = ceil_div(scalar_bits, 8);
    std::array<uint8_t, limbs> data;

    ScalarForExponent(const BigInt &scalar) {
        for (int i = 0; i < limbs; i++) {
            data[i] = (scalar >> (i * 8)).to_ulong() & 255;
        }
    }

    ScalarForExponent(const std::array<uint8_t, limbs> &raw_data) : data(raw_data) {}

    /*
    uint32_t rem;
    uint32_t rem_bits;
    uint32_t byte_pos;

    template<int num_bits>
    uint8_t get_next_bits() const {
        if (rem_bits >= num_bits) {
            result = rem & ((1 << num_bits) - 1);
            rem = rem >> num_bits;
            rem_bits -= num_bits;
            return result;
        } else {
            uint8_t next_byte;
            if (byte_pos >= limbs) {
                next_byte = 0;
            } else {
                next_byte = data[byte_pos];
            }
            result = rem | (next_byte << rem_bits);
            rem_bits += 8;
            byte_pos++;
            return get_next_bits<num_bits>();
        }
    }
    */

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
        // Match Python logic: if not carry_possible, assert carry==0 and r=0 (skip carry handling)
        // If carry_possible, handle carry and double
        if constexpr (carry_possible) {
            typename Curve::ProjPoint carry_point = curve.zero(ring).ctime_replace_if(carry, tables[i].data[(1 << wbits) - 1]);
            result = curve.add_point(result, carry_point, ring);
        }
        // When carry_possible is False, carry should be 0 (asserted in Python), so we skip carry handling
    }
    // Double after carry handling only when carry_possible (Python: r = r * 2 only in carry_possible branch)
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
    return result;
}

template<class Curve, class Ring, int scalar_bits, int table_log_size, int num_tables>
typename Curve::ProjPoint table_msm(const std::array<BigInt, num_tables> &scalars, const std::array<PrecomputedTable<Curve, Ring, table_log_size>, num_tables> &tables, const Curve &curve, const Ring &ring) {
    constexpr int wbits = table_log_size;
    constexpr int num_windows = (scalar_bits + wbits) / (wbits + 1);
    constexpr int carry_possible = scalar_bits % (wbits + 1) == 0;
    std::array<std::array<std::pair<uint8_t, bool>, num_windows>, num_tables> booth_indexes;
    typename Curve::ProjPoint result = curve.zero(ring);
    
    for (int i = 0; i < num_tables; i++) {
        bool carry = false;
        for (int j = 0; j < num_windows; j++) {
            uint64_t window_val = (scalars[i] >> (j * (wbits + 1))).to_ulong() & ((1ULL << (wbits + 1)) - 1);
            window_val += carry ? 1 : 0;
            booth_indexes[i][j] = tables[i].booth_encode(window_val);
            carry = booth_indexes[i][j].second;
        }
        // Match Python logic: if not carry_possible, assert carry==0 and r=0 (skip carry handling)
        // If carry_possible, handle carry and double
        if constexpr (carry_possible) {
            typename Curve::ProjPoint carry_point = curve.zero(ring).ctime_replace_if(carry, tables[i].data[(1 << wbits) - 1]);
            result = curve.add_point(result, carry_point, ring);
        }
        // When carry_possible is False, carry should be 0 (asserted in Python), so we skip carry handling
    }
    // Double after carry handling only when carry_possible (Python: r = r * 2 only in carry_possible branch)
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
    return result;
}

template<class Curve, class Ring, int scalar_bits, int table_log_size>
typename Curve::ProjPoint scalar_mult_table(const BigInt &scalar, const typename Curve::ProjPoint &P, const Curve &curve, const Ring &ring) {
    PrecomputedTable<Curve, Ring, table_log_size> table(P, curve, ring);
    return table_msm<Curve, Ring, scalar_bits, table_log_size, 1>(std::array<BigInt, 1>{scalar}, std::array<PrecomputedTable<Curve, Ring, table_log_size>, 1>{table}, curve, ring);
}


// Wrapper functions for blst division operations
// Returns {remainder, quotient} to match BigInt version ordering: {sk_zz_r, sk_zz}
inline std::pair<std::array<uint8_t, 16>, std::array<uint8_t, 16>> div_by_zz(const std::array<uint8_t, 32> &a) {
    union { vec256 l; pow256 s; } val;
    limbs_from_le_bytes(val.l, a.data(), 32);
    blst_div_by_zz(val.l);
    le_bytes_from_limbs(val.s, val.l, 32);
    
    std::array<uint8_t, 16> remainder;
    std::array<uint8_t, 16> quotient;
    
    // blst_div_by_zz stores: val.s is remainder, val.s+16 is quotient
    // But we return {remainder, quotient} to match BigInt version: {sk_zz_r, sk_zz}
    std::memcpy(remainder.data(), val.s, 16);
    std::memcpy(quotient.data(), val.s + 16, 16);
    
    return std::make_pair(remainder, quotient);
}

// Returns {remainder, quotient} to match BigInt version ordering
inline std::pair<std::array<uint8_t, 8>, std::array<uint8_t, 8>> div_by_z(const std::array<uint8_t, 16> &a) {
    // Define types for 128-bit operations (similar to pow256/vec256)
    typedef byte pow128[16];
    typedef limb_t vec128[NLIMBS(128)];
    
    union { vec128 l; pow128 s; } val;
    limbs_from_le_bytes(val.l, a.data(), 16);
    blst_div_by_z(val.l);
    le_bytes_from_limbs(val.s, val.l, 16);
    
    std::array<uint8_t, 8> remainder;
    std::array<uint8_t, 8> quotient;
    
    // blst_div_by_z stores: val.s is remainder, val.s+8 is quotient
    // But we return {remainder, quotient} to match BigInt version
    std::memcpy(remainder.data(), val.s, 8);
    std::memcpy(quotient.data(), val.s + 8, 8);
    
    return std::make_pair(remainder, quotient);
}

template<class Ring>
class GLVMult {
    public:
    using G1Curve = G1<Ring>;
    using ProjectivePoint = typename G1Curve::ProjPoint;
    using AffinePoint = typename G1Curve::AffPoint;

    typename Ring::StandardElement beta;
    typename Ring::StandardElement beta_inv;
    G1Curve curve;

    BigInt zz;

    GLVMult(const G1Curve &curve, const Ring &ring) : curve(curve) {
        beta = ring.from_hex("0x5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a00022e01fffffffefffe");
        zz = BigInt("ac45a4010001a4020000000100000000", 16);
    }

    template<int batch_size>
    std::array<ProjectivePoint, batch_size> batch_endomorphism(const std::array<ProjectivePoint, batch_size> &P, const Ring &ring) const {
        std::array<ProjectivePoint, batch_size> result;
        constexpr int chunk_size = std::min(batch_size, 8);
        static_assert(batch_size % chunk_size == 0, "batch_size must be divisible by chunk_size");
        // TODO: batch reduce
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
        auto [sk_zz_r, sk_zz] = div_by_zz(SK);
        ScalarForExponent<128> scalar_zz_r(sk_zz_r);
        ScalarForExponent<128> scalar_zz(sk_zz);
        std::array<ScalarForExponent<128>, 2> scalars = {scalar_zz_r, scalar_zz};
        PrecomputedTable<G1Curve, Ring, table_log_size> tableP(P, curve, ring);
        // TODO: If z is unchanged by endomorphism, it doesn't really need to be stored twice (helps cache efficiency).  Negation of y could also be more cache space-time trade-off.
        auto table_phi_data = batch_endomorphism<table_size>(tableP.data, ring);
        PrecomputedTable<G1Curve, Ring, table_log_size> tablePhi(table_phi_data);

        std::array<PrecomputedTable<G1Curve, Ring, table_log_size>, 2> tables = {tableP, tablePhi};
        return table_msm<G1Curve, Ring, 128, table_log_size, 2>(scalars, tables, curve, ring);
    }

    template<int table_log_size=4>
    ProjectivePoint multiply(const ProjectivePoint &P, const BigInt &SK, const Ring &ring) const {
        constexpr int table_size = 1 << table_log_size;
        BigInt sk_zz = SK / zz;
        BigInt sk_zz_r = SK % zz;
        std::array<BigInt, 2> scalars = {sk_zz_r, sk_zz};
        PrecomputedTable<G1Curve, Ring, table_log_size> tableP(P, curve, ring);
        // TODO: If z is unchanged by endomorphism, it doesn't really need to be stored twice (helps cache efficiency).  Negation of y could also be more cache space-time trade-off.
        auto table_phi_data = batch_endomorphism<table_size>(tableP.data, ring);
        PrecomputedTable<G1Curve, Ring, table_log_size> tablePhi(table_phi_data);

        std::array<PrecomputedTable<G1Curve, Ring, table_log_size>, 2> tables = {tableP, tablePhi};
        return table_msm<G1Curve, Ring, 128, table_log_size, 2>(scalars, tables, curve, ring);
    }
};

template<class Ring>
class GLSMult {
    public:
    using G2Curve = G2<Ring>;
    using FP2RingType = FP2Ring<Ring>;
    using ProjectivePoint = typename G2Curve::ProjPoint;
    using FP2StandardElement = typename FP2RingType::StandardElement;

    typename Ring::StandardElement psi_x_i;
    typename Ring::StandardElement psi_y_r;
    typename Ring::StandardElement psi_y_i;

    BigInt z;
    BigInt zz;

    GLSMult(const Ring &ring) {
        psi_x_i = ring.from_hex("0x1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaad");
        psi_y_r = ring.from_hex("0x135203e60180a68ee2e9c448d77a2cd91c3dedd930b1cf60ef396489f61eb45e304466cf3e67fa0af1ee7b04121bdea2");
        psi_y_i = ring.from_hex("0x06af0e0437ff400b6831e36d6bd17ffe48395dabc2d3435e77f76e17009241c5ee67992f72ec05f4c81084fbede3cc09");

        z = BigInt("d201000000010000", 16);
        zz = BigInt("ac45a4010001a4020000000100000000", 16);
    }
   
    ProjectivePoint minus_psi(ProjectivePoint P, const FP2RingType &/* fp2ring */, const Ring &ring) const {
        // Unroll; skip multiplication by 0
        // Combine negation step
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
        auto [sk_zz_r, sk_zz] = div_by_zz(SK);
        auto [sk_z_r_r, sk_z_r_q] = div_by_z(sk_zz_r);
        auto [sk_z_q_r, sk_z_q_q] = div_by_z(sk_zz);
        ScalarForExponent<64> scalar_z_r_r(sk_z_r_r);
        ScalarForExponent<64> scalar_z_r_q(sk_z_r_q);
        ScalarForExponent<64> scalar_z_q_r(sk_z_q_r);
        ScalarForExponent<64> scalar_z_q_q(sk_z_q_q);
        std::array<ScalarForExponent<64>, 4> scalars = {scalar_z_r_r, scalar_z_r_q, scalar_z_q_r, scalar_z_q_q};
        PrecomputedTable<G2Curve, FP2RingType, table_log_size> tableP(P, curve, fp2ring);
        std::array<ProjectivePoint, table_size> table_psi_data;
        std::array<ProjectivePoint, table_size> table_psi2_data;
        std::array<ProjectivePoint, table_size> table_psi3_data;
        for (int i = 0; i < table_size; i++) {
            table_psi_data[i] = minus_psi(tableP.data[i], fp2ring, ring);
            table_psi2_data[i] = minus_psi(table_psi_data[i], fp2ring, ring);
            table_psi3_data[i] = minus_psi(table_psi2_data[i], fp2ring, ring);
        }
        PrecomputedTable<G2Curve, FP2RingType, table_log_size> tablePsi(table_psi_data);
        PrecomputedTable<G2Curve, FP2RingType, table_log_size> tablePsi2(table_psi2_data);
        PrecomputedTable<G2Curve, FP2RingType, table_log_size> tablePsi3(table_psi3_data);
        std::array<PrecomputedTable<G2Curve, FP2RingType, table_log_size>, 4> tables = {tableP, tablePsi, tablePsi2, tablePsi3};
        return table_msm<G2Curve, FP2RingType, 64, table_log_size, 4>(scalars, tables, curve, fp2ring);
    }

    template<int table_log_size=4>
    ProjectivePoint multiply(const ProjectivePoint &P, const BigInt &SK, const FP2RingType &fp2ring, const G2Curve &curve, const Ring &ring) const {
        constexpr int table_size = 1 << table_log_size;
        BigInt sk_zz_high = SK / zz;
        BigInt sk_zz_low = SK % zz;
        BigInt sk_z_low_low = sk_zz_low % z;
        BigInt sk_z_low_high = sk_zz_low / z;
        BigInt sk_z_high_low = sk_zz_high % z;
        BigInt sk_z_high_high = sk_zz_high / z;
        std::array<BigInt, 4> scalars = {sk_z_low_low, sk_z_low_high, sk_z_high_low, sk_z_high_high};
        PrecomputedTable<G2Curve, FP2RingType, table_log_size> tableP(P, curve, fp2ring);
        std::array<ProjectivePoint, table_size> table_psi_data;
        std::array<ProjectivePoint, table_size> table_psi2_data;
        std::array<ProjectivePoint, table_size> table_psi3_data;
        for (int i = 0; i < table_size; i++) {
            table_psi_data[i] = minus_psi(tableP.data[i], fp2ring, ring);
            table_psi2_data[i] = minus_psi(table_psi_data[i], fp2ring, ring);
            table_psi3_data[i] = minus_psi(table_psi2_data[i], fp2ring, ring);
        }
        PrecomputedTable<G2Curve, FP2RingType, table_log_size> tablePsi(table_psi_data);
        PrecomputedTable<G2Curve, FP2RingType, table_log_size> tablePsi2(table_psi2_data);
        PrecomputedTable<G2Curve, FP2RingType, table_log_size> tablePsi3(table_psi3_data);
        std::array<PrecomputedTable<G2Curve, FP2RingType, table_log_size>, 4> tables = {tableP, tablePsi, tablePsi2, tablePsi3};
        return table_msm<G2Curve, FP2RingType, 64, table_log_size, 4>(scalars, tables, curve, fp2ring);
    }
};