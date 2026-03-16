#pragma once
#include "vector_impl.hpp"
#include "../precompute/precompute.hpp"
#include <type_traits>
#include <iostream>
typedef unsigned __int128 uint128_t;

template <int in_limbs, int out_limbs, int IN_DIGITS=1, int OUT_DIGITS=1, int INBITS = 52, int OUTBITS = 52>
class RNSMatrix {

    private:
        AVXVector<out_limbs*OUT_DIGITS> rns_mat[in_limbs*IN_DIGITS];
        uint64_t shifted_quotient_estimations[in_limbs*IN_DIGITS];
        AVXVector<out_limbs*OUT_DIGITS> correction;
        AVXVector<out_limbs*OUT_DIGITS> correction_shift;
        BigInt premult;
        BigInt postmult;

        static constexpr int precision = 64;
        static constexpr int mulbits = 52;

    public:

    // Helper to convert uint64_t array to BigInt array
    template<typename T, size_t N>
    static inline std::array<BigInt, N> to_bigint_array(const std::array<T, N> &arr) {
        if constexpr (std::is_same_v<T, BigInt>) {
            return arr; // No conversion needed
        } else {
            std::array<BigInt, N> result;
            for (size_t i = 0; i < N; i++) {
                result[i] = BigInt(arr[i]);
            }
            return result;
        }
    }

    // Template constructor that accepts both BigInt and uint64_t moduli arrays
    template<typename ModuliInT, typename ModuliOutT>
    inline RNSMatrix(
        const std::array<ModuliInT, in_limbs> &moduli_in,
        const std::array<ModuliOutT, out_limbs> &moduli_out,
        const BigInt &target_modulus,
        const BigInt &premult,
        const BigInt &postmult
    ) : premult(premult), postmult(postmult) {
        // Convert to BigInt arrays if needed
        std::array<BigInt, in_limbs> moduli_in_big = to_bigint_array(moduli_in);
        std::array<BigInt, out_limbs> moduli_out_big = to_bigint_array(moduli_out);

        std::array<std::array<uint64_t, out_limbs * OUT_DIGITS>, in_limbs * IN_DIGITS> rns_values_tf{};
        std::array<uint64_t, out_limbs * OUT_DIGITS> correction_vec{};
        std::array<uint64_t, out_limbs * OUT_DIGITS> correction_shift_vec{};
        std::array<BigInt, in_limbs * IN_DIGITS> sqe_tf{};

        gen_precompute_wide<in_limbs, out_limbs, IN_DIGITS, INBITS, OUT_DIGITS, OUTBITS, precision>(
            moduli_in_big, moduli_out_big, target_modulus, premult, postmult,
            rns_values_tf, correction_vec, sqe_tf, false  // Debug output disabled
        );

        // Load into members
        for (int i = 0; i < in_limbs * IN_DIGITS; i++) {
            rns_mat[i].load((uint64_t*)rns_values_tf[i].data());
            shifted_quotient_estimations[i] = static_cast<uint64_t>(sqe_tf[i].to_ulong());
        }
        correction.load((uint64_t*)correction_vec.data());
       
        // Cannot use shifted version if binary output.
        if (moduli_out_big[0] < BigInt(1) << 64) {
            for (int i = 0; i < out_limbs * OUT_DIGITS; i++) {
                correction_shift_vec[i] = ((uint128_t)(correction_vec[i]) * (1ull << mulbits)) % moduli_out_big[i].to_ulong();
            }
        }
        correction_shift.load((uint64_t*)correction_shift_vec.data());
    }

    inline void accumulate_loop(AVXVector<out_limbs * OUT_DIGITS> &out_hi, AVXVector<out_limbs * OUT_DIGITS> &out_lo, const AVXVector<in_limbs * IN_DIGITS> &residues, uint128_t &k_raw) const {
        // Need to allocate space with padding, although we should only iterate over the actual limbs.
        std::array<uint64_t, AVXVector<in_limbs * IN_DIGITS>::VEC_LIMBS * AVXVector<in_limbs * IN_DIGITS>::LIMBS_PER_VEC> scalars;
        residues.store(scalars.data());
        for (int i = 0; i < in_limbs * IN_DIGITS; i++) {
            uint64_t s = scalars[i];
            AVXScalar sv = Scalar(s);
            out_hi = out_hi.mulhi_scalar(rns_mat[i], sv);
            out_lo = out_lo.mullo_scalar(rns_mat[i], sv);
            k_raw += (uint128_t)(s) * (uint128_t)(shifted_quotient_estimations[i]);
        }
    }

    inline uint64_t add_correction(AVXVector<out_limbs * OUT_DIGITS> &out_hi, AVXVector<out_limbs * OUT_DIGITS> &out_lo, uint128_t raw, const AVXVector<out_limbs * OUT_DIGITS> &c_value) const {
        uint64_t k = raw >> precision;
        /*std::cout << "  [DEBUG C++ add_correction] k_raw=" << (uint64_t)(raw >> 64) << ":" << (uint64_t)raw << ", k=" << k << std::endl;*/
        AVXScalar sk = Scalar(k);
        out_hi = out_hi.mulhi_scalar(c_value, sk);
        out_lo = out_lo.mullo_scalar(c_value, sk);
        sk = scalar_shift(sk, mulbits);
        out_hi = out_hi.mullo_scalar(c_value, sk);
        AVXVector<out_limbs * OUT_DIGITS> lp = AVXVector<out_limbs * OUT_DIGITS>(0).mulhi_scalar(c_value, sk);
        out_hi = out_hi.add(lp.slli(mulbits));
        return k;
    }

    template<int batch>
    inline void add_correction_batch(std::array<AVXVector<out_limbs * OUT_DIGITS>, batch> &out_hi, std::array<AVXVector<out_limbs * OUT_DIGITS>, batch> &out_lo, std::array<uint128_t, batch> &k_raw, const AVXVector<out_limbs * OUT_DIGITS> &c_value) const {
        std::array<AVXScalar, batch> sk;
        for (int i = 0; i < batch; i++) {
            sk[i] = Scalar(k_raw[i] >> precision);
        }
        for (int i = 0; i < batch; i++) {
            out_hi[i] = out_hi[i].mulhi_scalar(c_value, sk[i]);
        }
        for (int i = 0; i < batch; i++) {
            out_lo[i] = out_lo[i].mullo_scalar(c_value, sk[i]);
        }
        for (int i = 0; i < batch; i++) {
            sk[i] = scalar_shift(sk[i], mulbits);
        }
        for (int i = 0; i < batch; i++) {
            out_hi[i] = out_hi[i].mullo_scalar(c_value, sk[i]);
        }
        std::array<AVXVector<out_limbs * OUT_DIGITS>, batch> lp;
        for (int i = 0; i < batch; i++) {
            lp[i] = AVXVector<out_limbs * OUT_DIGITS>(0).mulhi_scalar(c_value, sk[i]);
        }
        for (int i = 0; i < batch; i++) {
            lp[i] = lp[i].slli(mulbits);
        }
        for (int i = 0; i < batch; i++) {
            out_hi[i] = out_hi[i].add(lp[i]);
        }
    }

    
    template<int batch>
    inline void add_correction_batch2(std::array<AVXVector<out_limbs * OUT_DIGITS>, batch> &out_hi, std::array<AVXVector<out_limbs * OUT_DIGITS>, batch> &out_lo, std::array<uint128_t, batch> &k_raw, const AVXVector<out_limbs * OUT_DIGITS> &c_value, const AVXVector<out_limbs * OUT_DIGITS> &c_value_shift) const {
        std::array<AVXScalar, batch> sk;
        for (int i = 0; i < batch; i++) {
            sk[i] = Scalar(k_raw[i] >> precision);
        }
        for (int i = 0; i < batch; i++) {
            out_hi[i] = out_hi[i].mulhi_scalar(c_value, sk[i]);
        }
        for (int i = 0; i < batch; i++) {
            out_lo[i] = out_lo[i].mullo_scalar(c_value, sk[i]);
        }
        for (int i = 0; i < batch; i++) {
            sk[i] = scalar_shift(sk[i], mulbits);
        }
        // k = sk_shifted[i] * 2^52 + sk[i]
        // c_value_shift = (c_value * 2^52) % modulus
        for (int i = 0; i < batch; i++) {
            out_lo[i] = out_lo[i].mullo_scalar(c_value_shift, sk[i]);
        }
        for (int i = 0; i < batch; i++) {
            out_hi[i] = out_hi[i].mulhi_scalar(c_value_shift, sk[i]);
        }
    }

    inline void rns_reduce_raw(const AVXVector<in_limbs * IN_DIGITS> &residues1, AVXVector<out_limbs * OUT_DIGITS> &out_hi, AVXVector<out_limbs * OUT_DIGITS> &out_lo) const {
        uint128_t k_raw = 0;
        accumulate_loop(out_hi, out_lo, residues1, k_raw);
        (void)add_correction(out_hi, out_lo, k_raw, correction);
    }

    template<int batch>
    inline void accumulate_loop_batch_tree_1_level(
        const std::array<std::array<uint64_t, AVXVector<in_limbs * IN_DIGITS>::VEC_LIMBS * AVXVector<in_limbs * IN_DIGITS>::LIMBS_PER_VEC>, batch> scalars,
        std::array<uint128_t, batch> &k_raw,
        std::array<AVXVector<out_limbs * OUT_DIGITS>, batch> &out_hi,
        std::array<AVXVector<out_limbs * OUT_DIGITS>, batch> &out_lo
    ) const {
        std::array<AVXVector<out_limbs * OUT_DIGITS>, batch> out_hi_tmp{};
        std::array<AVXVector<out_limbs * OUT_DIGITS>, batch> out_lo_tmp{};
        static_assert(in_limbs * IN_DIGITS % 2 == 0, "in_limbs * IN_DIGITS must be even");
        int half_limbs = in_limbs * IN_DIGITS / 2;
        for (int i = 0; i < half_limbs; i++) {
            for (int j = 0; j < batch; j++) {
                uint64_t s = scalars[j][i];
                AVXScalar sv = Scalar(s);
                out_hi[j] = out_hi[j].mulhi_scalar(rns_mat[i], sv);
                out_lo[j] = out_lo[j].mullo_scalar(rns_mat[i], sv);
                k_raw[j] += (uint128_t)(s) * (uint128_t)(shifted_quotient_estimations[i]);

                s = scalars[j][i + half_limbs];
                sv = Scalar(s);
                out_hi_tmp[j] = out_hi_tmp[j].mulhi_scalar(rns_mat[i + half_limbs], sv);
                out_lo_tmp[j] = out_lo_tmp[j].mullo_scalar(rns_mat[i + half_limbs], sv);
                k_raw[j] += (uint128_t)(s) * (uint128_t)(shifted_quotient_estimations[i + half_limbs]);
            }
        }
        for (int j = 0; j < batch; j++) {
            out_hi[j] = out_hi[j].add(out_hi_tmp[j]);
            out_lo[j] = out_lo[j].add(out_lo_tmp[j]);
        }
    }

    template<int batch, int version=1>
    inline __attribute__((always_inline)) void rns_reduce_raw_batch(const std::array<AVXVector<in_limbs * IN_DIGITS>, batch> *residues1, std::array<AVXVector<out_limbs * OUT_DIGITS>, batch> &out_hi, std::array<AVXVector<out_limbs * OUT_DIGITS>, batch> &out_lo) const {
        using U64type = std::array<std::array<uint64_t, AVXVector<in_limbs * IN_DIGITS>::VEC_LIMBS * AVXVector<in_limbs * IN_DIGITS>::LIMBS_PER_VEC>, batch>;
        const U64type* scalars = reinterpret_cast<const U64type*>(residues1);
        std::array<uint128_t, batch> k_raw = {};
        if constexpr (version == 3) {
            accumulate_loop_batch_tree_1_level<batch>(scalars, k_raw, out_hi, out_lo);
        } else {
            for (int i = 0; i < in_limbs * IN_DIGITS; i++) {
                for (int j = 0; j < batch; j++) {
                    uint64_t s = (*scalars)[j][i];
                    AVXScalar sv = Scalar(s);
                    k_raw[j] += (uint128_t)(s) * (uint128_t)(shifted_quotient_estimations[i]);
                    out_hi[j] = out_hi[j].mulhi_scalar(rns_mat[i], sv);
                    out_lo[j] = out_lo[j].mullo_scalar(rns_mat[i], sv);
                }
            }
        }
        if constexpr (version == 2) {
            add_correction_batch2<batch>(out_hi, out_lo, k_raw, correction, correction_shift);
        } else {
            add_correction_batch<batch>(out_hi, out_lo, k_raw, correction);
        }
        // for (int i = 0; i < batch; i++) {
        //     (void)add_correction(out_hi[i], out_lo[i], k_raw[i], correction);
        // }
    }
    
    /*
    template<int batch, int version=1>
    inline __attribute__((always_inline)) std::pair<std::array<AVXVector<out_limbs * OUT_DIGITS>, batch>, std::array<AVXVector<out_limbs * OUT_DIGITS>, batch>> rns_reduce_raw_batch(
        const std::array<AVXVector<in_limbs * IN_DIGITS>, batch> &residues1,
        std::array<AVXVector<out_limbs * OUT_DIGITS>, batch> &out_hi_in,
        std::array<AVXVector<out_limbs * OUT_DIGITS>, batch> &out_lo_in) const {
        std::array<std::array<uint64_t, AVXVector<in_limbs * IN_DIGITS>::VEC_LIMBS * AVXVector<in_limbs * IN_DIGITS>::LIMBS_PER_VEC>, batch> scalars;
        // AVXBatch<in_limbs * IN_DIGITS, batch> scalars;
        std::array<AVXVector<out_limbs * OUT_DIGITS>, batch> out_hi;
        std::array<AVXVector<out_limbs * OUT_DIGITS>, batch> out_lo;
        std::array<uint64_t, batch> k;
        for (int i = 0; i < batch; i++) {
            residues1[i].store(scalars[i].data());
            out_hi[i] = out_hi_in[i];
            out_lo[i] = out_lo_in[i];
        }
        for (int j = 0; j < batch; j++) {
            uint128_t k_raw = 0;
            for (int i = 0; i < in_limbs * IN_DIGITS; i++) {
                uint64_t s = scalars[j][i];
                k_raw += (uint128_t)(s) * (uint128_t)(shifted_quotient_estimations[i]);
            }
            k[j] = k_raw >> precision;
        }
        #pragma clang loop unroll (enable)
        for (int j = 0; j < 12; j++) {
            AVXScalar sk = Scalar(k[j]);
            out_hi[j] = out_hi[j].mulhi_scalar(correction, sk);
            out_lo[j] = out_lo[j].mullo_scalar(correction, sk);
            sk = scalar_shift(sk, mulbits);
            out_lo[j] = out_lo[j].mullo_scalar(correction_shift, sk);
            out_hi[j] = out_hi[j].mulhi_scalar(correction_shift, sk);
        }
        for (int i = 0; i < in_limbs * IN_DIGITS; i++) {
            #pragma clang loop unroll (enable)
            for (int j = 0; j < batch; j++) {
                uint64_t s = scalars[j][i];
                AVXScalar sv = Scalar(s);
                out_hi[j] = out_hi[j].mulhi_scalar(rns_mat[i], sv);
                out_lo[j] = out_lo[j].mullo_scalar(rns_mat[i], sv);
            }
        }
        for (int i = 0; i < batch; i++) {
            out_hi_in[i] = out_hi[i];
            out_lo_in[i] = out_lo[i];
        }
    }
    */

    template<class Reducer>
    inline AVXVector<out_limbs * OUT_DIGITS> rns_reduce_with_acc(const AVXVector<in_limbs * IN_DIGITS> &residues1, const AVXVector<out_limbs * OUT_DIGITS> &acc1, const AVXVector<out_limbs * OUT_DIGITS> &acc2, const Reducer& reducer) const {
        AVXVector<out_limbs * OUT_DIGITS> out_hi = AVXVector<out_limbs * OUT_DIGITS>(0);
        AVXVector<out_limbs * OUT_DIGITS> out_lo = AVXVector<out_limbs * OUT_DIGITS>(0);
        out_hi = out_hi.mulhi(acc1, acc2);
        out_lo = out_lo.mullo(acc1, acc2);
        rns_reduce_raw(residues1, out_hi, out_lo);
        return reducer.reduce_full(out_hi, out_lo);
    }

    template<class Reducer>
    inline AVXVector<out_limbs * OUT_DIGITS> rns_reduce(const AVXVector<in_limbs * IN_DIGITS> &residues1, const Reducer& reducer) const {
        /*// Debug: print scalars after word_reinterpret
        std::array<uint64_t, in_limbs * IN_DIGITS> scalars;
        residues1.store((uint64_t*)scalars.data());
        std::cout << "  [DEBUG C++ rns_reduce] scalars=[";
        for (int i = 0; i < in_limbs * IN_DIGITS; i++) {
            std::cout << scalars[i];
            if (i < in_limbs * IN_DIGITS - 1) std::cout << ", ";
        }
        std::cout << "]" << std::endl;*/
        
        AVXVector<out_limbs * OUT_DIGITS> out_hi = AVXVector<out_limbs * OUT_DIGITS>(0);
        AVXVector<out_limbs * OUT_DIGITS> out_lo = AVXVector<out_limbs * OUT_DIGITS>(0);
        rns_reduce_raw(residues1, out_hi, out_lo);
        /*std::array<uint64_t, out_limbs * OUT_DIGITS> final_hi, final_lo;
        out_hi.store((uint64_t*)final_hi.data());
        out_lo.store((uint64_t*)final_lo.data());
        std::cout << "  [DEBUG C++ rns_reduce] before reducer.reduce_full: out_hi=[";
        for (int i = 0; i < out_limbs * OUT_DIGITS; i++) {
            std::cout << final_hi[i];
            if (i < out_limbs * OUT_DIGITS - 1) std::cout << ", ";
        }
        std::cout << "], out_lo=[";
        for (int i = 0; i < out_limbs * OUT_DIGITS; i++) {
            std::cout << final_lo[i];
            if (i < out_limbs * OUT_DIGITS - 1) std::cout << ", ";
        }
        std::cout << "]" << std::endl;*/
        AVXVector<out_limbs * OUT_DIGITS> result = reducer.reduce_full(out_hi, out_lo);
        /*std::array<uint64_t, out_limbs * OUT_DIGITS> result_arr;
        result.store((uint64_t*)result_arr.data());
        std::cout << "  [DEBUG C++ rns_reduce] final result=[";
        for (int i = 0; i < out_limbs * OUT_DIGITS; i++) {
            std::cout << result_arr[i];
            if (i < out_limbs * OUT_DIGITS - 1) std::cout << ", ";
        }
        std::cout << "]" << std::endl;*/
        return result;
    }

    // Test accessors
    inline void get_rns_mat(std::array<std::array<uint64_t, out_limbs * OUT_DIGITS>, in_limbs * IN_DIGITS> &out) const {
        for (int i = 0; i < in_limbs * IN_DIGITS; i++) {
            rns_mat[i].store((uint64_t*)out[i].data());
        }
    }
    
    inline void get_correction(std::array<uint64_t, out_limbs * OUT_DIGITS> &out) const {
        correction.store((uint64_t*)out.data());
    }
    
    inline void get_shifted_quotient_estimations(std::array<uint64_t, in_limbs * IN_DIGITS> &out) const {
        for (int i = 0; i < in_limbs * IN_DIGITS; i++) {
            out[i] = shifted_quotient_estimations[i];
        }
    }
    
    inline const BigInt& get_premult() const {
        return premult;
    }
    
    inline const BigInt& get_postmult() const {
        return postmult;
    }
};