#pragma once
#include "changebase.hpp"
#include "vector_impl.hpp"
#include "../reduction/montgomery.hpp"
#include "conversion.hpp"
#include <tuple>


namespace blst::avx {

template <template <int> class Reduction, int limbs>
inline AVXVector<limbs> ele_mult(const AVXVector<limbs> &a, const AVXVector<limbs> & b, const Reduction<limbs> &r) {
    auto ab_lo = AVXVector<limbs>(0).mullo(a, b);
    auto ab_hi = AVXVector<limbs>(0).mulhi(a, b);
    return r.reduce_small(ab_hi, ab_lo);
}

template<template <int> class Reduction, int limbs1, int limbs2>
class IntRNS2 {
    
    public:

    const RNSMatrix<limbs1, limbs2> r1;
    const RNSMatrix<limbs2, limbs1> r2;
    const Reduction<limbs1> m1;
    const Reduction<limbs2> m2;


    inline IntRNS2(
        const std::array<uint64_t, limbs1>& moduli_in_m1,
        const std::array<uint64_t, limbs2>& moduli_in_m2,
        const std::array<std::array<uint64_t, limbs2>, limbs1>& r1_rns_mat,
        const std::array<uint64_t, limbs2>& r1_correction,
        const std::array<uint64_t, limbs1>& r1_sqe,
        const std::array<std::array<uint64_t, limbs1>, limbs2>& r2_rns_mat,
        const std::array<uint64_t, limbs1>& r2_correction,
        const std::array<uint64_t, limbs2>& r2_sqe
    )
    : r1(r1_rns_mat, r1_correction, r1_sqe, moduli_in_m2)
    , r2(r2_rns_mat, r2_correction, r2_sqe, moduli_in_m1)
    , m1(moduli_in_m1)
    , m2(moduli_in_m2)
    {}

    inline void mul_reduce(const AVXVector<limbs1> &a1, const AVXVector<limbs2> &a2, const AVXVector<limbs1> &b1, const AVXVector<limbs2> &b2, AVXVector<limbs1> &out1, AVXVector<limbs2> &out2) const {
        auto ab_m1 = ele_mult(a1, b1, m1);
        out2 = r1.rns_reduce_with_acc(ab_m1, a2, b2, m2);
        out1 = r2.rns_reduce(out2, m1);
    }

    inline void mul_wide(const AVXVector<limbs1> &a1, const AVXVector<limbs2> &a2, const AVXVector<limbs1> &b1, const AVXVector<limbs2> &b2, AVXVector<limbs1> &out1_lo, AVXVector<limbs1> &out1_hi, AVXVector<limbs2> &out2_lo, AVXVector<limbs2> &out2_hi) const {
        out1_lo = AVXVector<limbs1>(0).mullo(a1, b1);
        out1_hi = AVXVector<limbs1>(0).mulhi(a1, b1);
        out2_lo = AVXVector<limbs2>(0).mullo(a2, b2);
        out2_hi = AVXVector<limbs2>(0).mulhi(a2, b2);
    }

    inline void mul_acc_wide(const AVXVector<limbs1> &a1, const AVXVector<limbs2> &a2, const AVXVector<limbs1> &b1, const AVXVector<limbs2> &b2, AVXVector<limbs1> &out1_lo, AVXVector<limbs1> &out1_hi, AVXVector<limbs2> &out2_lo, AVXVector<limbs2> &out2_hi) const {
        out1_lo = out1_lo.mullo(a1, b1);
        out1_hi = out1_hi.mulhi(a1, b1);
        out2_lo = out2_lo.mullo(a2, b2);
        out2_hi = out2_hi.mulhi(a2, b2);
    }

    inline void reduce_wide(const AVXVector<limbs1> &a1_lo, const AVXVector<limbs1> &a1_hi, const AVXVector<limbs2> &a2_lo, const AVXVector<limbs2> &a2_hi, AVXVector<limbs1> &out1, AVXVector<limbs2> &out2) const {
        auto a1 = m1.reduce_small(a1_hi, a1_lo);
        // need to make copy of a2_hi and a2_lo to avoid modifying the input
        auto a2_hi_copy = a2_hi;
        auto a2_lo_copy = a2_lo;
        r1.rns_reduce_raw(a1, a2_hi_copy, a2_lo_copy);
        out2 = m2.reduce_full(a2_hi_copy, a2_lo_copy);
        out1 = r2.rns_reduce(out2, m1);
    }

    template<int batch>
    inline void reduce_wide_batch(const std::array<AVXVector<limbs1>, batch> &a1_lo, const std::array<AVXVector<limbs1>, batch> &a1_hi, const std::array<AVXVector<limbs2>, batch> &a2_lo, const std::array<AVXVector<limbs2>, batch> &a2_hi, std::array<AVXVector<limbs2>, batch> &out2) const {
        std::array<AVXVector<limbs1>, batch> a1;
        std::array<AVXVector<limbs2>, batch> a2_hi_copy;
        std::array<AVXVector<limbs2>, batch> a2_lo_copy;
        a1 = m1.template reduce_small_batch<batch>(a1_hi, a1_lo);
        for (int i = 0; i < batch; i++) {
            a2_hi_copy[i] = a2_hi[i];
            a2_lo_copy[i] = a2_lo[i];
        }
        r1.template rns_reduce_raw_batch<batch>(&a1, a2_hi_copy, a2_lo_copy);
        out2 = m2.template reduce_full_batch<batch>(a2_hi_copy, a2_lo_copy);
    }

    template<int batch>
    inline void batch_expand_wide(const std::array<AVXVector<limbs2>, batch> &out2, std::array<AVXVector<limbs1>, batch> &out1) const {
        std::array<AVXVector<limbs1>, batch> z_hi;
        std::array<AVXVector<limbs1>, batch> z_lo;
        for (int i = 0; i < batch; i++) {
            z_hi[i] = AVXVector<limbs1>(0);
            z_lo[i] = AVXVector<limbs1>(0);
        }
        r2.template rns_reduce_raw_batch<batch>(&out2, z_hi, z_lo);
        // Do I also need reduce_full_batch to amortize loading of Montgomery factors?
        out1 = m1.template reduce_full_batch<batch>(z_hi, z_lo);
    }
};

template<class T, std::size_t l1, std::size_t l2>
std::array<T, l1 + l2> concat(const std::array<T, l1> &a, const std::array<T, l2> &b) {
    std::array<T, l1 + l2> result;
    for (std::size_t i = 0; i < l1; i++) {
        result[i] = a[i];
    }
    for (std::size_t i = 0; i < l2; i++) {
        result[l1 + i] = b[i];
    }
    return result;
}

}  // namespace blst::avx
