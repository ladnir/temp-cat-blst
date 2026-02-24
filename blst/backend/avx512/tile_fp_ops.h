//#ifndef BLST_BACKEND_AVX512_TILE_FP_OPS_H
//#define BLST_BACKEND_AVX512_TILE_FP_OPS_H
//
//#include "blst/backend/avx512/tile_types.h"
//#include "blst/backend/avx512/fp_ifma_ops.h"
//#include "blst/fields/Fp.h"
//
//namespace blst::tile
//{
//#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)
//    constexpr u64 kAvxWidth = kTileWidth;
//    constexpr limb_t kAvxLimbMask = (((limb_t)1) << kTileLimbBits) - 1;
//
//    BLST_HD inline void conv_48to64_fp_8x1w(cc512i* r, const cc512i* a)
//    {
//        const cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
//        const cc512i a4 = a[4], a5 = a[5], a6 = a[6], a7 = a[7];
//        r[0] = _cc512_or_si512(_cc512_slli_epi64(a1, 48), a0);
//        r[1] = _cc512_or_si512(_cc512_slli_epi64(a2, 32), _cc512_srli_epi64(a1, 16));
//        r[2] = _cc512_or_si512(_cc512_slli_epi64(a3, 16), _cc512_srli_epi64(a2, 32));
//        r[3] = _cc512_or_si512(_cc512_slli_epi64(a5, 48), a4);
//        r[4] = _cc512_or_si512(_cc512_slli_epi64(a6, 32), _cc512_srli_epi64(a5, 16));
//        r[5] = _cc512_or_si512(_cc512_slli_epi64(a7, 16), _cc512_srli_epi64(a6, 32));
//    }
//
//    BLST_HD inline void conv_64to48_fp_8x1w(cc512i* r, const cc512i* a)
//    {
//        const cc512i a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
//        const cc512i a4 = a[4], a5 = a[5];
//        const cc512i bmask = _cc512_set1_epi64((long long)kAvxLimbMask);
//        r[0] = _cc512_and_si512(a0, bmask);
//        r[1] = _cc512_and_si512(_cc512_or_si512(_cc512_srli_epi64(a0, 48), _cc512_slli_epi64(a1, 16)), bmask);
//        r[2] = _cc512_and_si512(_cc512_or_si512(_cc512_srli_epi64(a1, 32), _cc512_slli_epi64(a2, 32)), bmask);
//        r[3] = _cc512_and_si512(_cc512_srli_epi64(a2, 16), bmask);
//        r[4] = _cc512_and_si512(a3, bmask);
//        r[5] = _cc512_and_si512(_cc512_or_si512(_cc512_srli_epi64(a3, 48), _cc512_slli_epi64(a4, 16)), bmask);
//        r[6] = _cc512_and_si512(_cc512_or_si512(_cc512_srli_epi64(a4, 32), _cc512_slli_epi64(a5, 32)), bmask);
//        r[7] = _cc512_and_si512(_cc512_srli_epi64(a5, 16), bmask);
//    }
//
//    BLST_HD inline void fp_tile_unpack_all(Fp out[kAvxWidth], const FpTileRegAvx512& in)
//    {
//        cc512i w64[NLIMBS(384)];
//        alignas(64) u64 lanes[NLIMBS(384)][kAvxWidth];
//        conv_48to64_fp_8x1w(w64, in.data);
//        for (u64 limb = 0; limb < NLIMBS(384); ++limb) {
//            _cc512_store_si512((void*)lanes[limb], w64[limb]);
//        }
//        for (u64 lane = 0; lane < kAvxWidth; ++lane) {
//            for (u64 limb = 0; limb < NLIMBS(384); ++limb) {
//                out[lane][limb] = lanes[limb][lane];
//            }
//        }
//    }
//
//    BLST_HD inline void fp_tile_pack_all(FpTileRegAvx512& out, const Fp in[kAvxWidth])
//    {
//        cc512i w64[NLIMBS(384)];
//        alignas(64) u64 lanes[NLIMBS(384)][kAvxWidth];
//        for (u64 lane = 0; lane < kAvxWidth; ++lane) {
//            for (u64 limb = 0; limb < NLIMBS(384); ++limb) {
//                lanes[limb][lane] = in[lane][limb];
//            }
//        }
//        for (u64 limb = 0; limb < NLIMBS(384); ++limb) {
//            w64[limb] = _cc512_load_si512((const void*)lanes[limb]);
//        }
//        conv_64to48_fp_8x1w(out.data, w64);
//    }
//
//    BLST_HD inline void fp_tile_pack_lanes(FpTileRegAvx512& out, const Fp* lanes_in, u64 lane_count)
//    {
//        cc512i w64[NLIMBS(384)];
//        alignas(64) u64 lanes[NLIMBS(384)][kAvxWidth] = {};
//        for (u64 lane = 0; lane < lane_count; ++lane) {
//            for (u64 limb = 0; limb < NLIMBS(384); ++limb) {
//                lanes[limb][lane] = lanes_in[lane][limb];
//            }
//        }
//        for (u64 limb = 0; limb < NLIMBS(384); ++limb) {
//            w64[limb] = _cc512_load_si512((const void*)lanes[limb]);
//        }
//        conv_64to48_fp_8x1w(out.data, w64);
//    }
//
//    BLST_HD inline void fp_tile_pack_splat(FpTileRegAvx512& out, const Fp& in)
//    {
//        cc512i w64[NLIMBS(384)];
//        for (u64 limb = 0; limb < NLIMBS(384); ++limb) {
//            w64[limb] = _cc512_set1_epi64((long long)in[limb]);
//        }
//        conv_64to48_fp_8x1w(out.data, w64);
//    }
//
//    BLST_HD inline void fp_tile_mem_zero(FpTileMemViewAvx512& mem)
//    {
//        for (u64 limb = 0; limb < FpTileMemViewAvx512::kLimbs; ++limb) {
//            limb_t* row = mem.row(limb);
//            for (u64 lane = 0; lane < kAvxWidth; ++lane) {
//                row[lane] = 0;
//            }
//        }
//    }
//
//    BLST_HD inline void fp_tile_reg_zero(FpTileRegAvx512& reg)
//    {
//        const cc512i z = _cc512_setzero_si512();
//        for (u64 limb = 0; limb < FpTileRegAvx512::kLimbs; ++limb) {
//            reg[limb] = z;
//        }
//    }
//
//    BLST_HD inline void fp_tile_load(FpTileRegAvx512& reg, const FpTileMemViewAvx512& mem)
//    {
//        const cc512i mask = _cc512_set1_epi64((long long)kAvxLimbMask);
//        for (u64 limb = 0; limb < FpTileMemViewAvx512::kLimbs; ++limb) {
//            reg[limb] = _cc512_and_si512(_cc512_loadu_si512((const void*)mem.row(limb)), mask);
//        }
//    }
//
//    BLST_HD inline void fp_tile_store(FpTileMemViewAvx512& mem, const FpTileRegAvx512& reg)
//    {
//        const cc512i mask = _cc512_set1_epi64((long long)kAvxLimbMask);
//        for (u64 limb = 0; limb < FpTileMemViewAvx512::kLimbs; ++limb) {
//            _cc512_storeu_si512((void*)mem.row(limb), _cc512_and_si512(reg[limb], mask));
//        }
//    }
//
//    BLST_HD inline void fp_lane_unpack(Fp& out, const FpTileRegAvx512& in, u64 lane)
//    {
//        Fp lanes[kAvxWidth];
//        fp_tile_unpack_all(lanes, in);
//        for (u64 i = 0; i < NLIMBS(384); ++i) {
//            out[i] = lanes[lane][i];
//        }
//    }
//
//    BLST_HD inline void fp_lane_pack(FpTileRegAvx512& out, u64 lane, const Fp& in)
//    {
//        Fp lanes[kAvxWidth];
//        fp_tile_unpack_all(lanes, out);
//        for (u64 i = 0; i < NLIMBS(384); ++i) {
//            lanes[lane][i] = in[i];
//        }
//        fp_tile_pack_all(out, lanes);
//    }
//
//    // Keep lane selection in-vector to avoid scalar unpack/repack in hot loops.
//    BLST_HD inline void fp_tile_broadcast_lane(FpTileRegAvx512& out, const FpTileRegAvx512& in, u64 lane)
//    {
//        const cc512i idx = _cc512_set1_epi64((long long)lane);
//        for (u64 limb = 0; limb < FpTileRegAvx512::kLimbs; ++limb) {
//            out[limb] = _cc512_permutexvar_epi64(idx, in[limb]);
//        }
//    }
//
//    BLST_HD inline void fp2_tile_broadcast_lane(Fp2TileRegAvx512& out, const Fp2TileRegAvx512& in, u64 lane)
//    {
//        fp_tile_broadcast_lane(out[0], in[0], lane);
//        fp_tile_broadcast_lane(out[1], in[1], lane);
//    }
//
//    BLST_HD inline void fp6_tile_broadcast_lane(Fp6TileRegAvx512& out, const Fp6TileRegAvx512& in, u64 lane)
//    {
//        fp2_tile_broadcast_lane(out[0], in[0], lane);
//        fp2_tile_broadcast_lane(out[1], in[1], lane);
//        fp2_tile_broadcast_lane(out[2], in[2], lane);
//    }
//
//    BLST_HD inline void fp_tile_add(FpTileRegAvx512& out, const FpTileRegAvx512& a, const FpTileRegAvx512& b)
//    {
//        add_fp_8x1w(out.data, a.data, b.data);
//    }
//
//    BLST_HD inline void fp_tile_sub(FpTileRegAvx512& out, const FpTileRegAvx512& a, const FpTileRegAvx512& b)
//    {
//        sub_fp_8x1w(out.data, a.data, b.data);
//    }
//
//    BLST_HD inline void fp_tile_add_mod(FpTileRegAvx512& out, const FpTileRegAvx512& a, const FpTileRegAvx512& b, const limb_t*)
//    {
//        // IFMA add kernel already reduces modulo p.
//        fp_tile_add(out, a, b);
//    }
//
//    BLST_HD inline void fp_tile_sub_mod(FpTileRegAvx512& out, const FpTileRegAvx512& a, const FpTileRegAvx512& b, const limb_t*)
//    {
//        // IFMA sub kernel already reduces modulo p.
//        fp_tile_sub(out, a, b);
//    }
//
//    BLST_HD inline void fp_tile_neg_mod(FpTileRegAvx512& out, const FpTileRegAvx512& a, const limb_t*)
//    {
//        FpTileRegAvx512 zero = {};
//        fp_tile_sub(out, zero, a);
//    }
//
//    BLST_HD inline void fp_tile_mul_mont(FpTileRegAvx512& out, const FpTileRegAvx512& a, const FpTileRegAvx512& b, const limb_t*, limb_t)
//    {
//        mul_redc_fp_8x1w(out.data, a.data, b.data);
//    }
//
//    BLST_HD inline void fp_tile_sqr_mont(FpTileRegAvx512& out, const FpTileRegAvx512& a, const limb_t*, limb_t)
//    {
//        mul_redc_fp_8x1w(out.data, a.data, a.data);
//    }
//
//    BLST_HD inline void fp2_tile_add(Fp2TileRegAvx512& out, const Fp2TileRegAvx512& a, const Fp2TileRegAvx512& b) { fp_tile_add(out[0], a[0], b[0]); fp_tile_add(out[1], a[1], b[1]); }
//    BLST_HD inline void fp2_tile_sub(Fp2TileRegAvx512& out, const Fp2TileRegAvx512& a, const Fp2TileRegAvx512& b) { fp_tile_sub(out[0], a[0], b[0]); fp_tile_sub(out[1], a[1], b[1]); }
//    BLST_HD inline void fp2_tile_add_mod(Fp2TileRegAvx512& out, const Fp2TileRegAvx512& a, const Fp2TileRegAvx512& b, const limb_t* m) { fp_tile_add_mod(out[0], a[0], b[0], m); fp_tile_add_mod(out[1], a[1], b[1], m); }
//    BLST_HD inline void fp2_tile_sub_mod(Fp2TileRegAvx512& out, const Fp2TileRegAvx512& a, const Fp2TileRegAvx512& b, const limb_t* m) { fp_tile_sub_mod(out[0], a[0], b[0], m); fp_tile_sub_mod(out[1], a[1], b[1], m); }
//    BLST_HD inline void fp2_tile_neg(Fp2TileRegAvx512& out, const Fp2TileRegAvx512& a, const limb_t* m) { fp_tile_neg_mod(out[0], a[0], m); fp_tile_neg_mod(out[1], a[1], m); }
//    BLST_HD inline void fp2_tile_conjugate(Fp2TileRegAvx512& out, const Fp2TileRegAvx512& a, const limb_t* m) { out[0] = a[0]; fp_tile_neg_mod(out[1], a[1], m); }
//    BLST_HD inline void fp2_tile_mul_by_u_plus_1(Fp2TileRegAvx512& out, const Fp2TileRegAvx512& a, const limb_t* m)
//    {
//        // Alias-safe: callers frequently use out == a in fp6/fp12 formulas.
//        FpTileRegAvx512 t0;
//        FpTileRegAvx512 t1;
//        fp_tile_sub_mod(t0, a[0], a[1], m);
//        fp_tile_add_mod(t1, a[0], a[1], m);
//        out[0] = t0;
//        out[1] = t1;
//    }
//
//    BLST_HD inline void fp2_tile_mul_mont(Fp2TileRegAvx512& out, const Fp2TileRegAvx512& a, const Fp2TileRegAvx512& b, const limb_t* m, limb_t n0)
//    {
//        FpTileRegAvx512 a0b0, a1b1, sum_a, sum_b, cross;
//        fp_tile_mul_mont(a0b0, a[0], b[0], m, n0);
//        fp_tile_mul_mont(a1b1, a[1], b[1], m, n0);
//        fp_tile_add_mod(sum_a, a[0], a[1], m);
//        fp_tile_add_mod(sum_b, b[0], b[1], m);
//        fp_tile_mul_mont(cross, sum_a, sum_b, m, n0);
//        fp_tile_sub_mod(out[0], a0b0, a1b1, m);
//        fp_tile_sub_mod(cross, cross, a0b0, m);
//        fp_tile_sub_mod(out[1], cross, a1b1, m);
//    }
//
//    BLST_HD inline void fp2_tile_sqr_mont(Fp2TileRegAvx512& out, const Fp2TileRegAvx512& a, const limb_t* m, limb_t n0)
//    {
//        FpTileRegAvx512 sum, diff, prod;
//        fp_tile_add_mod(sum, a[0], a[1], m);
//        fp_tile_sub_mod(diff, a[0], a[1], m);
//        fp_tile_mul_mont(out[0], sum, diff, m, n0);
//        fp_tile_mul_mont(prod, a[0], a[1], m, n0);
//        fp_tile_add_mod(out[1], prod, prod, m);
//    }
//
//    BLST_HD inline void fp2_tile_mem_zero(Fp2TileMemViewAvx512& mem) { fp_tile_mem_zero(mem[0]); fp_tile_mem_zero(mem[1]); }
//    BLST_HD inline void fp2_tile_reg_zero(Fp2TileRegAvx512& reg) { fp_tile_reg_zero(reg[0]); fp_tile_reg_zero(reg[1]); }
//    BLST_HD inline void fp2_tile_load(Fp2TileRegAvx512& reg, const Fp2TileMemViewAvx512& mem) { fp_tile_load(reg[0], mem[0]); fp_tile_load(reg[1], mem[1]); }
//    BLST_HD inline void fp2_tile_store(Fp2TileMemViewAvx512& mem, const Fp2TileRegAvx512& reg) { fp_tile_store(mem[0], reg[0]); fp_tile_store(mem[1], reg[1]); }
//    BLST_HD inline void fp4_tile_mem_zero(Fp4TileMemViewAvx512& mem) { fp2_tile_mem_zero(mem[0]); fp2_tile_mem_zero(mem[1]); }
//#endif
//} // namespace blst::tile
//
//#endif // BLST_BACKEND_AVX512_TILE_FP_OPS_H
