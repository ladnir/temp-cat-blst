#pragma once
#include "Fp.h" 
#include "Fp12.h" 
#include <cuda_runtime.h>

namespace blst {

    // 32 lanes * 16 bytes = 512 bytes. 
    // Aligned to 128 ensures every LimbTile starts on a cache sector boundary.
    struct alignas(128) LimbTile {
        ulong2 data[32];
    };

    struct FpTile {
        // 384-bit field = 6 limbs = 3x ulong2.
        // So we need 3 "rows" of data.
        LimbTile data[3];

        // --------------------------------------------------------
        // Converters: Global Memory <-> Local Registers
        // --------------------------------------------------------

        // Force inline to ensure the compiler registers-allocates 'ret' 
        // instead of spilling to local stack.
        BLST_DEVICE __forceinline__ Fp load() const {
            Fp ret;
            const int lane = threadIdx.x & 31; // Masking is cheaper than modulo

            // Read 128-bit chunks. 
            // These generate LD.E.128 (Global Load 128-bit)
            ulong2 l0 = data[0].data[lane];
            ulong2 l1 = data[1].data[lane];
            ulong2 l2 = data[2].data[lane];

            // Unpack ulong2 (x,y) into flat limbs [0..5]
            ret.data[0] = l0.x;
            ret.data[1] = l0.y;
            ret.data[2] = l1.x;
            ret.data[3] = l1.y;
            ret.data[4] = l2.x;
            ret.data[5] = l2.y;

            return ret;
        }

        BLST_DEVICE __forceinline__ void store(const Fp& in) {
            const int lane = threadIdx.x & 31;

            // Repack flat limbs into ulong2 for vector stores
            ulong2 l0 = make_ulong2(in.data[0], in.data[1]);
            ulong2 l1 = make_ulong2(in.data[2], in.data[3]);
            ulong2 l2 = make_ulong2(in.data[4], in.data[5]);

            // These generate ST.E.128 (Global Store 128-bit)
            data[0].data[lane] = l0;
            data[1].data[lane] = l1;
            data[2].data[lane] = l2;
        }

        BLST_DEVICE __forceinline__ auto& operator[](size_t idx) {
            return data[idx];
        }

        BLST_DEVICE __forceinline__ const auto& operator[](size_t idx) const {
            return data[idx];
        }

        // --------------------------------------------------------
        // Math Wrappers (Load -> Op -> Store)
        // --------------------------------------------------------
        // These are convenience wrappers. 
        // For heavy loops, manually call load(), do math on Fp, then store().

        BLST_DEVICE __forceinline__ static void add(FpTile& ret, const FpTile& a, const FpTile& b) {
            Fp la = a.load();
            Fp lb = b.load();
            Fp::add(la, la, lb); // Result in la
            ret.store(la);
        }

        BLST_DEVICE __forceinline__ static void sub(FpTile& ret, const FpTile& a, const FpTile& b) {
            Fp la = a.load();
            Fp lb = b.load();
            Fp::sub(la, la, lb);
            ret.store(la);
        }

        BLST_DEVICE __forceinline__ static void mul(FpTile& ret, const FpTile& a, const FpTile& b) {
            Fp la = a.load();
            Fp lb = b.load();
            Fp::mul(la, la, lb);
            ret.store(la);
        }

        BLST_DEVICE __forceinline__ static void sqr(FpTile& ret, const FpTile& a) {
            Fp la = a.load();
            Fp::sqr(la, la);
            ret.store(la);
        }

        // --------------------------------------------------------
        // Accessors
        // --------------------------------------------------------

        //// Direct access if you need to manually poke a specific limb
        //BLST_DEVICE __forceinline__ ulong2& operator()(size_t limb_idx_div_2) {
        //    return data[limb_idx_div_2].data[threadIdx.x & 31];
        //}

        //BLST_DEVICE __forceinline__ const ulong2& operator()(size_t limb_idx_div_2) const {
        //    return lidatambs[limb_idx_div_2].data[threadIdx.x & 31];
        //}
    };


    // ==================================================================================
// Fp2 Tile
// ==================================================================================
    struct Fp2Tile {
        FpTile data[2]; // c0, c1

        // --------------------------------------------------------
        // Load / Store
        // --------------------------------------------------------
        BLST_DEVICE __forceinline__ Fp2 load() const {
            Fp2 ret;
            // Assuming Fp2 has members like fp[0], fp[1] or similar.
            // Adjust .fp[i] to match your local Fp2 struct definition.
            ret[0] = data[0].load();
            ret[1] = data[1].load();
            return ret;
        }

        BLST_DEVICE __forceinline__ void store(const Fp2& in) {
            data[0].store(in[0]);
            data[1].store(in[1]);
        }
        BLST_DEVICE __forceinline__ auto& operator[](size_t idx) {
            return data[idx];
        }

        BLST_DEVICE __forceinline__ const auto& operator[](size_t idx) const {
            return data[idx];
        }

        // --------------------------------------------------------
        // Core Math Wrappers
        // --------------------------------------------------------
        BLST_DEVICE __forceinline__ static void add(Fp2Tile& ret, const Fp2Tile& a, const Fp2Tile& b) {
            Fp2 la = a.load();
            Fp2 lb = b.load();
            Fp2::add(la, la, lb);
            ret.store(la);
        }

        BLST_DEVICE __forceinline__ static void sub(Fp2Tile& ret, const Fp2Tile& a, const Fp2Tile& b) {
            Fp2 la = a.load();
            Fp2 lb = b.load();
            Fp2::sub(la, la, lb);
            ret.store(la);
        }

        BLST_DEVICE __forceinline__ static void mul(Fp2Tile& ret, const Fp2Tile& a, const Fp2Tile& b) {
            Fp2 la = a.load();
            Fp2 lb = b.load();
            Fp2::mul(la, la, lb);
            ret.store(la);
        }

        BLST_DEVICE __forceinline__ static void sqr(Fp2Tile& ret, const Fp2Tile& a) {
            Fp2 la = a.load();
            Fp2::sqr(la, la);
            ret.store(la);
        }

        BLST_DEVICE __forceinline__ static void neg(Fp2Tile& ret, const Fp2Tile& a) {
            Fp2 la = a.load();
            Fp2::neg(la, la);
            ret.store(la);
        }
    };

    // ==================================================================================
    // Fp4 Tile
    // Note: Fp4 is often just an array of 2 Fp2s in BLST towers, 
    // but provided here as a distinct type per your prototype.
    // ==================================================================================
    struct Fp4Tile {
        Fp2Tile data[2];

        BLST_DEVICE __forceinline__ Fp4 load() const {
            Fp4 ret;
            // Assuming Fp4 is struct { Fp2 data[2]; } or similar
            ret.data[0] = data[0].load();
            ret.data[1] = data[1].load();
            return ret;
        }

        BLST_DEVICE __forceinline__ void store(const Fp4& in) {
            data[0].store(in.data[0]);
            data[1].store(in.data[1]);
        }

        BLST_DEVICE __forceinline__ auto& operator[](size_t idx) {
            return data[idx];
        }

        BLST_DEVICE __forceinline__ const auto& operator[](size_t idx) const {
            return data[idx];
        }
        //// Most libraries don't expose explicit Fp4 math, but if you have it:
        //BLST_DEVICE __forceinline__ static void add(Fp4Tile& ret, const Fp4Tile& a, const Fp4Tile& b) {
        //    Fp4 la = a.load();
        //    Fp4 lb = b.load();
        //    Fp4::add(la, la, lb);
        //    ret.store(la);
        //}

        // ... sub, mul, sqr as needed
    };

    // ==================================================================================
    // Fp6 Tile
    // Structure: 3 Fp2 elements (usually for Sextic Twist)
    // ==================================================================================
    struct Fp6Tile {
        Fp2Tile data[3];

        BLST_DEVICE __forceinline__ Fp6 load() const {
            Fp6 ret;
            // Assuming Fp6 is struct { Fp2 data[3]; } or similar
            ret.data[0] = data[0].load();
            ret.data[1] = data[1].load();
            ret.data[2] = data[2].load();
            return ret;
        }

        BLST_DEVICE __forceinline__ void store(const Fp6& in) {
            data[0].store(in.data[0]);
            data[1].store(in.data[1]);
            data[2].store(in.data[2]);
        }

        BLST_DEVICE __forceinline__ auto& operator[](size_t idx) {
            return data[idx];
        }

        BLST_DEVICE __forceinline__ const auto& operator[](size_t idx) const {
            return data[idx];
        }
        //// --------------------------------------------------------
        //// Core Math Wrappers
        //// --------------------------------------------------------
        //BLST_DEVICE __forceinline__ static void add(Fp6Tile& ret, const Fp6Tile& a, const Fp6Tile& b) {
        //    Fp6 la = a.load();
        //    Fp6 lb = b.load();
        //    Fp6::add(la, la, lb);
        //    ret.store(la);
        //}

        //BLST_DEVICE __forceinline__ static void sub(Fp6Tile& ret, const Fp6Tile& a, const Fp6Tile& b) {
        //    Fp6 la = a.load();
        //    Fp6 lb = b.load();
        //    Fp6::sub(la, la, lb);
        //    ret.store(la);
        //}

        //BLST_DEVICE __forceinline__ static void mul(Fp6Tile& ret, const Fp6Tile& a, const Fp6Tile& b) {
        //    Fp6 la = a.load();
        //    Fp6 lb = b.load();
        //    Fp6::mul(la, la, lb);
        //    ret.store(la);
        //}

        //BLST_DEVICE __forceinline__ static void sqr(Fp6Tile& ret, const Fp6Tile& a) {
        //    Fp6 la = a.load();
        //    Fp6::sqr(la, la);
        //    ret.store(la);
        //}
    };

    // ==================================================================================
    // Fp12 Tile
    // Structure: 2 Fp6 elements
    // ==================================================================================
    struct Fp12Tile {
        Fp6Tile data[2];

        BLST_DEVICE __forceinline__ Fp12 load() const {
            Fp12 ret;
            // Assuming Fp12 is struct { Fp6 data[2]; }
            ret.data[0] = data[0].load();
            ret.data[1] = data[1].load();
            return ret;
        }

        BLST_DEVICE __forceinline__ void store(const Fp12& in) {
            data[0].store(in.data[0]);
            data[1].store(in.data[1]);
        }

        BLST_DEVICE __forceinline__ auto& operator[](size_t idx) {
            return data[idx];
		}

        BLST_DEVICE __forceinline__ const auto& operator[](size_t idx) const {
			return data[idx];
		}

        // --------------------------------------------------------
        // Core Math Wrappers
        // --------------------------------------------------------
        BLST_DEVICE __forceinline__ static void add(Fp12Tile& ret, const Fp12Tile& a, const Fp12Tile& b) {
            Fp12 la = a.load();
            Fp12 lb = b.load();
            Fp12::add(la, la, lb);
            ret.store(la);
        }

        BLST_DEVICE __forceinline__ static void sub(Fp12Tile& ret, const Fp12Tile& a, const Fp12Tile& b) {
            Fp12 la = a.load();
            Fp12 lb = b.load();
            Fp12::sub(la, la, lb);
            ret.store(la);
        }

        BLST_DEVICE __forceinline__ static void mul(Fp12Tile& ret, const Fp12Tile& a, const Fp12Tile& b) {
            Fp12 la = a.load();
            Fp12 lb = b.load();
            Fp12::mul(la, la, lb);
            ret.store(la);
        }

        BLST_DEVICE __forceinline__ static void sqr(Fp12Tile& ret, const Fp12Tile& a) {
            Fp12 la = a.load();
            Fp12::sqr(la, la);
            ret.store(la);
        }

        // Cyclotomic squaring is critical for pairings (Miller loop / Final Exp)
        BLST_DEVICE __forceinline__ static void cyclotomic_sqr(Fp12Tile& ret, const Fp12Tile& a) {
            Fp12 la = a.load();
            // Assuming your Fp12 class has this specialized function
            Fp12::cyclotomic_sqr(la, la);
            ret.store(la);
        }

        BLST_DEVICE __forceinline__ static void inverse(Fp12Tile& ret, const Fp12Tile& a) {
            Fp12 la = a.load();
            Fp12::inverse(la, la);
            ret.store(la);
        }
    };
} // namespace blst

