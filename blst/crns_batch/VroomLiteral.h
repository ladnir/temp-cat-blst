#ifndef BLST_CRNS_BATCH_VROOM_LITERAL_H
#define BLST_CRNS_BATCH_VROOM_LITERAL_H

#include "blst/fields/Fp.h"
#include "blst/fields/Fp12.h"

#include <array>

#if defined(BLST_CRNS_BATCH_VROOM_VENDOR)

#include "blst/crns_batch/vroom_port/src/bounded_ring.hpp"
#include "blst/crns_batch/vroom_port/src/fp12.hpp"

namespace blst::crns_batch::vroom_literal
{
    using VroomRing = BoundedRing<>;
    using VroomFpElement = typename VroomRing::StandardElement;
    using VroomFp12Ring = FP12Ring<VroomRing>;
    using VroomFp12Element = typename VroomFp12Ring::StandardElement;

    inline constexpr const char* kBls12381ModulusDec =
        "4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787";

    struct VroomFpPack8
    {
        std::array<VroomFpElement, 8> lanes{};
    };

    inline auto ring() -> const VroomRing&
    {
        static const VroomRing kRing(BigInt(kBls12381ModulusDec, 10));
        return kRing;
    }

    inline auto fp12_ring() -> const VroomFp12Ring&
    {
        static const VroomFp12Ring kFp12Ring(ring());
        return kFp12Ring;
    }

    inline auto fp_to_bigint(const Fp& in) -> BigInt
    {
        Fp canonical{};
        Fp::from_mont(canonical, in);

        BigInt out(0);
        for (int i = static_cast<int>(Fp::num_limb()) - 1; i >= 0; --i) {
            out <<= 64;
            out += BigInt(static_cast<unsigned long>(canonical[static_cast<size_t>(i)]));
        }
        return out;
    }

    inline auto bigint_to_fp(const BigInt& in) -> Fp
    {
        BigInt tmp = in;
        const BigInt mask = (BigInt(1) << 64) - BigInt(1);

        Fp canonical{};
        for (size_t i = 0; i < Fp::num_limb(); ++i) {
            canonical[i] = static_cast<u64>((tmp & mask).to_ulong());
            tmp >>= 64;
        }

        Fp out{};
        Fp::to_mont(out, canonical);
        return out;
    }

    inline auto fp12_coeff(const Fp12& in, size_t idx) -> const Fp&
    {
        const size_t fp6_idx = idx / 6;
        const size_t fp2_idx = (idx % 6) / 2;
        const size_t fp_idx = idx % 2;
        return in[fp6_idx][fp2_idx][fp_idx];
    }

    inline auto fp12_coeff(Fp12& in, size_t idx) -> Fp&
    {
        const size_t fp6_idx = idx / 6;
        const size_t fp2_idx = (idx % 6) / 2;
        const size_t fp_idx = idx % 2;
        return in[fp6_idx][fp2_idx][fp_idx];
    }

    inline auto fp12_to_bigint_array(const Fp12& in) -> std::array<BigInt, 12>
    {
        std::array<BigInt, 12> out{};
        for (size_t i = 0; i < out.size(); ++i) {
            out[i] = fp_to_bigint(fp12_coeff(in, i));
        }
        return out;
    }

    inline auto bigint_array_to_fp12(const std::array<BigInt, 12>& in) -> Fp12
    {
        Fp12 out{};
        for (size_t i = 0; i < in.size(); ++i) {
            fp12_coeff(out, i) = bigint_to_fp(in[i]);
        }
        return out;
    }

    inline void from_aos(VroomFpPack8& out, const Fp* in8)
    {
        const VroomRing& r = ring();
        for (size_t i = 0; i < out.lanes.size(); ++i) {
            out.lanes[i] = r.from_bigint(fp_to_bigint(in8[i]));
        }
    }

    inline void to_aos(Fp* out8, const VroomFpPack8& in)
    {
        const VroomRing& r = ring();
        for (size_t i = 0; i < in.lanes.size(); ++i) {
            out8[i] = bigint_to_fp(r.to_bigint(in.lanes[i]));
        }
    }

    inline void mul(VroomFpPack8& out, const VroomFpPack8& a, const VroomFpPack8& b)
    {
        out.lanes = ring().batch_modmul<8>(a.lanes, b.lanes);
    }

    inline auto mul(const VroomFpPack8& a, const VroomFpPack8& b) -> VroomFpPack8
    {
        VroomFpPack8 out{};
        out.lanes = ring().batch_modmul<8>(a.lanes, b.lanes);
        return out;
    }

    inline auto fp12_from_blst(const Fp12& in) -> VroomFp12Element
    {
        VroomFp12Element out{};
        const auto coeffs = fp12_to_bigint_array(in);
        for (size_t i = 0; i < coeffs.size(); ++i) {
            out[i] = ring().from_bigint(coeffs[i]);
        }
        return out;
    }

    inline auto fp12_to_blst(const VroomFp12Element& in) -> Fp12
    {
        std::array<BigInt, 12> coeffs{};
        for (size_t i = 0; i < coeffs.size(); ++i) {
            coeffs[i] = ring().to_bigint(in[i]);
        }
        return bigint_array_to_fp12(coeffs);
    }

    inline auto fp12_mul(const VroomFp12Element& a, const VroomFp12Element& b) -> VroomFp12Element
    {
        return fp12_ring().fp12_flat_mul(a, b, ring());
    }
}

#endif

#endif
