#pragma once

#include "blst/curve/point.h"
#include "blst/curve/GT.h"

namespace blst
{
    /**
     * Number of affine pairing terms accumulated before evaluating one Miller
     * batch. This is an implementation constant, not a configurable ABI knob.
     */
    inline constexpr size_t kPairingBatchCapacity = 8;

    struct AggregatedSignature {
        G1Point e1;
        G2Point e2;
    };

    /**
     * Native storage backing the opaque C `blst_pairing` context.
     *
     * Callers allocate this storage through the public size/alignment queries
     * and manipulate it only through the pairing API. The state machine is:
     * initialize, aggregate zero or more same-orientation items, commit, then
     * final-verify, export the Miller result, or merge with another committed
     * context. Aggregation may resume after a commit.
     *
     * `gt` is a selected-backend Miller result, not a GT element. `Q` and `P`
     * hold the final partial batch. The signature slots are independent real
     * objects rather than union alternatives so their C++ lifetimes are
     * unambiguous. No field is a supported direct-call interface.
     */
    struct PAIRING {
        unsigned int ctrl;
        unsigned int nelems;
        const void* DST;
        size_t DST_len;
        MillerResult gt;
        AggregatedSignature AggrSign;
        G2AffinePoint Q[kPairingBatchCapacity];
        G1AffinePoint P[kPairingBatchCapacity];
    };

    BLST_HD int blst_fp12_finalverify(const MillerResult& GT1, const MillerResult& GT2);
    /** Commits pending pairs and exports the pre-final-exponentiation Miller product. */
    BLST_HD void blst_pairing_get_miller_result(Fp12* out, PAIRING* ctx);

}
