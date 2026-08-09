/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "blst/sign/aggregate.h"
#include "blst/curve/map_to_G.h"
#include "blst/curve/GT.h"

#include <new>
#include <type_traits>

namespace blst
{
    namespace {

    template<class Point, class Affine>
    BLST_HD void point_from_affine(Point& out, const Affine& in)
    {
        Point::from_affine(out, in);
    }

    template<class Affine, class Point>
    BLST_HD Affine point_to_affine(const Point& in)
    {
        Affine out = {};
        Point::to_affine(out, in);
        return out;
    }

    BLST_HD void miller_from_fp12_boundary(MillerResult& out, const Fp12& in)
    {
#if BLST_DEFAULT_AVX
        out.data = avx::avx_gt_detail::from_scalar(in);
#else
        out.data = in;
#endif
    }

    BLST_HD void fp12_from_miller_boundary(Fp12& out, const MillerResult& in)
    {
#if BLST_DEFAULT_AVX
        out = avx::avx_gt_detail::to_scalar(in.data);
#else
        out = in.data;
#endif
    }

    }  // namespace


    /*
     * PAIRING is a batched Miller-product state machine.
     *
     * Normal verification:
     *   1. Initialize with the hash/encode mode and borrowed DST.
     *   2. Add each public-key/message pair. The aggregate signature may be
     *      supplied once; later calls use a null signature.
     *   3. Commit to flush the final partial batch.
     *   4. Final-verify with a null external signature contribution.
     *
     * External signature aggregation:
     *   Aggregate encoded signatures with blst_aggregate_in_g1/g2, turn the
     *   final affine signature into a Miller contribution with
     *   blst_aggregated_in_g1/g2, aggregate all PK/message pairs with null
     *   signatures, commit, and pass the external contribution to finalverify.
     *
     * Parallel verification:
     *   Each worker owns a context, aggregates a slice, and commits it. Only
     *   committed contexts of the same PK/signature orientation may be merged.
     *   The merged context can then be final-verified normally.
     *
     * Expert raw pairing:
     *   raw_aggregate queues trusted on-curve affine (G2,G1) pairs without
     *   validation. get_miller_result commits and exports the accumulated
     *   pre-final-exponentiation value; final_exp is a separate operation.
     *
     * The default aggregate entry points validate the signature subgroup and
     * reject an infinite PK, but assume the PK was already subgroup-validated.
     * The chk_n variants make both subgroup checks explicit. All routines are
     * allocation-free and no exception may cross their C boundary.
     */


    enum {
        AGGR_UNDEFINED = 0,
        AGGR_MIN_SIG = 1,
        AGGR_MIN_PK = 2,
        AGGR_SIGN_SET = 0x10,
        AGGR_GT_SET = 0x20,
        AGGR_HASH_OR_ENCODE = 0x40,
        AGGR_DST_AFTER_CTX = 0x80
    };
    constexpr unsigned int kMinSigOrPk = AGGR_MIN_SIG | AGGR_MIN_PK;

    constexpr size_t sizeof_pairing = (sizeof(PAIRING) + 7) & ~(size_t)7;
    static_assert(sizeof_pairing == sizeof(PAIRING));

    size_t blst_pairing_sizeof(void)
    {
        return sizeof_pairing;
    }

    void blst_pairing_init(PAIRING* ctx, int hash_or_encode,
        const void* DST, size_t DST_len)
    {
        static_assert(std::is_trivially_default_constructible_v<MillerResult>);
        static_assert(std::is_trivially_default_constructible_v<AggregatedSignature>);
        static_assert(std::is_trivially_default_constructible_v<G2AffinePoint>);
        static_assert(std::is_trivially_default_constructible_v<G1AffinePoint>);
        static_assert(std::is_trivially_default_constructible_v<PAIRING>);
        static_assert(std::is_trivially_destructible_v<PAIRING>);
        ctx = ::new (static_cast<void*>(ctx)) PAIRING;
        const void* trailing_dst = reinterpret_cast<const byte*>(ctx) + sizeof_pairing;
        const bool dst_after_ctx = DST == trailing_dst;
        ctx->ctrl = AGGR_UNDEFINED | (hash_or_encode ? AGGR_HASH_OR_ENCODE : 0)
            | (dst_after_ctx ? AGGR_DST_AFTER_CTX : 0);
        ctx->nelems = 0;
        ctx->DST = dst_after_ctx ? nullptr : DST;
        ctx->DST_len = DST_len;
    }

    BLST_HD const void* pairing_get_dst(const PAIRING* ctx)
    {
        return (ctx->ctrl & AGGR_DST_AFTER_CTX) != 0
            ? reinterpret_cast<const byte*>(ctx) + sizeof_pairing
            : ctx->DST;
    }

    const void* blst_pairing_get_dst(const PAIRING* ctx)
    {
        return pairing_get_dst(ctx);
    }

    /*
     * Optional |nbits|-wide |scalar| is used to facilitate multiple aggregated
     * signature verification as discussed at
     * https://ethresear.ch/t/fast-verification-of-multiple-bls-signatures/5407.
     * Usage pattern is not finalized yet, because (sig != NULL) is better and
     * will be handled separately...
     */
    BLST_HD BLST_ERROR PAIRING_Aggregate_PK_in_G2(PAIRING* ctx,
        const G2AffinePoint* PK,
        size_t pk_groupcheck,
        const G1AffinePoint* sig,
        size_t sig_groupcheck,
        const byte* scalar, size_t nbits,
        const void* msg, size_t msg_len,
        const void* aug, size_t aug_len)
    {
		if (nbits > 256 || (nbits != 0 && scalar == NULL))
			return BLST_BAD_SCALAR;

        if (ctx->ctrl & AGGR_MIN_PK)
            return BLST_AGGR_TYPE_MISMATCH;

        const bool aggregate_signature = sig != NULL && !G1AffinePoint::is_inf(*sig);
        G1Point signature_point;
        if (aggregate_signature) {
            point_from_affine(signature_point, *sig);
            if (sig_groupcheck && !G1Point::in_subgroup(signature_point))
                return BLST_POINT_NOT_IN_GROUP;
        }

        if (PK != NULL) {
            if (G2AffinePoint::is_inf(*PK))
                return BLST_PK_IS_INFINITY;
            if (pk_groupcheck && !G2AffinePoint::in_subgroup(*PK))
                return BLST_POINT_NOT_IN_GROUP;
        }

        /* All fallible validation is complete before context publication. */
        ctx->ctrl |= AGGR_MIN_SIG;

        /*
         * Since we don't know if the signature is individual or aggregated,
         * the only sensible thing to do is to skip over infinite one and
         * count on the corresponding infinite public key to be rejected,
         * in case the signature is non-aggregated that is.
         */
        if (aggregate_signature) {
            G1Point& S = ctx->AggrSign.e1;

            if (ctx->ctrl & AGGR_SIGN_SET) {
                if (nbits != 0 && scalar != NULL) {
                    G1Point::mult_vartime(signature_point, signature_point, scalar, nbits);
                    G1Point::dadd(S, S, signature_point);
                }
                else {
                    G1Point::dadd(S, S, *sig);
                }
            }
            else {
                ctx->ctrl |= AGGR_SIGN_SET;
                if (nbits != 0 && scalar != NULL)
                    G1Point::mult_vartime(S, signature_point, scalar, nbits);
                else
                    S = signature_point;
            }
        }

        if (PK != NULL) {
            unsigned int n;
            G1Point H;
            const void* DST = pairing_get_dst(ctx);

            if (ctx->ctrl & AGGR_HASH_OR_ENCODE)
                Hash_to_curve(H,
                    (const unsigned char*)msg, msg_len,
                    (const unsigned char*)DST, ctx->DST_len,
                    (const unsigned char*)aug, aug_len);
            else
                Encode_to_curve(H,
                    (const unsigned char*)msg, msg_len,
                    (const unsigned char*)DST, ctx->DST_len,
                    (const unsigned char*)aug, aug_len);

            if (nbits != 0 && scalar != NULL)
                G1Point::mult_vartime(H, H, scalar, nbits);

            n = ctx->nelems;
            ctx->Q[n] = *PK;
            ctx->P[n] = point_to_affine<G1AffinePoint>(H);
            if (++n == kPairingBatchCapacity) {
                if (ctx->ctrl & AGGR_GT_SET) {
                    MillerResult gt;
                    GT::miller_loop_n(gt, ctx->Q, ctx->P, n);
                    MillerResult::mul(ctx->gt, ctx->gt, gt);
                }
                else {
                    GT::miller_loop_n(ctx->gt, ctx->Q, ctx->P, n);
                    ctx->ctrl |= AGGR_GT_SET;
                }
                n = 0;
            }
            ctx->nelems = n;
        }

        return BLST_SUCCESS;
    }

    BLST_HD BLST_ERROR blst_pairing_aggregate_pk_in_g2(PAIRING* ctx,
        const G2AffinePoint* PK,
        const G1AffinePoint* signature,
        const void* msg, size_t msg_len,
        const void* aug, size_t aug_len)
    {
        return PAIRING_Aggregate_PK_in_G2(ctx, PK, 0, signature, 1, NULL, 0,
            msg, msg_len, aug, aug_len);
    }

    BLST_HD BLST_ERROR blst_pairing_mul_n_aggregate_pk_in_g2(PAIRING* ctx,
        const G2AffinePoint* PK,
        const G1AffinePoint* sig,
        const byte* scalar,
        size_t nbits,
        const void* msg,
        size_t msg_len,
        const void* aug,
        size_t aug_len)
    {
        return PAIRING_Aggregate_PK_in_G2(ctx, PK, 0, sig, 1, scalar, nbits,
            msg, msg_len, aug, aug_len);
    }

    BLST_HD BLST_ERROR blst_pairing_chk_n_aggr_pk_in_g2(PAIRING* ctx,
        const G2AffinePoint* PK,
        size_t pk_grpchk,
        const G1AffinePoint* signature,
        size_t sig_grpchk,
        const void* msg, size_t msg_len,
        const void* aug, size_t aug_len)
    {
        return PAIRING_Aggregate_PK_in_G2(ctx, PK, pk_grpchk, signature, sig_grpchk,
            NULL, 0, msg, msg_len, aug, aug_len);
    }

    BLST_HD BLST_ERROR blst_pairing_chk_n_mul_n_aggr_pk_in_g2(PAIRING* ctx,
        const G2AffinePoint* PK,
        size_t pk_grpchk,
        const G1AffinePoint* sig,
        size_t sig_grpchk,
        const byte* scalar,
        size_t nbits,
        const void* msg,
        size_t msg_len,
        const void* aug,
        size_t aug_len)
    {
        return PAIRING_Aggregate_PK_in_G2(ctx, PK, pk_grpchk, sig, sig_grpchk,
            scalar, nbits,
            msg, msg_len, aug, aug_len);
    }

    BLST_HD BLST_ERROR PAIRING_Aggregate_PK_in_G1(PAIRING* ctx,
        const G1AffinePoint* PK,
        size_t pk_groupcheck,
        const G2AffinePoint* sig,
        size_t sig_groupcheck,
        const byte* scalar, size_t nbits,
        const void* msg, size_t msg_len,
        const void* aug, size_t aug_len)
    {
		if (nbits > 256 || (nbits != 0 && scalar == NULL))
			return BLST_BAD_SCALAR;

        if (ctx->ctrl & AGGR_MIN_SIG)
            return BLST_AGGR_TYPE_MISMATCH;

        const bool aggregate_signature = sig != NULL && !G2AffinePoint::is_inf(*sig);
        G2Point signature_point;
        if (aggregate_signature) {
            point_from_affine(signature_point, *sig);
            if (sig_groupcheck && !G2Point::in_subgroup(signature_point))
                return BLST_POINT_NOT_IN_GROUP;
        }

        if (PK != NULL) {
            if (G1AffinePoint::is_inf(*PK))
                return BLST_PK_IS_INFINITY;
            if (pk_groupcheck && !G1AffinePoint::in_subgroup(*PK))
                return BLST_POINT_NOT_IN_GROUP;
        }

        /* All fallible validation is complete before context publication. */
        ctx->ctrl |= AGGR_MIN_PK;

        /*
         * Since we don't know if the signature is individual or aggregated,
         * the only sensible thing to do is to skip over infinite one and
         * count on the corresponding infinite public key to be rejected,
         * in case the signature is non-aggregated that is.
         */
        if (aggregate_signature) {
            G2Point& S = ctx->AggrSign.e2;

            if (ctx->ctrl & AGGR_SIGN_SET) {
                if (nbits != 0 && scalar != NULL) {
                    G2Point::mult_vartime(signature_point, signature_point, scalar, nbits);
                    G2Point::dadd(S, S, signature_point);
                }
                else {
                    G2Point::dadd(S, S, *sig);
                }
            }
            else {
                ctx->ctrl |= AGGR_SIGN_SET;
                if (nbits != 0 && scalar != NULL)
                    G2Point::mult_vartime(S, signature_point, scalar, nbits);
                else
                    S = signature_point;
            }
        }

        if (PK != NULL) {
            unsigned int n;
            G2Point H;
            G1Point pk;
            const void* DST = pairing_get_dst(ctx);

            if (ctx->ctrl & AGGR_HASH_OR_ENCODE)
                Hash_to_curve(H,
                    (const unsigned char*)msg, msg_len,
                    (const unsigned char*)DST, ctx->DST_len,
                    (const unsigned char*)aug, aug_len);
            else
                Encode_to_curve(H,
                    (const unsigned char*)msg, msg_len,
                    (const unsigned char*)DST, ctx->DST_len,
                    (const unsigned char*)aug, aug_len);

            G1AffinePoint scaled_pk = {};
            if (nbits != 0 && scalar != NULL) {
                point_from_affine(pk, *PK);
                G1Point::mult_vartime(pk, pk, scalar, nbits);
                scaled_pk = point_to_affine<G1AffinePoint>(pk);
                PK = &scaled_pk;
            }

            n = ctx->nelems;
            ctx->Q[n] = point_to_affine<G2AffinePoint>(H);
            ctx->P[n] = *PK;
            if (++n == kPairingBatchCapacity) {
                if (ctx->ctrl & AGGR_GT_SET) {
                    MillerResult gt;
                    GT::miller_loop_n(gt, ctx->Q, ctx->P, n);
                    MillerResult::mul(ctx->gt, ctx->gt, gt);
                }
                else {
                    GT::miller_loop_n(ctx->gt, ctx->Q, ctx->P, n);
                    ctx->ctrl |= AGGR_GT_SET;
                }
                n = 0;
            }
            ctx->nelems = n;
        }

        return BLST_SUCCESS;
    }

    BLST_HD BLST_ERROR blst_pairing_aggregate_pk_in_g1(PAIRING* ctx,
        const G1AffinePoint* PK,
        const G2AffinePoint* signature,
        const void* msg, size_t msg_len,
        const void* aug, size_t aug_len)
    {
        return PAIRING_Aggregate_PK_in_G1(ctx, PK, 0, signature, 1, NULL, 0,
            msg, msg_len, aug, aug_len);
    }

    BLST_HD BLST_ERROR blst_pairing_mul_n_aggregate_pk_in_g1(PAIRING* ctx,
        const G1AffinePoint* PK,
        const G2AffinePoint* sig,
        const byte* scalar,
        size_t nbits,
        const void* msg,
        size_t msg_len,
        const void* aug,
        size_t aug_len)
    {
        return PAIRING_Aggregate_PK_in_G1(ctx, PK, 0, sig, 1, scalar, nbits,
            msg, msg_len, aug, aug_len);
    }

    BLST_HD BLST_ERROR blst_pairing_chk_n_aggr_pk_in_g1(PAIRING* ctx,
        const G1AffinePoint* PK,
        size_t pk_grpchk,
        const G2AffinePoint* signature,
        size_t sig_grpchk,
        const void* msg, size_t msg_len,
        const void* aug, size_t aug_len)
    {
        return PAIRING_Aggregate_PK_in_G1(ctx, PK, pk_grpchk, signature, sig_grpchk,
            NULL, 0, msg, msg_len, aug, aug_len);
    }

    BLST_HD BLST_ERROR blst_pairing_chk_n_mul_n_aggr_pk_in_g1(PAIRING* ctx,
        const G1AffinePoint* PK,
        size_t pk_grpchk,
        const G2AffinePoint* sig,
        size_t sig_grpchk,
        const byte* scalar,
        size_t nbits,
        const void* msg,
        size_t msg_len,
        const void* aug,
        size_t aug_len)
    {
        return PAIRING_Aggregate_PK_in_G1(ctx, PK, pk_grpchk, sig, sig_grpchk,
            scalar, nbits,
            msg, msg_len, aug, aug_len);
    }

    BLST_HD void PAIRING_Commit(PAIRING* ctx)
    {
        unsigned int n;

        if ((n = ctx->nelems) != 0) {
            if (ctx->ctrl & AGGR_GT_SET) {
                MillerResult gt;
                GT::miller_loop_n(gt, ctx->Q, ctx->P, n);
                MillerResult::mul(ctx->gt, ctx->gt, gt);
            }
            else {
                GT::miller_loop_n(ctx->gt, ctx->Q, ctx->P, n);
                ctx->ctrl |= AGGR_GT_SET;
            }
            ctx->nelems = 0;
        }
    }

    BLST_HD void blst_pairing_commit(PAIRING* ctx)
    {
        PAIRING_Commit(ctx);
    }

    BLST_HD BLST_ERROR blst_pairing_merge(PAIRING* ctx, const PAIRING* ctx1)
    {
        if ((ctx->ctrl & kMinSigOrPk) != AGGR_UNDEFINED
            && (ctx1->ctrl & kMinSigOrPk) != AGGR_UNDEFINED
            && (ctx->ctrl & ctx1->ctrl & kMinSigOrPk) == 0)
            return BLST_AGGR_TYPE_MISMATCH;

        /* context producers are expected to have called blst_pairing_commit */
        if (ctx->nelems || ctx1->nelems)
            return BLST_AGGR_TYPE_MISMATCH;

        ctx->ctrl |= ctx1->ctrl & kMinSigOrPk;

        switch (ctx->ctrl & kMinSigOrPk) {
        case AGGR_MIN_SIG:
            if (ctx->ctrl & ctx1->ctrl & AGGR_SIGN_SET) {
                G1Point::dadd(ctx->AggrSign.e1, ctx->AggrSign.e1,
                    ctx1->AggrSign.e1);
            }
            else if (ctx1->ctrl & AGGR_SIGN_SET) {
                ctx->ctrl |= AGGR_SIGN_SET;
                ctx->AggrSign.e1 = ctx1->AggrSign.e1;
            }
            break;
        case AGGR_MIN_PK:
            if (ctx->ctrl & ctx1->ctrl & AGGR_SIGN_SET) {
                G2Point::dadd(ctx->AggrSign.e2, ctx->AggrSign.e2,
                    ctx1->AggrSign.e2);
            }
            else if (ctx1->ctrl & AGGR_SIGN_SET) {
                ctx->ctrl |= AGGR_SIGN_SET;
                ctx->AggrSign.e2 = ctx1->AggrSign.e2;
            }
            break;
        case AGGR_UNDEFINED:
            break;
        default:
            return BLST_AGGR_TYPE_MISMATCH;
        }

        if (ctx->ctrl & ctx1->ctrl & AGGR_GT_SET) {
            MillerResult::mul(ctx->gt, ctx->gt, ctx1->gt);
        }
        else if (ctx1->ctrl & AGGR_GT_SET) {
            ctx->ctrl |= AGGR_GT_SET;
            ctx->gt = ctx1->gt;
        }

        return BLST_SUCCESS;
    }

    BLST_HD bool_t PAIRING_FinalVerify(const PAIRING* ctx, const Fp12* GTsig)
    {
        MillerResult gt;

        if (!(ctx->ctrl & AGGR_GT_SET))
            return 0;

        if (GTsig != NULL) {
            miller_from_fp12_boundary(gt, *GTsig);
        }
        else if (ctx->ctrl & AGGR_SIGN_SET) {
            switch (ctx->ctrl & kMinSigOrPk) {
            case AGGR_MIN_SIG:
            {
                const G2AffinePoint& g2 = G2AffinePoint::generator();
                const G1AffinePoint sig = point_to_affine<G1AffinePoint>(ctx->AggrSign.e1);
                GT::miller_loop_n(gt, &g2, &sig, 1);
                break;
            }
            case AGGR_MIN_PK:
            {
                const G1AffinePoint& g1 = G1AffinePoint::generator();
                const G2AffinePoint sig = point_to_affine<G2AffinePoint>(ctx->AggrSign.e2);
                GT::miller_loop_n(gt, &sig, &g1, 1);
                break;
            }
            default:
                return 0;
            }
        }
        else {
            /*
             * The aggregated signature was infinite, relation between the
             * hashes and the public keys has to be VERY special...
             */
            gt = MillerResult::identity();
        }

		MillerResult::conjugate(gt);
        MillerResult::mul(gt, gt, ctx->gt);
        GT result = {};
        GT::final_exp(result, gt);

        return GT::is_identity(result);
    }

    BLST_HD int blst_pairing_finalverify(const PAIRING* ctx, const Fp12* GTsig)
    {
        return (int)PAIRING_FinalVerify(ctx, GTsig);
    }

    BLST_HD int blst_fp12_finalverify(const MillerResult& GT1, const MillerResult& GT2)
    {
        MillerResult value = GT1;
        MillerResult::conjugate(value);
        MillerResult::mul(value, value, GT2);
        GT result = {};
        GT::final_exp(result, value);

        return (int)GT::is_identity(result);
    }

    BLST_HD void blst_pairing_raw_aggregate(PAIRING* ctx, const G2AffinePoint* q,
        const G1AffinePoint* p)
    {
        unsigned int n;

        if (G2AffinePoint::is_inf(*q) & G1AffinePoint::is_inf(*p))
            return;

        n = ctx->nelems;
        ctx->Q[n] = *q;
        ctx->P[n] = *p;
        if (++n == kPairingBatchCapacity) {
            if (ctx->ctrl & AGGR_GT_SET) {
                MillerResult gt;
                GT::miller_loop_n(gt, ctx->Q, ctx->P, n);
                MillerResult::mul(ctx->gt, ctx->gt, gt);
            }
            else {
                GT::miller_loop_n(ctx->gt, ctx->Q, ctx->P, n);
                ctx->ctrl |= AGGR_GT_SET;
            }
            n = 0;
        }
        ctx->nelems = n;
    }

    BLST_HD void blst_pairing_get_miller_result(Fp12* out, PAIRING* ctx)
    {
        PAIRING_Commit(ctx);
        if (ctx->ctrl & AGGR_GT_SET)
            fp12_from_miller_boundary(*out, ctx->gt);
        else
            *out = Fp12::one();
    }

    /*
     * PAIRING context-free entry points.
     *
     * To perform FastAggregateVerify, aggregate all public keys and
     * signatures with corresponding blst_aggregate_in_g{12}, convert
     * result to affine and call suitable blst_core_verify_pk_in_g{12}
     * or blst_aggregated_in_g{12}...
     */
    BLST_HD BLST_ERROR blst_aggregate_in_g1(G1Point& out, const G1Point* in,
        const unsigned char* zwire)
    {
        G1AffinePoint affine = {};
        G1Point P = {};
        BLST_ERROR ret;

        ret = G1AffinePoint::deserialize(affine, zwire);

        if (ret != BLST_SUCCESS)
            return ret;

        if (G1AffinePoint::is_inf(affine)) {
            out = in == NULL ? G1Point::identity() : *in;
            return BLST_SUCCESS;
        }

        G1Point::from_affine(P, affine);

        if (!G1Point::in_subgroup(P))
            return BLST_POINT_NOT_IN_GROUP;

        if (in == NULL)
            out = P;
        else
            G1Point::dadd(out, *in, affine);

        return BLST_SUCCESS;
    }

    BLST_HD BLST_ERROR blst_aggregate_in_g2(G2Point& out, const G2Point* in,
        const unsigned char* zwire)
    {
        G2AffinePoint affine = {};
        G2Point P = {};
        BLST_ERROR ret;

        ret = G2AffinePoint::deserialize(affine, zwire);

        if (ret != BLST_SUCCESS)
            return ret;

        if (G2AffinePoint::is_inf(affine)) {
            out = in == NULL ? G2Point::identity() : *in;
            return BLST_SUCCESS;
        }

        G2Point::from_affine(P, affine);

        if (!G2Point::in_subgroup(P))
            return BLST_POINT_NOT_IN_GROUP;

        if (in == NULL) {
            out = P;
        }
        else {
            G2Point::dadd(out, *in, affine);
        }
        return BLST_SUCCESS;
    }

    BLST_HD void blst_aggregated_in_g1(MillerResult& ret, const G1AffinePoint& sig)
    {
        const G2AffinePoint& g2 = G2AffinePoint::generator();
        GT::miller_loop_n(ret, &g2, &sig, 1);
    }

    BLST_HD void blst_aggregated_in_g2(MillerResult& ret, const G2AffinePoint& sig)
    {
        const G1AffinePoint& g1 = G1AffinePoint::generator();
        GT::miller_loop_n(ret, &sig, &g1, 1);
    }

    BLST_HD BLST_ERROR blst_core_verify_pk_in_g1(const G1AffinePoint* pk,
        const G2AffinePoint* signature,
        int hash_or_encode,
        const void* msg, size_t msg_len,
        const void* DST, size_t DST_len,
        const void* aug, size_t aug_len)
    {
        PAIRING ctx;
        BLST_ERROR ret;

        ctx.ctrl = AGGR_UNDEFINED | (hash_or_encode ? AGGR_HASH_OR_ENCODE : 0);
        ctx.nelems = 0;
        ctx.DST = DST;
        ctx.DST_len = DST_len;

        ret = PAIRING_Aggregate_PK_in_G1(&ctx, pk, 1, signature, 1, NULL, 0,
            msg, msg_len, aug, aug_len);
        if (ret != BLST_SUCCESS)
            return ret;

        PAIRING_Commit(&ctx);

        return PAIRING_FinalVerify(&ctx, NULL) ? BLST_SUCCESS : BLST_VERIFY_FAIL;
    }

    BLST_HD BLST_ERROR blst_core_verify_pk_in_g2(const G2AffinePoint* pk,
        const G1AffinePoint* signature,
        int hash_or_encode,
        const void* msg, size_t msg_len,
        const void* DST, size_t DST_len,
        const void* aug, size_t aug_len)
    {
        PAIRING ctx;
        BLST_ERROR ret;

        ctx.ctrl = AGGR_UNDEFINED | (hash_or_encode ? AGGR_HASH_OR_ENCODE : 0);
        ctx.nelems = 0;
        ctx.DST = DST;
        ctx.DST_len = DST_len;

        ret = PAIRING_Aggregate_PK_in_G2(&ctx, pk, 1, signature, 1, NULL, 0,
            msg, msg_len, aug, aug_len);
        if (ret != BLST_SUCCESS)
            return ret;

        PAIRING_Commit(&ctx);

        return PAIRING_FinalVerify(&ctx, NULL) ? BLST_SUCCESS : BLST_VERIFY_FAIL;
    }
}
