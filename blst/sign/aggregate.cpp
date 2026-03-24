/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "blst/sign/aggregate.h"
#include "blst/curve/map_to_G.h"
#include "blst/curve/GT.h"
//#include "blst/curve/map_to_G.h"

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

    BLST_HD void gt_from_fp12_boundary(GT& out, const Fp12& in)
    {
        unsigned char serialized[12 * 48] = {};
        scalar::GT::serialize(serialized, scalar::GT(in));
        GT::deserialize(out, serialized);
    }

    BLST_HD void fp12_from_gt_boundary(Fp12& out, const GT& in)
    {
        unsigned char serialized[12 * 48] = {};
        GT::serialize(serialized, in);
        scalar::GT scalar_gt = {};
        scalar::GT::deserialize(scalar_gt, serialized);
        out = scalar_gt.data;
    }

    }  // namespace


    /*
     * Usage pattern on single-processor system is
     *
     * blst_pairing_init(ctx, hash_or_encode, DST);
     * blst_pairing_aggregate_pk_in_g1(ctx, PK[0], aggregated_signature, msg[0]);
     * blst_pairing_aggregate_pk_in_g1(ctx, PK[1], NULL, msg[1]);
     * ...
     * blst_pairing_commit(ctx);
     * blst_pairing_finalverify(ctx, NULL);
     *
     ***********************************************************************
     * Usage pattern on multi-processor system is
     *
     *   blst_pairing_init(pk[0], hash_or_encode, DST);
     *   blst_pairing_init(pk[1], hash_or_encode, DST);
     *   ...
     * start threads each processing an N/nthreads slice of PKs and messages:
     *     blst_pairing_aggregate_pk_in_g1(pk[i], PK[i*n+0], NULL, msg[i*n+0]);
     *     blst_pairing_aggregate_pk_in_g1(pk[i], PK[i*n+1], NULL, msg[i*n+1]);
     *     ...
     *     blst_pairing_commit(pkx);
     *   ...
     * meanwhile in main thread
     *   blst_fp12 gtsig;
     *   blst_aggregated_in_g2(&gtsig, aggregated_signature);
     * join threads and merge their contexts:
     *   blst_pairing_merge(pk[0], pk[1]);
     *   blst_pairing_merge(pk[0], pk[2]);
     *   ...
     *   blst_pairing_finalverify(pk[0], gtsig);
     */


    enum {
        AGGR_UNDEFINED = 0,
        AGGR_MIN_SIG = 1,
        AGGR_MIN_PK = 2,
        AGGR_SIGN_SET = 0x10,
        AGGR_GT_SET = 0x20,
        AGGR_HASH_OR_ENCODE = 0x40
    };
#define MIN_SIG_OR_PK (AGGR_MIN_SIG | AGGR_MIN_PK)

    const size_t sizeof_pairing = (sizeof(PAIRING) + 7) & ~(size_t)7;

    size_t blst_pairing_sizeof(void)
    {
        return sizeof_pairing;
    }

    void blst_pairing_init(PAIRING* ctx, int hash_or_encode,
        const void* DST, size_t DST_len)
    {
        ctx->ctrl = AGGR_UNDEFINED | (hash_or_encode ? AGGR_HASH_OR_ENCODE : 0);
        ctx->nelems = 0;
        ctx->DST = (uptr_t)DST == (uptr_t)((byte*)ctx + sizeof_pairing) ? (void*)42
            : DST;
        ctx->DST_len = DST_len;
    }

    BLST_HD const void* pairing_get_dst(const PAIRING* ctx)
    {
        return (uptr_t)ctx->DST == (uptr_t)42 ? (const byte*)ctx + sizeof_pairing
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
        if (ctx->ctrl & AGGR_MIN_PK)
            return BLST_AGGR_TYPE_MISMATCH;

        ctx->ctrl |= AGGR_MIN_SIG;

        /*
         * Since we don't know if the signature is individual or aggregated,
         * the only sensible thing to do is to skip over infinite one and
         * count on the corresponding infinite public key to be rejected,
         * in case the signature is non-aggregated that is.
         */
        if (sig != NULL && !G1AffinePoint::is_inf(*sig)) {
            G1Point& S = ctx->AggrSign.e1;
            G1Point P;

            point_from_affine(P, *sig);

            if (sig_groupcheck && !G1Point::in_subgroup(P))
                return BLST_POINT_NOT_IN_GROUP;

            if (ctx->ctrl & AGGR_SIGN_SET) {
                if (nbits != 0 && scalar != NULL) {
                    G1Point::mult(P, P, scalar, nbits);
                    G1Point::dadd(S, S, P);
                }
                else {
                    G1Point::dadd(S, S, *sig);
                }
            }
            else {
                ctx->ctrl |= AGGR_SIGN_SET;
                if (nbits != 0 && scalar != NULL)
                    G1Point::mult(S, P, scalar, nbits);
                else
                    S = P;
            }
        }

        if (PK != NULL) {
            unsigned int n;
            G1Point H;
            const void* DST = pairing_get_dst(ctx);

            /*
             * Reject infinite public keys.
             */
            if (G2AffinePoint::is_inf(*PK))
                return BLST_PK_IS_INFINITY;

            if (pk_groupcheck) {
                G2Point P;

                point_from_affine(P, *PK);
                if (!G2Point::in_subgroup(P))
                    return BLST_POINT_NOT_IN_GROUP;
            }

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
                G1Point::mult(H, H, scalar, nbits);

            n = ctx->nelems;
            ctx->Q[n] = *PK;
            ctx->P[n] = point_to_affine<G1AffinePoint>(H);
            if (++n == N_MAX) {
                if (ctx->ctrl & AGGR_GT_SET) {
                    GT gt;
                    GT::miller_loop_n(gt, ctx->Q, ctx->P, n);
                    GT::add(ctx->gt, ctx->gt, gt);
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
        if (ctx->ctrl & AGGR_MIN_SIG)
            return BLST_AGGR_TYPE_MISMATCH;

        ctx->ctrl |= AGGR_MIN_PK;

        /*
         * Since we don't know if the signature is individual or aggregated,
         * the only sensible thing to do is to skip over infinite one and
         * count on the corresponding infinite public key to be rejected,
         * in case the signature is non-aggregated that is.
         */
        if (sig != NULL && !G2AffinePoint::is_inf(*sig)) {
            G2Point& S = ctx->AggrSign.e2;
            G2Point P;

            point_from_affine(P, *sig);

            if (sig_groupcheck && !G2Point::in_subgroup(P))
                return BLST_POINT_NOT_IN_GROUP;

            if (ctx->ctrl & AGGR_SIGN_SET) {
                if (nbits != 0 && scalar != NULL) {
                    G2Point::mult(P, P, scalar, nbits);
                    G2Point::dadd(S, S, P);
                }
                else {
                    G2Point::dadd(S, S, *sig);
                }
            }
            else {
                ctx->ctrl |= AGGR_SIGN_SET;
                if (nbits != 0 && scalar != NULL)
                    G2Point::mult(S, P, scalar, nbits);
                else
                    S = P;
            }
        }

        if (PK != NULL) {
            unsigned int n;
            G2Point H;
            G1Point pk;
            const void* DST = pairing_get_dst(ctx);

            /*
             * Reject infinite public keys.
             */
            if (G1AffinePoint::is_inf(*PK))
                return BLST_PK_IS_INFINITY;

            if (pk_groupcheck) {
                G1Point P;

                point_from_affine(P, *PK);
                if (!G1Point::in_subgroup(P))
                    return BLST_POINT_NOT_IN_GROUP;
            }

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
                G1Point::mult(pk, pk, scalar, nbits);
                scaled_pk = point_to_affine<G1AffinePoint>(pk);
                PK = &scaled_pk;
            }

            n = ctx->nelems;
            ctx->Q[n] = point_to_affine<G2AffinePoint>(H);
            ctx->P[n] = *PK;
            if (++n == N_MAX) {
                if (ctx->ctrl & AGGR_GT_SET) {
                    GT gt;
                    GT::miller_loop_n(gt, ctx->Q, ctx->P, n);
                    GT::add(ctx->gt, ctx->gt, gt);
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
                GT gt;
                GT::miller_loop_n(gt, ctx->Q, ctx->P, n);
                GT::add(ctx->gt, ctx->gt, gt);
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
        if ((ctx->ctrl & MIN_SIG_OR_PK) != AGGR_UNDEFINED
            && (ctx1->ctrl & MIN_SIG_OR_PK) != AGGR_UNDEFINED
            && (ctx->ctrl & ctx1->ctrl & MIN_SIG_OR_PK) == 0)
            return BLST_AGGR_TYPE_MISMATCH;

        /* context producers are expected to have called blst_pairing_commit */
        if (ctx->nelems || ctx1->nelems)
            return BLST_AGGR_TYPE_MISMATCH;

        ctx->ctrl |= ctx1->ctrl & MIN_SIG_OR_PK;

        switch (ctx->ctrl & MIN_SIG_OR_PK) {
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
            GT::add(ctx->gt, ctx->gt, ctx1->gt);
        }
        else if (ctx1->ctrl & AGGR_GT_SET) {
            ctx->ctrl |= AGGR_GT_SET;
            ctx->gt = ctx1->gt;
        }

        return BLST_SUCCESS;
    }

    BLST_HD bool_t PAIRING_FinalVerify(const PAIRING* ctx, const Fp12* GTsig)
    {
        GT gt;

        if (!(ctx->ctrl & AGGR_GT_SET))
            return 0;

        if (GTsig != NULL) {
            gt_from_fp12_boundary(gt, *GTsig);
        }
        else if (ctx->ctrl & AGGR_SIGN_SET) {
            switch (ctx->ctrl & MIN_SIG_OR_PK) {
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
            gt = GT::identity();
        }

		GT::conjugate(gt);
        GT::add(gt, gt, ctx->gt);
        GT::final_mult(gt, gt);

        return GT::is_identity(gt);
    }

    BLST_HD int blst_pairing_finalverify(const PAIRING* ctx, const Fp12* GTsig)
    {
        return (int)PAIRING_FinalVerify(ctx, GTsig);
    }

    BLST_HD int blst_fp12_finalverify(const GT& GT1, const GT& GT2)
    {
        GT gt = GT1;
        GT::conjugate(gt);
        GT::add(gt, gt, GT2);
        GT::final_mult(gt, gt);

        return (int)GT::is_identity(gt);
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
        if (++n == N_MAX) {
            if (ctx->ctrl & AGGR_GT_SET) {
                GT gt;
                GT::miller_loop_n(gt, ctx->Q, ctx->P, n);
                GT::add(ctx->gt, ctx->gt, gt);
            }
            else {
                GT::miller_loop_n(ctx->gt, ctx->Q, ctx->P, n);
                ctx->ctrl |= AGGR_GT_SET;
            }
            n = 0;
        }
        ctx->nelems = n;
    }

    BLST_HD Fp12 blst_pairing_as_fp12(PAIRING* ctx)
    {
        PAIRING_Commit(ctx);
        Fp12 out = {};
        fp12_from_gt_boundary(out, ctx->gt);
        return out;
    }

    BLST_HD void blst_pairing_as_fp12_to(Fp12* out, PAIRING* ctx)
    {
        *out = blst_pairing_as_fp12(ctx);
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
            if (in == NULL)
                out = G1Point::identity();
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
            if (in == NULL)
                out = G2Point::identity();
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

    BLST_HD void blst_aggregated_in_g1(GT& ret, const G1AffinePoint& sig)
    {
        const G2AffinePoint& g2 = G2AffinePoint::generator();
        GT::miller_loop_n(ret, &g2, &sig, 1);
    }

    BLST_HD void blst_aggregated_in_g2(GT& ret, const G2AffinePoint& sig)
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
