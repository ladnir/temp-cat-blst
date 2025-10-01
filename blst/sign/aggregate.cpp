/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#include "blst/sign/aggregate.h"
#include "blst/curve/map_to_G.h"
//#include "blst/curve/map_to_G.h"

namespace blst
{

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

    BLST_HD void FROM_AFFINE(auto& out, const auto& in) 
    { 
        vec_copy(&out.X, &in.X, 2*sizeof(in.X)), 
        vec_select(&out.Z, &in.X, &get_BLST_Rx().p, sizeof(in.X), 
                             vec_is_zero(&in.X, 2*sizeof(in.X))); 
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
        if (sig != NULL && !vec_is_zero(sig, sizeof(*sig))) {
            G1Point& S = ctx->AggrSign.e1;
            G1Point P;

            FROM_AFFINE(P, *sig);

            if (sig_groupcheck && !G1Point::in_subgroup(P))
                return BLST_POINT_NOT_IN_GROUP;

            if (ctx->ctrl & AGGR_SIGN_SET) {
                if (nbits != 0 && scalar != NULL) {
                    mult_w<G1Point,5>(P, P, scalar, nbits);
                    G1Point::dadd(S, S, P, NULL);
                }
                else {
                    G1Point::dadd_affine(S, S, *sig);
                }
            }
            else {
                ctx->ctrl |= AGGR_SIGN_SET;
                if (nbits != 0 && scalar != NULL)
                    mult_w<G1Point,5>(S, P, scalar, nbits);
                else
                    vec_copy(&S, &P, sizeof(P));
            }
        }

        if (PK != NULL) {
            unsigned int n;
            G1Point H;
            const void* DST = pairing_get_dst(ctx);

            /*
             * Reject infinite public keys.
             */
            if (vec_is_zero(PK, sizeof(*PK)))
                return BLST_PK_IS_INFINITY;

            if (pk_groupcheck) {
                G2Point P;

                FROM_AFFINE(P, *PK);
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
                mult_w<G1Point,5>(H, H, scalar, nbits);

            G1Point::from_jacobian(H, H);

            n = ctx->nelems;
            vec_copy(ctx->Q + n, PK, sizeof(G2AffinePoint));
            vec_copy(ctx->P + n, &H, sizeof(G1AffinePoint));
            if (++n == N_MAX) {
                if (ctx->ctrl & AGGR_GT_SET) {
                    Fp12 GT;
                    miller_loop_n(GT, ctx->Q, ctx->P, n);
                    Fp12::mul(ctx->GT, ctx->GT, GT);
                }
                else {
                    miller_loop_n(ctx->GT, ctx->Q, ctx->P, n);
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
        if (sig != NULL && !vec_is_zero(sig, sizeof(*sig))) {
            G2Point& S = ctx->AggrSign.e2;
            G2Point P;

            FROM_AFFINE(P, *sig);

            if (sig_groupcheck && !G2Point::in_subgroup(P))
                return BLST_POINT_NOT_IN_GROUP;

            if (ctx->ctrl & AGGR_SIGN_SET) {
                if (nbits != 0 && scalar != NULL) {

                    mult_w< G2Point,5>(P, P, scalar, nbits);
                    G2Point::dadd(S, S, P, NULL);
                }
                else {
                    G2Point::dadd_affine(S, S, *sig);
                }
            }
            else {
                ctx->ctrl |= AGGR_SIGN_SET;
                if (nbits != 0 && scalar != NULL)
                    mult_w<G2Point,5>(S, P, scalar, nbits);
                else
                    vec_copy(&S, &P, sizeof(P));
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
            if (vec_is_zero(PK, sizeof(*PK)))
                return BLST_PK_IS_INFINITY;

            if (pk_groupcheck) {
                G1Point P;

                FROM_AFFINE(P, *PK);
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

            G2Point::from_jacobian(H, H);

            if (nbits != 0 && scalar != NULL) {
                FROM_AFFINE(pk, *PK);
                mult_w<G1Point,5>(pk, pk, scalar, nbits);
                G1Point::from_jacobian(pk, pk);
                PK = (const G1AffinePoint*)&pk;
            }

            n = ctx->nelems;
            vec_copy(ctx->Q + n, &H, sizeof(G2AffinePoint));
            vec_copy(ctx->P + n, &PK, sizeof(G1AffinePoint));
            if (++n == N_MAX) {
                if (ctx->ctrl & AGGR_GT_SET) {
                    Fp12 GT;
                    miller_loop_n(GT, ctx->Q, ctx->P, n);
                    Fp12::mul(ctx->GT, ctx->GT, GT);
                }
                else {
                    miller_loop_n(ctx->GT, ctx->Q, ctx->P, n);
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
                Fp12 GT;
                miller_loop_n(GT, ctx->Q, ctx->P, n);
                Fp12::mul(ctx->GT, ctx->GT, GT);
            }
            else {
                miller_loop_n(ctx->GT, ctx->Q, ctx->P, n);
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
                    ctx1->AggrSign.e1, NULL);
            }
            else if (ctx1->ctrl & AGGR_SIGN_SET) {
                ctx->ctrl |= AGGR_SIGN_SET;
                vec_copy(&ctx->AggrSign.e1, &ctx1->AggrSign.e1,
                    sizeof(ctx->AggrSign.e1));
            }
            break;
        case AGGR_MIN_PK:
            if (ctx->ctrl & ctx1->ctrl & AGGR_SIGN_SET) {
                G2Point::dadd(ctx->AggrSign.e2, ctx->AggrSign.e2,
                    ctx1->AggrSign.e2, NULL);
            }
            else if (ctx1->ctrl & AGGR_SIGN_SET) {
                ctx->ctrl |= AGGR_SIGN_SET;
                vec_copy(&ctx->AggrSign.e2, &ctx1->AggrSign.e2,
                    sizeof(ctx->AggrSign.e2));
            }
            break;
        case AGGR_UNDEFINED:
            break;
        default:
            return BLST_AGGR_TYPE_MISMATCH;
        }

        if (ctx->ctrl & ctx1->ctrl & AGGR_GT_SET) {
            Fp12::mul(ctx->GT, ctx->GT, ctx1->GT);
        }
        else if (ctx1->ctrl & AGGR_GT_SET) {
            ctx->ctrl |= AGGR_GT_SET;
            vec_copy(&ctx->GT, &ctx1->GT, sizeof(ctx->GT));
        }

        return BLST_SUCCESS;
    }

    BLST_HD bool_t PAIRING_FinalVerify(const PAIRING* ctx, const Fp12* GTsig)
    {
        Fp12 GT;

        if (!(ctx->ctrl & AGGR_GT_SET))
            return 0;

        if (GTsig != NULL) {
            vec_copy(&GT, GTsig, sizeof(GT));
        }
        else if (ctx->ctrl & AGGR_SIGN_SET) {
            AggregatedSignature AggrSign;

            switch (ctx->ctrl & MIN_SIG_OR_PK) {
            case AGGR_MIN_SIG:
                G1Point::from_jacobian(AggrSign.e1, ctx->AggrSign.e1);
                miller_loop_n(GT, (const G2AffinePoint*)&get_BLST_G2(),
                    (const G1AffinePoint*)&AggrSign.e1, 1);
                break;
            case AGGR_MIN_PK:
                G2Point::from_jacobian(AggrSign.e2, ctx->AggrSign.e2);
                miller_loop_n(GT, (const G2AffinePoint*)&AggrSign.e2,
                    (const G1AffinePoint*)&get_BLST_G1(), 1);
                break;
            default:
                return 0;
            }
        }
        else {
            /*
             * The aggregated signature was infinite, relation between the
             * hashes and the public keys has to be VERY special...
             */
            vec_copy(&GT, get_BLST_Rx().p12, sizeof(GT));
        }

        Fp12::conjugate(GT);
        Fp12::mul(GT, GT, ctx->GT);
        Fp12::final_exp(GT, GT);

        /* return GT==1 */
        return vec_is_equal(&GT[0][0], get_BLST_Rx().p2, sizeof(GT[0][0])) &
            vec_is_zero(&GT[0][1], sizeof(GT) - sizeof(GT[0][0]));
    }

    BLST_HD int blst_pairing_finalverify(const PAIRING* ctx, const Fp12* GTsig)
    {
        return (int)PAIRING_FinalVerify(ctx, GTsig);
    }

    BLST_HD int blst_fp12_finalverify(const Fp12& GT1, const Fp12& GT2)
    {
        Fp12 GT;

        vec_copy(&GT, &GT1, sizeof(GT));
        Fp12::conjugate(GT);
        Fp12::mul(GT, GT, GT2);
        Fp12::final_exp(GT, GT);

        /* return GT==1 */
        return (int)(vec_is_equal(&GT[0][0], &get_BLST_Rx().p2, sizeof(GT[0][0])) &
            vec_is_zero(&GT[0][1], sizeof(GT) - sizeof(GT[0][0])));
    }

    BLST_HD void blst_pairing_raw_aggregate(PAIRING* ctx, const G2AffinePoint* q,
        const G1AffinePoint* p)
    {
        unsigned int n;

        if (vec_is_zero(q, sizeof(*q)) & vec_is_zero(p, sizeof(*p)))
            return;

        n = ctx->nelems;
        vec_copy(ctx->Q + n, q, sizeof(*q));
        vec_copy(ctx->P + n, p, sizeof(*p));
        if (++n == N_MAX) {
            if (ctx->ctrl & AGGR_GT_SET) {
                Fp12 GT;
                miller_loop_n(GT, ctx->Q, ctx->P, n);
                Fp12::mul(ctx->GT, ctx->GT, GT);
            }
            else {
                miller_loop_n(ctx->GT, ctx->Q, ctx->P, n);
                ctx->ctrl |= AGGR_GT_SET;
            }
            n = 0;
        }
        ctx->nelems = n;
    }

    BLST_HD Fp12* blst_pairing_as_fp12(PAIRING* ctx)
    {
        PAIRING_Commit(ctx);
        return (Fp12*)&ctx->GT;
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
        G1Point P;
        BLST_ERROR ret;

        ret = G1AffinePoint::deserialize(*(G1AffinePoint*)&P, zwire);

        if (ret != BLST_SUCCESS)
            return ret;

        if (vec_is_zero(&P, sizeof(G1AffinePoint))) {
            if (in == NULL)
                vec_zero(&out, sizeof(out));
            return BLST_SUCCESS;
        }

        vec_copy(P.Z, get_BLST_Rx().p, sizeof(P.Z));

        if (!G1Point::in_subgroup(P))
            return BLST_POINT_NOT_IN_GROUP;

        if (in == NULL)
            vec_copy(&out, &P, sizeof(P));
        else
            G1Point::dadd_affine(out, *in, *(G1AffinePoint*)&P);

        return BLST_SUCCESS;
    }

    BLST_HD BLST_ERROR blst_aggregate_in_g2(G2Point& out, const G2Point* in,
        const unsigned char* zwire)
    {
        G2Point P[1];
        BLST_ERROR ret;

        ret = G2AffinePoint::deserialize(*(G2AffinePoint*)P, zwire);

        if (ret != BLST_SUCCESS)
            return ret;

        if (vec_is_zero(P, sizeof(G2AffinePoint))) {
            if (in == NULL)
                vec_zero(&out, sizeof(out));
            return BLST_SUCCESS;
        }

        vec_copy(&P->Z, get_BLST_Rx().p, sizeof(P->Z));

        if (!G2Point::in_subgroup(*P))
            return BLST_POINT_NOT_IN_GROUP;

        if (in == NULL) {
            vec_copy(&out, P, sizeof(P));
        }
        else {
            G2Point::dadd_affine(out, *in, *(G2AffinePoint*)P);
        }
        return BLST_SUCCESS;
    }

    BLST_HD void blst_aggregated_in_g1(Fp12& ret, const G1AffinePoint& sig)
    {
        miller_loop_n(ret, (const G2AffinePoint*)&get_BLST_G2(), &sig, 1);
    }

    BLST_HD void blst_aggregated_in_g2(Fp12& ret, const G2AffinePoint& sig)
    {
        miller_loop_n(ret, &sig, (const G1AffinePoint*)&get_BLST_G1(), 1);
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