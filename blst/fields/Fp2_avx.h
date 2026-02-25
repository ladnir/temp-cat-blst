#ifndef BLST_FIELDS_FP2_AVX_H
#define BLST_FIELDS_FP2_AVX_H

#include "blst/fields/Fp_avx.h"

#if defined(BLST_AVX512) || defined(BLST_AVX512_PORTABLE)

namespace blst
{
    /**
     * AVX512 Fp2 helpers ported from avxbls.
     *
     * All inputs/outputs are Montgomery-form field elements packed into AVX lanes.
     * Layout-specific function names are kept so kernels can be compared line-for-line
     * with avxbls.
     */

    template <class PackT, std::size_t N>
    struct AvxPackTuple
    {
        PackT packs[N];

        BLST_FORCEINLINE PackT& operator[](std::size_t i)
        {
            BLST_ASSERT(i < N);
            return packs[i];
        }

        BLST_FORCEINLINE const PackT& operator[](std::size_t i) const
        {
            BLST_ASSERT(i < N);
            return packs[i];
        }

        BLST_FORCEINLINE PackT* data()
        {
            return packs;
        }

        BLST_FORCEINLINE const PackT* data() const
        {
            return packs;
        }
    };

    /// Two Fp vectors (c0,c1), each in `fp_8x1w`.
    using Fp2Pack8x1x1 = AvxPackTuple<fp_8x1w, 2>;
    using fp2_8x1x1w = Fp2Pack8x1x1;
    /// Widened Fp2 pair for `fp2_8x1x1w`.
    using UFp2Pack8x1x1 = AvxPackTuple<fpx2_8x1w, 2>;
    using fp2x2_8x1x1w = UFp2Pack8x1x1;

    /// Four packed Fp2 values in single-length schedule (`8` vector words).
    using Fp2Pack4x2x1 = fp_8x1w;
    using fp2_4x2x1w = Fp2Pack4x2x1;
    /// Widened companion layout for `fp2_4x2x1w`.
    using UFp2Pack4x2x1 = fpx2_8x1w;
    using fp2x2_4x2x1w = UFp2Pack4x2x1;

    /// Single-length schedule used by `mul_fp2_2x4x1w`.
    using Fp2Pack2x4x1 = fp_8x1w;
    using fp2_2x4x1w = Fp2Pack2x4x1;
    /// Widened companion layout for `fp2_2x4x1w`.
    using UFp2Pack2x4x1 = fpx2_8x1w;
    using fp2x2_2x4x1w = UFp2Pack2x4x1;

    /// Two packed Fp2 values in half-length schedule (`4` vector words).
    using Fp2Pack2x2x2 = fp_4x2w;
    using fp2_2x2x2w = Fp2Pack2x2x2;
    /// Widened companion layout for `fp2_2x2x2w`.
    using UFp2Pack2x2x2 = fpx2_4x2w;
    using fp2x2_2x2x2w = UFp2Pack2x2x2;

    /// Single-length schedule used by `mul_fp2_1x4x2w`.
    using Fp2Pack1x4x2 = fp_4x2w;
    using fp2_1x4x2w = Fp2Pack1x4x2;
    /// Widened companion layout for `fp2_1x4x2w`.
    using UFp2Pack1x4x2 = fpx2_4x2w;
    using fp2x2_1x4x2w = UFp2Pack1x4x2;

    static_assert(std::is_trivially_default_constructible_v<Fp2Pack8x1x1>);
    static_assert(std::is_trivially_copyable_v<Fp2Pack8x1x1>);
    static_assert(std::is_trivially_default_constructible_v<UFp2Pack8x1x1>);
    static_assert(std::is_trivially_copyable_v<UFp2Pack8x1x1>);

    /** @name Fp2 single-length helpers */
    ///@{
    /// Component-wise add for `fp2_8x1x1w`.
    void add_fp2_8x1x1w(fp2_8x1x1w& r, const fp2_8x1x1w& a, const fp2_8x1x1w& b);
    /// Component-wise add for `fp2_4x2x1w`.
    void add_fp2_4x2x1w(fp2_4x2x1w& r, const fp2_4x2x1w& a, const fp2_4x2x1w& b);
    /// Component-wise subtract for `fp2_4x2x1w`.
    void sub_fp2_4x2x1w(fp2_4x2x1w& r, const fp2_4x2x1w& a, const fp2_4x2x1w& b);

    /// Mixed add/sub schedule used by fp12 kernels (`aass`).
    void aass_fp2_4x2x1w(fp2_4x2x1w& r, const fp2_4x2x1w& a, const fp2_4x2x1w& b);
    /// Mixed add/sub schedule used by fp12 kernels (`assa`).
    void assa_fp2_4x2x1w(fp2_4x2x1w& r, const fp2_4x2x1w& a, const fp2_4x2x1w& b);
    /// Mixed add/sub schedule used by fp12 kernels (`sass`).
    void sass_fp2_4x2x1w(fp2_4x2x1w& r, const fp2_4x2x1w& a, const fp2_4x2x1w& b);
    /// Mixed add/sub schedule used by fp12 kernels (`ssaa`).
    void ssaa_fp2_4x2x1w(fp2_4x2x1w& r, const fp2_4x2x1w& a, const fp2_4x2x1w& b);

    /// Multiply by `(u+1)` in `fp2_4x2x1w` layout.
    void mul_by_u_plus_1_fp2_4x2x1w(fp2_4x2x1w& r, const fp2_4x2x1w& a);

    /// Schoolbook Fp2 multiply in `fp2_4x2x1w` schedule.
    void mul_fp2_4x2x1w(fp2_4x2x1w& r, const fp2_4x2x1w& a, const fp2_4x2x1w& b);
    /// Specialized Fp2 square in `fp2_4x2x1w` schedule.
    void sqr_fp2_4x2x1w(fp2_4x2x1w& r, const fp2_4x2x1w& a);

    /// Schoolbook Fp2 multiply in `fp2_2x4x1w` schedule.
    void mul_fp2_2x4x1w(fp2_2x4x1w& r, const fp2_2x4x1w& a, const fp2_2x4x1w& b);

    /// Component-wise add for `fp2_2x2x2w`.
    void add_fp2_2x2x2w(fp2_2x2x2w& r, const fp2_2x2x2w& a, const fp2_2x2x2w& b);
    /// Component-wise subtract for `fp2_2x2x2w`.
    void sub_fp2_2x2x2w(fp2_2x2x2w& r, const fp2_2x2x2w& a, const fp2_2x2x2w& b);
    /// Mixed add/sub schedule for `fp2_2x2x2w` (`as`).
    void as_fp2_2x2x2w(fp2_2x2x2w& r, const fp2_2x2x2w& a, const fp2_2x2x2w& b);
    /// Mixed add/sub schedule for `fp2_2x2x2w` (`sa`).
    void sa_fp2_2x2x2w(fp2_2x2x2w& r, const fp2_2x2x2w& a, const fp2_2x2x2w& b);
    /// Multiply by `(u+1)` in `fp2_2x2x2w` layout.
    void mul_by_u_plus_1_fp2_2x2x2w(fp2_2x2x2w& r, const fp2_2x2x2w& a);
    /// Schoolbook Fp2 multiply in `fp2_2x2x2w` schedule.
    void mul_fp2_2x2x2w(fp2_2x2x2w& r, const fp2_2x2x2w& a, const fp2_2x2x2w& b);
    /// Specialized Fp2 square in `fp2_2x2x2w` schedule.
    void sqr_fp2_2x2x2w(fp2_2x2x2w& r, const fp2_2x2x2w& a);

    /// Schoolbook Fp2 multiply in `fp2_1x4x2w` schedule.
    void mul_fp2_1x4x2w(fp2_1x4x2w& r, const fp2_1x4x2w& a, const fp2_1x4x2w& b);
    ///@}

    /** @name Fp2 widened helpers and kernels */
    ///@{
    /// Component-wise widened add for `fp2x2_8x1x1w`.
    void add_fp2x2_8x1x1w(fp2x2_8x1x1w& r, const fp2x2_8x1x1w& a, const fp2x2_8x1x1w& b);
    /// Component-wise widened subtract for `fp2x2_8x1x1w`.
    void sub_fp2x2_8x1x1w(fp2x2_8x1x1w& r, const fp2x2_8x1x1w& a, const fp2x2_8x1x1w& b);
    /// Mixed widened schedule (`asx4`) for `fp2x2_8x1x1w`.
    void asx4_fp2x2_8x1x1w(fp2x2_8x1x1w& r, const fp2x2_8x1x1w& a, const fp2x2_8x1x1w& b);
    /// Multiply by `(u+1)` in widened `fp2x2_8x1x1w`.
    void mul_by_u_plus_1_fp2x2_8x1x1w(fp2x2_8x1x1w& r, const fp2x2_8x1x1w& a);

    /// Karatsuba widened Fp2 multiply in `fp2x2_8x1x1w`.
    void mul_fp2x2_8x1x1w(fp2x2_8x1x1w& r, const fp2_8x1x1w& a, const fp2_8x1x1w& b);

    /// Component-wise widened add for `fp2x2_4x2x1w`.
    void add_fp2x2_4x2x1w(fp2x2_4x2x1w& r, const fp2x2_4x2x1w& a, const fp2x2_4x2x1w& b);
    /// Component-wise widened subtract for `fp2x2_4x2x1w`.
    void sub_fp2x2_4x2x1w(fp2x2_4x2x1w& r, const fp2x2_4x2x1w& a, const fp2x2_4x2x1w& b);
    /// Mixed widened schedule (`sax2`) for `fp2x2_4x2x1w`.
    void sax2_fp2x2_4x2x1w(fp2x2_4x2x1w& r, const fp2x2_4x2x1w& a, const fp2x2_4x2x1w& b);
    /// Multiply by `(u+1)` in widened `fp2x2_4x2x1w`.
    void mul_by_u_plus_1_fp2x2_4x2x1w(fp2x2_4x2x1w& r, const fp2x2_4x2x1w& a);
    /// Schoolbook widened multiply in `fp2x2_4x2x1w`.
    void mul_fp2x2_4x2x1w(fp2x2_4x2x1w& r, const fp2_4x2x1w& a, const fp2_4x2x1w& b);
    /// Specialized widened square in `fp2x2_4x2x1w`.
    void sqr_fp2x2_4x2x1w(fp2x2_4x2x1w& r, const fp2_4x2x1w& a);
    /// Reduce widened `fp2x2_4x2x1w` to single-length `fp2_4x2x1w`.
    void redc_fp2x2_4x2x1w(fp2_4x2x1w& r, const fp2x2_4x2x1w& a);

    /// Schoolbook widened multiply in `fp2x2_2x4x1w`.
    void mul_fp2x2_2x4x1w(fp2x2_2x4x1w& r, const fp2_2x4x1w& a, const fp2_2x4x1w& b);

    /// Reduce widened `fp2x2_2x2x2w` to single-length `fp2_2x2x2w`.
    void redc_fp2x2_2x2x2w(fp2_2x2x2w& r, const fp2x2_2x2x2w& a);
    /// Component-wise widened add for `fp2x2_2x2x2w`.
    void add_fp2x2_2x2x2w(fp2x2_2x2x2w& r, const fp2x2_2x2x2w& a, const fp2x2_2x2x2w& b);
    /// Component-wise widened subtract for `fp2x2_2x2x2w`.
    void sub_fp2x2_2x2x2w(fp2x2_2x2x2w& r, const fp2x2_2x2x2w& a, const fp2x2_2x2x2w& b);
    /// Multiply by `(u+1)` in widened `fp2x2_2x2x2w`.
    void mul_by_u_plus_1_fp2x2_2x2x2w(fp2x2_2x2x2w& r, const fp2x2_2x2x2w& a);
    /// Schoolbook widened multiply in `fp2x2_2x2x2w`.
    void mul_fp2x2_2x2x2w(fp2x2_2x2x2w& r, const fp2_2x2x2w& a, const fp2_2x2x2w& b);
    /// Specialized widened square in `fp2x2_2x2x2w`.
    void sqr_fp2x2_2x2x2w(fp2x2_2x2x2w& r, const fp2_2x2x2w& a);

    /// Schoolbook widened multiply in `fp2x2_1x4x2w`.
    void mul_fp2x2_1x4x2w(fp2x2_1x4x2w& r, const fp2_1x4x2w& a, const fp2_1x4x2w& b);
    ///@}
}

#endif

#endif // BLST_FIELDS_FP2_AVX_H
