#include "blst/avx/common/preprocess.h"
#include "blst/avx/cpu/vector/vector_impl.hpp"
#include "blst/avx/rns/elementwise_reduce.h"
#include "blst/avx/rns/matrixnok_change_base.h"
#include "blst/avx/rns/matrixnok_fp.h"
#include "blst/avx/fields/matrixnok_fp2.h"
#include "blst/avx/fields/matrixnok_fp12.h"
#include "blst/avx/curve/matrixnok_g1.h"
#include "blst/avx/curve/matrixnok_g2.h"
#include "blst/avx/gt/glv4.h"

#include <cstdint>

#if !defined(BLST_AVX512)
#error "MatrixNoK code-generation probes require native AVX-512"
#endif

using blst::avx::AVXVector;
using blst::avx::BoundedElement;
using blst::avx::Bounds;
using blst::avx::MontgomeryReduce2;
using blst::avx::ShortVector;

extern "C" NOINLINE uint64_t audit_avx_extract_lane0(const AVXVector<8>* value)
{
    return value->extract0();
}

extern "C" NOINLINE void audit_avx_broadcast_reduction_quotient(
    uint64_t out[8], const AVXVector<8>* high, const AVXVector<8>* low)
{
    const ShortVector low_part = ShortVector::from_avx_lo(*low).srli(52);
    const ShortVector quotient = ShortVector::from_avx_lo(*high).add(low_part);
    quotient.broadcast_to_avx<8>().store(out);
}

extern "C" NOINLINE void audit_avx_reduce_bound40(
    uint64_t out[8], const MontgomeryReduce2<8, 50, 6>* reducer, const uint64_t in[8])
{
    using Input = BoundedElement<Bounds<0, 40>, 8, 50>;
    AVXVector<8> value;
    value.load(in);
    const auto reduced = reducer->min_sub_reduce<40, 1>(Input(value));
    reduced.data.store(out);
}

extern "C" NOINLINE void audit_avx_matrixnok_r1_change_base(
    uint64_t out_hi[8], uint64_t out_lo[8], const uint64_t residues_in[8],
    const uint64_t initial_hi[8], const uint64_t initial_lo[8])
{
    std::array<AVXVector<8>, 1> residues;
    std::array<AVXVector<8>, 1> high;
    std::array<AVXVector<8>, 1> low;
    residues[0].load(residues_in);
    high[0].load(initial_hi);
    low[0].load(initial_lo);
    blst::avx::matrixnok::r1_change_base_raw_batch(residues, high, low);
    high[0].store(out_hi);
    low[0].store(out_lo);
}

extern "C" NOINLINE void audit_avx_matrixnok_fp_multiply(
    blst::avx::matrixnok::FpRingCandidate::Element* out,
    const blst::avx::matrixnok::FpRingCandidate::Element* left,
    const blst::avx::matrixnok::FpRingCandidate::Element* right)
{
    *out = blst::avx::matrixnok::FpRingCandidate::modmul(*left, *right);
}

extern "C" NOINLINE void audit_avx_matrixnok_fp2_multiply(
    blst::avx::matrixnok::Fp2RingCandidate::Element* out,
    const blst::avx::matrixnok::Fp2RingCandidate::Element* left,
    const blst::avx::matrixnok::Fp2RingCandidate::Element* right)
{
    const blst::avx::matrixnok::Fp2RingCandidate fp2ring;
    fp2ring.multiply(*out, *left, *right);
}

extern "C" NOINLINE void audit_avx_matrixnok_fp2_square(
    blst::avx::matrixnok::Fp2RingCandidate::Element* out,
    const blst::avx::matrixnok::Fp2RingCandidate::Element* value)
{
    const blst::avx::matrixnok::Fp2RingCandidate fp2ring;
    fp2ring.square(*out, *value);
}

extern "C" NOINLINE void audit_avx_matrixnok_fp12_multiply(
    blst::avx::matrixnok::Fp12RingCandidate::StandardElement* out,
    const blst::avx::matrixnok::Fp12RingCandidate::StandardElement* left,
    const blst::avx::matrixnok::Fp12RingCandidate::StandardElement* right)
{
    const blst::avx::matrixnok::FpRingCandidate ring;
    const blst::avx::matrixnok::Fp12RingCandidate fp12ring(ring);
    *out = fp12ring.fp12_flat_mul(*left, *right, ring);
}

extern "C" NOINLINE void audit_avx_matrixnok_fp12_cyclotomic_square(
    blst::avx::matrixnok::Fp12RingCandidate::StandardElement* out,
    const blst::avx::matrixnok::Fp12RingCandidate::StandardElement* value)
{
    const blst::avx::matrixnok::FpRingCandidate ring;
    const blst::avx::matrixnok::Fp12RingCandidate fp12ring(ring);
    *out = fp12ring.narrow_to_standard(
        fp12ring.fp12_flat_cyclotomic_sqr(*value, ring), ring);
}

extern "C" NOINLINE void audit_avx_matrixnok_gt_exp_ct(
    blst::avx::matrixnok::Fp12RingCandidate::StandardElement* out,
    const blst::avx::matrixnok::Fp12RingCandidate::StandardElement* value,
    const uint64_t scalar[4])
{
    const blst::avx::matrixnok::FpRingCandidate ring;
    const blst::avx::matrixnok::Fp12RingCandidate fp12ring(ring);
    *out = blst::avx::avx_cyclotomic_exp_glv4_ct(
        *value, scalar, fp12ring, ring);
}

extern "C" NOINLINE void audit_avx_matrixnok_g1_add(
    blst::avx::matrixnok::G1ProjectiveCandidate* out,
    const blst::avx::matrixnok::G1ProjectiveCandidate* left,
    const blst::avx::matrixnok::G1ProjectiveCandidate* right)
{
    const blst::avx::matrixnok::FpRingCandidate ring;
    const blst::avx::matrixnok::G1CurveCandidate curve;
    *out = curve.add_point(*left, *right, ring);
}

extern "C" NOINLINE void audit_avx_matrixnok_g1_double(
    blst::avx::matrixnok::G1ProjectiveCandidate* out,
    const blst::avx::matrixnok::G1ProjectiveCandidate* value)
{
    const blst::avx::matrixnok::FpRingCandidate ring;
    const blst::avx::matrixnok::G1CurveCandidate curve;
    *out = curve.double_point(*value, ring);
}

extern "C" NOINLINE void audit_avx_matrixnok_g2_add(
    blst::avx::matrixnok::G2ProjectiveCandidate* out,
    const blst::avx::matrixnok::G2ProjectiveCandidate* left,
    const blst::avx::matrixnok::G2ProjectiveCandidate* right)
{
    const blst::avx::matrixnok::FpRingCandidate ring;
    const blst::avx::matrixnok::G2Fp2RingCandidate fp2ring(ring);
    const blst::avx::matrixnok::G2CurveCandidate curve;
    *out = curve.add_point(*left, *right, fp2ring);
}

extern "C" NOINLINE void audit_avx_matrixnok_g2_double(
    blst::avx::matrixnok::G2ProjectiveCandidate* out,
    const blst::avx::matrixnok::G2ProjectiveCandidate* value)
{
    const blst::avx::matrixnok::FpRingCandidate ring;
    const blst::avx::matrixnok::G2Fp2RingCandidate fp2ring(ring);
    const blst::avx::matrixnok::G2CurveCandidate curve;
    *out = curve.double_point(*value, fp2ring);
}

extern "C" NOINLINE void audit_avx_matrixnok_g2_gls_ct(
    blst::avx::matrixnok::G2ProjectiveCandidate* out,
    const blst::avx::matrixnok::G2ProjectiveCandidate* value,
    const std::array<uint8_t, 32>* scalar,
    const blst::avx::matrixnok::G2GLSCandidate* gls)
{
    const blst::avx::matrixnok::FpRingCandidate ring;
    const blst::avx::matrixnok::G2Fp2RingCandidate fp2ring(ring);
    const blst::avx::matrixnok::G2CurveCandidate curve;
    *out = gls->multiply(*value, *scalar, fp2ring, curve, ring);
}
