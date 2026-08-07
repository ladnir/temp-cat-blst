#include "blst/avx/common/preprocess.h"
#include "blst/avx/cpu/vector/vector_impl.hpp"
#include "blst/avx/rns/elementwise_reduce.h"

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
