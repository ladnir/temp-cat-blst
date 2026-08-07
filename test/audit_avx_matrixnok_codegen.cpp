#include "blst/avx/common/preprocess.h"
#include "blst/avx/cpu/vector/vector_impl.hpp"

#include <cstdint>

#if !defined(BLST_AVX512)
#error "MatrixNoK code-generation probes require native AVX-512"
#endif

using blst::avx::AVXVector;
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
