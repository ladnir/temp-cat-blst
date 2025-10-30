
#pragma once
#include <cstdint>
#include "blst/common/defines.h"

#ifdef BLST_CUDA

// ptx_add_mod_384_u64
//#include "cuda/add_mod_384_cuda.h"
template<typename Field>
BLST_DEVICE 
void ptx_add_mod_u64(
    uint64_t* BLST_restrict r,
    const uint64_t* BLST_restrict a,
    const uint64_t* BLST_restrict b);

template<typename Field>
BLST_HD 
void ptx_mul_mont_u64(
    uint64_t* BLST_restrict r,
    const uint64_t* BLST_restrict a,
    const uint64_t* BLST_restrict b
);

#endif // BLST_CUDA