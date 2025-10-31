
#pragma once
#include <cstdint>
#include "blst/common/defines.h"

#ifdef BLST_CUDA

// Addition modulo p of two field elements a and b.
template<typename Field>
BLST_DEVICE 
void ptx_add_mod_u64(
    uint64_t* BLST_restrict r,
    const uint64_t* BLST_restrict a,
    const uint64_t* BLST_restrict b);


template<typename Field>
BLST_DEVICE 
void ptx_sub_mod_u64(
    uint64_t* BLST_restrict r,
    const uint64_t* BLST_restrict a,
    const uint64_t* BLST_restrict b);

// Montgomery multiplication of two field elements a and b.
template<typename Field>
BLST_HD 
void ptx_mul_mont_u64(
    uint64_t* BLST_restrict r,
    const uint64_t* BLST_restrict a,
    const uint64_t* BLST_restrict b
);

// Montgomery squaring of a field element a.
template<typename Field>
BLST_HD 
void ptx_sqr_mont_u64(
    uint64_t* BLST_restrict r,
    const uint64_t* BLST_restrict a);




#endif // BLST_CUDA