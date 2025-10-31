#include "Fp.h"

#ifdef BLST_CUDA
#include "blst/backend/cuda/add_mod_cuda.h"
#include "blst/backend/cuda/mul_mont_cuda.h"

template 
BLST_DEVICE void ptx_add_mod_u64<blst::Fp>(
	uint64_t* BLST_restrict r,
	const uint64_t* BLST_restrict a,
	const uint64_t* BLST_restrict b);

template
BLST_DEVICE void ptx_sub_mod_u64<blst::Fp>(
	uint64_t* BLST_restrict r,
	const uint64_t* BLST_restrict a,
	const uint64_t* BLST_restrict b);

template
BLST_HD void ptx_mul_mont_u64<blst::Fp>(
	uint64_t* BLST_restrict r,
	const uint64_t* BLST_restrict a,
	const uint64_t* BLST_restrict b
);
template
BLST_HD void ptx_sqr_mont_u64<blst::Fp>(
	uint64_t* BLST_restrict r,
	const uint64_t* BLST_restrict a);

#endif
