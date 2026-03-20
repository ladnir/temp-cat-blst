#include "blst/scalar/fields/Fp.h"

#ifdef BLST_CUDA
#include "blst/backend/cuda/add_mod_cuda.h"
#include "blst/backend/cuda/mul_mont_cuda.h"

namespace blst::scalar
{
	template
		BLST_DEVICE void ptx_add_mod_u64<blst::scalar::Fp>(
			uint64_t* r,
			const uint64_t* a,
			const uint64_t* b);

	template
		BLST_DEVICE void ptx_sub_mod_u64<blst::scalar::Fp>(
			uint64_t* r,
			const uint64_t* a,
			const uint64_t* b);

	template
		BLST_HD void ptx_mul_mont_u64<blst::scalar::Fp>(
			uint64_t* r,
			const uint64_t* a,
			const uint64_t* b
		);
	template
		BLST_HD void ptx_sqr_mont_u64<blst::scalar::Fp>(
			uint64_t* r,
			const uint64_t* a);
}

#endif
