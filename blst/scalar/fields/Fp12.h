#ifndef BLST_SCALAR_FP12_H
#define BLST_SCALAR_FP12_H


#include "blst/common/defines.h"
#include "blst/backend/backend.h"
#include "blst/scalar/fields/Fp2.h"



namespace blst::scalar
{

	struct Fp4 {
		Fp2 data[2];
		BLST_HD Fp2& operator[](size_t idx) {
			return data[idx];
		}
		BLST_HD const Fp2& operator[](size_t idx)const {
			return data[idx];
		}
	};


	/// @brief Precomputed schedule for constant-time GLV-4 exponentiation.
	/// Contains a compact bytecode sequence of operations (squares/multiplies).

	/// @brief Precomputed schedule for constant-time GLV-2 JSF exponentiation.
	/// Contains a compact bytecode sequence of operations (squares/multiplies).
	/// The ecnoding has two signed digits per index. Then we pack two such
	/// pairs into a byte for storage efficiency.



}

#endif /* BLST_SCALAR_FP12_H */



