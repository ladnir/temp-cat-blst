/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef BLST_FP_H
#define BLST_FP_H

#include "blst/scalar/fields/Fp.h"

namespace blst
{
	/**
	 * @brief Main native C++ BLS12-381 base-field type.
	 *
	 * This is an exact alias of `scalar::Fp`; its representation, preconditions,
	 * aliasing, and timing contracts are documented on that implementation type.
	 */
	using scalar::Fp;

	using scalar::BLST_NEG_P_HOST;
	using scalar::BLST_P_HOST;
	using scalar::BLST_RR_HOST;
	using scalar::BLST_RRRR_HOST;

	using scalar::ct_inverse_mod_v384;
	using scalar::get_BLST_negP;
	using scalar::get_BLST_P;
	using scalar::get_BLST_RR;
	using scalar::get_BLST_RRRR;
}

#endif /* BLST_FP_H */
