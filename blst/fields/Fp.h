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
	using scalar::Fp;
	using scalar::radix384;

	using scalar::BLST_NEG_P_HOST;
	using scalar::BLST_P_HOST;
	using scalar::BLST_RR_HOST;
	using scalar::BLST_RRRR_HOST;
	using scalar::BLST_Rx_HOST;

	using scalar::ct_inverse_mod_v384;
	using scalar::get_BLST_negP;
	using scalar::get_BLST_P;
	using scalar::get_BLST_RR;
	using scalar::get_BLST_RRRR;
	using scalar::get_BLST_Rx;
}

#endif /* BLST_FP_H */
