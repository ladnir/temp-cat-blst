/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef BLST_BTX_H
#define BLST_BTX_H

#include "blst.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifndef BLST_FILL_RANDOM_FN_DEFINED
#define BLST_FILL_RANDOM_FN_DEFINED
typedef bool (*blst_fill_random_fn)(void* user_data, blst_byte* out, size_t out_len);
#endif

typedef struct blst_btx_setup_result blst_btx_setup_result;
typedef struct blst_btx_public_parameters blst_btx_public_parameters;
typedef struct blst_btx_prepared_public_parameters blst_btx_prepared_public_parameters;
typedef struct blst_btx_server_secret_share blst_btx_server_secret_share;
typedef struct blst_btx_batch_precomputation blst_btx_batch_precomputation;

typedef enum blst_btx_cross_term_algorithm {
    BLST_BTX_CROSS_TERM_NAIVE = 0,
    BLST_BTX_CROSS_TERM_FFT = 1,
    BLST_BTX_CROSS_TERM_MULTIPAIR = 2,
    BLST_BTX_CROSS_TERM_AUTO = 3,
} blst_btx_cross_term_algorithm;

typedef struct blst_btx_params {
    size_t max_batch_size;
    size_t num_servers;
    size_t threshold;
    size_t key_size;
    const char* kdf_salt;
    const char* kdf_info;
} blst_btx_params;

typedef struct blst_btx_partial_decryption_share {
    uint32_t index;
    blst_p1 value;
} blst_btx_partial_decryption_share;

blst_btx_setup_result* blst_btx_setup(
    const blst_btx_params* params,
    blst_fill_random_fn fill_random,
    void* user_data);

const blst_btx_public_parameters* blst_btx_setup_result_public_parameters(
    const blst_btx_setup_result* setup_result);

size_t blst_btx_setup_result_secret_share_count(
    const blst_btx_setup_result* setup_result);

const blst_btx_server_secret_share* blst_btx_setup_result_secret_share(
    const blst_btx_setup_result* setup_result,
    size_t index);

void blst_btx_setup_result_free(blst_btx_setup_result* setup_result);

blst_btx_prepared_public_parameters* blst_btx_prepare_public_parameters(
    const blst_btx_public_parameters* public_parameters);

void blst_btx_prepared_public_parameters_free(
    blst_btx_prepared_public_parameters* prepared_public_parameters);

bool blst_btx_encrypt(
    const blst_btx_params* params,
    const blst_btx_public_parameters* public_parameters,
    blst_fill_random_fn fill_random,
    void* user_data,
    blst_p1* out_c1,
    blst_byte* out_key,
    size_t out_key_len);

bool blst_btx_create_partial_decryption_share(
    const blst_btx_server_secret_share* secret_share,
    const blst_p1* ordered_c1,
    size_t n,
    blst_btx_partial_decryption_share* out_share);

bool blst_btx_combine_partial_decryption_shares(
    const blst_btx_partial_decryption_share* shares,
    size_t n,
    const blst_btx_public_parameters* public_parameters,
    blst_p1* out_aggregate_predecryption_key);

blst_btx_batch_precomputation* blst_btx_precompute_batch_decryption(
    const blst_btx_prepared_public_parameters* prepared_public_parameters,
    const blst_p1* ciphertexts,
    size_t n,
    blst_btx_cross_term_algorithm algorithm);

size_t blst_btx_batch_precomputation_batch_size(
    const blst_btx_batch_precomputation* batch_precomputation);

void blst_btx_batch_precomputation_free(
    blst_btx_batch_precomputation* batch_precomputation);

bool blst_btx_finalize_batch_decryption(
    const blst_btx_params* params,
    const blst_btx_prepared_public_parameters* prepared_public_parameters,
    const blst_btx_batch_precomputation* batch_precomputation,
    const blst_p1* aggregate_predecryption_key,
    blst_byte* out_keys,
    size_t out_keys_len);

bool blst_btx_recover_batch_keys(
    const blst_btx_params* params,
    const blst_btx_prepared_public_parameters* prepared_public_parameters,
    const blst_p1* ciphertexts,
    size_t n,
    const blst_p1* aggregate_predecryption_key,
    blst_btx_cross_term_algorithm algorithm,
    blst_byte* out_keys,
    size_t out_keys_len);

#ifdef __cplusplus
}
#endif

#endif
