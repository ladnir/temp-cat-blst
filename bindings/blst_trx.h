/*
 * Copyright Supranational LLC
 * Licensed under the Apache License, Version 2.0, see LICENSE for details.
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef __BLST_TRX_H__
#define __BLST_TRX_H__

#include "blst.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifndef BLST_FILL_RANDOM_FN_DEFINED
#define BLST_FILL_RANDOM_FN_DEFINED
typedef bool (*blst_fill_random_fn)(void* user_data, byte* out, size_t out_len);
#endif

typedef struct blst_trx_crs blst_trx_crs;
typedef struct blst_trx_public_key_material blst_trx_public_key_material;
typedef struct blst_trx_encryption_key blst_trx_encryption_key;
typedef struct blst_trx_encrypted_ciphertext blst_trx_encrypted_ciphertext;
typedef struct blst_trx_verified_ciphertext blst_trx_verified_ciphertext;
typedef struct blst_trx_batch_digest blst_trx_batch_digest;

typedef struct blst_trx_params {
    size_t max_batch_size;
    size_t num_contexts;
    size_t session_key_size;
    uint8_t serialization_version;
    const char* h1_dst;
    const char* hf_dst;
    const char* kdf_salt;
    const char* kdf_info;
    const char* ots_sign_dst;
} blst_trx_params;

typedef struct blst_trx_validator_secret_share {
    uint32_t index;
    blst_scalar sk_share;
} blst_trx_validator_secret_share;

typedef struct blst_trx_validator_public_share {
    uint32_t index;
    blst_p2 pk_share;
    blst_p2 pk_tau_share;
} blst_trx_validator_public_share;

typedef struct blst_trx_partial_decryption_share {
    uint32_t index;
    blst_p1 value;
} blst_trx_partial_decryption_share;

blst_trx_crs* blst_trx_generate_crs(
    const blst_trx_params* params,
    blst_fill_random_fn fill_random,
    void* user_data);

void blst_trx_crs_free(blst_trx_crs* crs);

bool blst_trx_derive_validator_public_share(
    const blst_trx_validator_secret_share* secret_share,
    const blst_trx_crs* crs,
    blst_trx_validator_public_share* out_public_share);

blst_trx_public_key_material* blst_trx_derive_public_key_material(
    const blst_trx_validator_public_share* public_shares,
    size_t public_share_count,
    const uint32_t* subset_indices,
    size_t subset_count);

void blst_trx_public_key_material_free(
    blst_trx_public_key_material* public_key_material);

blst_trx_encryption_key* blst_trx_to_encryption_key(
    const blst_trx_public_key_material* public_key_material);

void blst_trx_encryption_key_free(
    blst_trx_encryption_key* encryption_key);

blst_trx_encrypted_ciphertext* blst_trx_encrypt(
    const blst_trx_params* params,
    const blst_trx_encryption_key* encryption_key,
    const byte* associated_data,
    size_t associated_data_len,
    blst_fill_random_fn fill_random,
    void* user_data);

size_t blst_trx_encrypted_ciphertext_session_key_size(
    const blst_trx_encrypted_ciphertext* encrypted_ciphertext);

bool blst_trx_encrypted_ciphertext_copy_session_key(
    const blst_trx_encrypted_ciphertext* encrypted_ciphertext,
    byte* out_session_key,
    size_t out_session_key_len);

void blst_trx_encrypted_ciphertext_free(
    blst_trx_encrypted_ciphertext* encrypted_ciphertext);

blst_trx_verified_ciphertext* blst_trx_verify_ciphertext(
    const blst_trx_params* params,
    const blst_trx_encrypted_ciphertext* encrypted_ciphertext);

void blst_trx_verified_ciphertext_free(
    blst_trx_verified_ciphertext* verified_ciphertext);

blst_trx_batch_digest* blst_trx_build_batch_digest(
    const blst_trx_params* params,
    const blst_trx_verified_ciphertext* const* verified_ciphertexts,
    size_t verified_ciphertext_count,
    size_t context_id,
    const blst_trx_crs* crs);

size_t blst_trx_batch_digest_valid_ciphertext_count(
    const blst_trx_batch_digest* batch_digest);

bool blst_trx_compute_eval_proofs(
    const blst_trx_batch_digest* batch_digest,
    const blst_trx_crs* crs,
    blst_p1* out_proofs,
    size_t out_proof_count);

bool blst_trx_verify_evaluation_proof(
    const blst_trx_batch_digest* batch_digest,
    size_t proof_index,
    const blst_p1* proof,
    const blst_trx_crs* crs);

bool blst_trx_create_partial_decryption(
    const blst_trx_params* params,
    const blst_trx_validator_secret_share* secret_share,
    const blst_trx_batch_digest* batch_digest,
    const blst_trx_public_key_material* public_key_material,
    blst_trx_partial_decryption_share* out_share);

bool blst_trx_verify_partial_decryption(
    const blst_trx_params* params,
    const blst_trx_partial_decryption_share* share,
    const blst_trx_validator_public_share* public_share,
    const blst_trx_batch_digest* batch_digest,
    const blst_trx_public_key_material* public_key_material);

bool blst_trx_combine_shares(
    const blst_trx_partial_decryption_share* shares,
    size_t share_count,
    const blst_trx_batch_digest* batch_digest,
    const blst_trx_public_key_material* public_key_material,
    blst_p1* out_sigma);

bool blst_trx_decapsulate_one(
    const blst_trx_params* params,
    const blst_trx_batch_digest* batch_digest,
    size_t ciphertext_index,
    const blst_p1* proof,
    const blst_p1* sigma,
    const blst_trx_crs* crs,
    const blst_trx_public_key_material* public_key_material,
    byte* out_session_key,
    size_t out_session_key_len);

void blst_trx_batch_digest_free(blst_trx_batch_digest* batch_digest);

#ifdef __cplusplus
}
#endif

#endif
