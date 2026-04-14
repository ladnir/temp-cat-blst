#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include "blst/curve/point.h"
#include "blst/fields/Fr.h"

namespace bte::trx
{
    using Bytes = std::vector<std::uint8_t>;
    using RandomBytesSource = std::function<void(std::span<std::uint8_t>)>;

    struct TrxParameters
    {
        std::size_t max_batch_size = 0;
        std::size_t num_contexts = 0;
        std::size_t session_key_size = 32;
        std::uint8_t serialization_version = 1;

        std::string h1_dst = "BTE_TRX_H1_V1";
        std::string hf_dst = "BTE_TRX_HF_V1";
        std::string kdf_salt = "BTE_TRX_KDF_SALT_V1";
        std::string kdf_info = "BTE_TRX_KDF_INFO_V1";
        std::string ots_sign_dst = "BTE_TRX_OTS_BLS_SIG_V1";
    };

    struct CrsContext
    {
        std::vector<blst::G1Point> g1_powers;
    };

    struct Crs
    {
        std::size_t max_batch_size = 0;
        std::vector<CrsContext> contexts;
        blst::G2Point h_tau = blst::G2Point::infinity();
    };

    struct ValidatorSecretShare
    {
        std::uint32_t index = 0;
        blst::Fr sk_share{};
    };

    struct ValidatorPublicShare
    {
        std::uint32_t index = 0;
        blst::G2Point pk_share = blst::G2Point::infinity();
        blst::G2Point pk_tau_share = blst::G2Point::infinity();
    };

    struct PublicKeyMaterial
    {
        std::vector<ValidatorPublicShare> validator_public_shares;
        blst::G2Point pk = blst::G2Point::infinity();
        blst::G2Point pk_tau = blst::G2Point::infinity();
        std::size_t threshold = 0;
        mutable std::string cached_h1_dst;
        mutable bool has_cached_h1_pk = false;
        mutable blst::G1Point cached_h1_pk = blst::G1Point::infinity();
    };

    struct EncryptionKey
    {
        blst::G2Point pk = blst::G2Point::infinity();
        blst::G2Point pk_tau = blst::G2Point::infinity();
        mutable std::string cached_h1_dst;
        mutable bool has_cached_h1_pk = false;
        mutable blst::G1Point cached_h1_pk = blst::G1Point::infinity();
    };

    struct Ciphertext
    {
        std::array<std::uint8_t, 96> ct1{};
        std::array<std::uint8_t, 96> ct2{};
        Bytes ct3;
        std::array<std::uint8_t, 48> vk_sig{};
        std::array<std::uint8_t, 96> sig{};
    };

    struct CiphertextEnvelope
    {
        Ciphertext ciphertext;
        Bytes associated_data;
    };

    struct EncryptedCiphertext
    {
        Bytes session_key;
        Ciphertext ciphertext;
        Bytes associated_data;
    };

    struct VerifiedCiphertext
    {
        Ciphertext ciphertext;
        Bytes associated_data;
        blst::G2Point ct1 = blst::G2Point::infinity();
        blst::G2Point ct2 = blst::G2Point::infinity();
        blst::Fr tag{};
    };

    struct BatchDigest
    {
        std::size_t context_id = 0;
        std::vector<std::size_t> valid_indices;
        std::vector<VerifiedCiphertext> valid_ciphertexts;
        std::vector<blst::Fr> tags;
        std::vector<blst::Fr> polynomial;
        blst::G1Point commitment = blst::G1Point::infinity();
        mutable bool has_cached_eval_proofs = false;
        mutable std::vector<blst::G1Point> cached_eval_proofs;
    };

    struct PartialDecryptionShare
    {
        std::uint32_t index = 0;
        blst::G1Point value = blst::G1Point::infinity();
    };
}
